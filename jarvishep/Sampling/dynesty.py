#!/usr/bin/env python3

import json
import os
import threading
from uuid import uuid4

import numpy as np
import pandas as pd

from jarvishep.Sampling.Source.Dynesty.py.dynesty.pool import JarvisFactoryAsyncPool
from jarvishep.Sampling.sampler import SamplingVirtial
from jarvishep.sample import Sample

def prior_transform(u):
    uuid = str(uuid4())
    ret = np.append(u, [uuid])
    return ret


DynestyFactoryPool = JarvisFactoryAsyncPool


class Dynesty(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Dynesty"
        self.sampler = None
        self.max_workers = os.cpu_count() or 1
        self._dynesty_pool = None
        self._dynesty_workers = 1
        self._factory_submit_limit = 1
        self._submit_gate = threading.BoundedSemaphore(value=1)
        self._execution_profile = {}

    def load_schema_file(self):
        self.schema = self.path['DynestySchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.set_bucket_alloc()
        self.init_generator()

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Dynesty sampler")

    def set_execution_profile(self, dynesty_workers=None, max_pending_factory=None):
        """Optional runtime override for dynesty/factory coordination."""
        profile = {}
        if dynesty_workers is not None:
            profile["dynesty_workers"] = int(dynesty_workers)
        if max_pending_factory is not None:
            profile["max_pending_factory"] = int(max_pending_factory)
        self._execution_profile = profile

    def __next__(self):
        # for dynesty sampling method, next() method is only for testing the assembing line, not the real sampling method
        u = np.random.random(self._dimensions)
        param = self.map_point_into_distribution(u)
        return param
        

    def map_point_into_distribution(self, row) -> np.ndarray:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii])
        return result

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        self.load_variable()
        # print(self.vars[0]._parameters)
        self._nlive = self.config['Sampling']['Bounds']['nlive']
        self._rstate = np.random.default_rng(self.config['Sampling']['Bounds']['rseed'])
        self._runnested = self.config['Sampling']['Bounds']['run_nested']
        self._dimensions = len(self.vars)

    def initialize(self):
        self.logger.warning("Initializing the Dynesty Sampling")

    def init_sampler_db(self):
        self.info['db'] = {
            "nested result":    os.path.join(self.info['sample']['task_result_dir'], "DATABASE", "dynesty_result.csv")
        }

    def _shutdown_dynesty_pool(self):
        if self._dynesty_pool is not None:
            self._dynesty_pool.shutdown(wait_for_tasks=True, cancel_futures=True)
            self._dynesty_pool = None

    def _ensure_bucket_allocator(self):
        if getattr(self, "bucket_alloc", None) is not None:
            return

        sample_cfg = self.info.get("sample", {}) if isinstance(self.info, dict) else {}
        base_path = sample_cfg.get("sample_dirs")
        if not base_path:
            task_root = sample_cfg.get("task_result_dir", os.getcwd())
            base_path = os.path.join(task_root, "SAMPLE")
            sample_cfg["sample_dirs"] = base_path

        limit = 200
        width = 6
        cfg = getattr(self, "config", None)
        if isinstance(cfg, dict):
            scan_cfg = cfg.get("Scan", {})
            directory_cfg = None
            if isinstance(scan_cfg, dict):
                directory_cfg = scan_cfg.get("sample_directory", None)
            if not isinstance(directory_cfg, dict):
                directory_cfg = cfg.get("Directory_Setting", None)
            if isinstance(directory_cfg, dict):
                limit = int(directory_cfg.get("limit", limit))
                width = int(directory_cfg.get("width", width))

        from jarvishep.Sampling.bucketallocator import BucketAllocator

        self.bucket_alloc = BucketAllocator(
            base_path=base_path,
            limit=limit,
            width=width,
            start_bucket=1,
            on_bucket_sealed=self._on_bucket_sealed if self._archive_enabled() else None,
        )
        self.bucket_alloc.check_and_update()
        if self._archive_enabled():
            self._ensure_archive_manager()

    def _resolve_execution_profile(self):
        factory_workers = int(getattr(self.factory, "_max_workers", self.max_workers) or self.max_workers or 1)
        nlive = int(getattr(self, "_nlive", 1) or 1)
        base_workers = int(self.max_workers or 1)

        requested_dynesty_workers = int(self._execution_profile.get("dynesty_workers", base_workers))
        dynesty_workers = max(1, min(requested_dynesty_workers, base_workers, factory_workers, nlive))

        requested_pending = int(self._execution_profile.get("max_pending_factory", dynesty_workers))
        max_pending_factory = max(1, min(requested_pending, factory_workers))

        self._dynesty_workers = dynesty_workers
        self._factory_submit_limit = max_pending_factory
        self._submit_gate = threading.BoundedSemaphore(value=max_pending_factory)

        if self._dynesty_pool is None or self._dynesty_pool.size != dynesty_workers:
            self._shutdown_dynesty_pool()
            self._dynesty_pool = DynestyFactoryPool(dynesty_workers)

        self._log_execution_profile(
            dynesty_workers=dynesty_workers,
            max_pending_factory=max_pending_factory,
            factory_workers=factory_workers,
            nlive=nlive,
        )

    @staticmethod
    def _format_worker_summary(rows):
        metric_header = "Metric"
        value_header = "Value"
        metric_width = max([len(metric_header)] + [len(str(k)) for k, _ in rows])
        value_width = max([len(value_header)] + [len(str(v)) for _, v in rows])

        lines = [
            f"{metric_header:<{metric_width}}\t{value_header:<{value_width}}",
            f"{'-' * metric_width}\t{'-' * value_width}",
        ]
        for metric, value in rows:
            lines.append(f"{str(metric):<{metric_width}}\t{str(value):<{value_width}}")

        return "\n\t" + "\n\t".join(lines)

    def _log_execution_profile(self, dynesty_workers, max_pending_factory, factory_workers, nlive):
        table = self._format_worker_summary(
            [
                ("dynesty_workers", dynesty_workers),
                ("factory_submit_limit", max_pending_factory),
                ("factory_workers", factory_workers),
                ("nlive", nlive),
            ]
        )
        self.logger.warning(f"Dynesty Worker Summary ->{table}")

    def _submit_with_backpressure(self, sample):
        if not self._submit_gate.acquire(blocking=False):
            self.logger.info(
                "Dynesty factory submit window full -> wait -> limit={} | uuid={}".format(
                    self._factory_submit_limit,
                    sample.uuid,
                )
            )
            self._submit_gate.acquire()
        try:
            future = self.factory.submit_task(sample.info)
            return future.result()
        finally:
            self._submit_gate.release()

    def run_nested(self):
        self._ensure_bucket_allocator()
        self._resolve_execution_profile()
        base_sample_cfg = self.info['sample']

        def log_likelihood(params):
            param = params[0:-1].astype(np.float64, copy=False)
            uuid = params[-1]
            pars = self.map_point_into_distribution(param)
            sample = Sample(pars)
            sample.update_uuid(uuid)
            sample_cfg = self.build_sample_config(
                base_sample_cfg,
                save_dir=self.bucket_alloc.next_bucket_dir(),
            )
            sample.set_config(sample_cfg)
            try:
                return self._submit_with_backpressure(sample)
            except Exception as exc:
                self.logger.error(f"[WorkerFactory] future exception consumed: uuid={sample.uuid} error={exc}")
                raise
            finally:
                sample.close()

        
        self.init_sampler_db()

        from jarvishep.Sampling.Source.Dynesty.py.dynesty import DynamicNestedSampler
        try:
            self.sampler = DynamicNestedSampler(
                loglikelihood=log_likelihood, 
                prior_transform=prior_transform,
                ndim=self._dimensions,
                nlive=self._nlive,
                pool=self._dynesty_pool,
                rstate=self._rstate,
                queue_size=self._dynesty_workers,
                log_file_path=self.info['logfile']
            )
            self.sampler.run_nested(**self._runnested)
        finally:
            self._shutdown_dynesty_pool()

    def finalize(self):
        self.save_dynesty_results_to_csv()
        # self.plot_dynesty_results()

    def _build_logl_to_logvol_interpolator(self):
        from scipy.interpolate import interp1d

        results = self.sampler.results
        try:
            x_raw = results["logl"]
        except Exception:
            x_raw = getattr(results, "logl", [])
        try:
            y_raw = results["logvol"]
        except Exception:
            y_raw = getattr(results, "logvol", [])

        x = np.asarray(x_raw, dtype=float)
        y = np.asarray(y_raw, dtype=float)
        finite = np.isfinite(x) & np.isfinite(y)
        x = x[finite]
        y = y[finite]
        if x.size < 2:
            self.logger.warning("Dynesty logl/logvol interpolation disabled -> insufficient finite points")
            return None

        order = np.argsort(x, kind="mergesort")
        x_sorted = x[order]
        y_sorted = y[order]

        unique_x = []
        unique_y = []
        for xv, yv in zip(x_sorted, y_sorted):
            if unique_x and xv == unique_x[-1]:
                unique_y[-1] = yv
            else:
                unique_x.append(xv)
                unique_y.append(yv)

        x_unique = np.asarray(unique_x, dtype=float)
        y_unique = np.asarray(unique_y, dtype=float)
        if x_unique.size < 2:
            self.logger.warning("Dynesty logl/logvol interpolation disabled -> no unique logl points")
            return None

        kind = "cubic" if x_unique.size >= 4 else "linear"
        return interp1d(
            x_unique,
            y_unique,
            kind=kind,
            bounds_error=False,
            fill_value="extrapolate",
            assume_sorted=True,
        )

    def save_dynesty_results_to_csv(self):
        data = {
            "uuid":             self.sampler.results['samples_uid'],
            "log_weight":       self.sampler.results['logwt'],
            "log_Like":         self.sampler.results['logl'],
            "log_PriorVolume":  self.sampler.results['logvol'],
            "log_Evidence":     self.sampler.results['logz'],
            "log_Evidence_err": self.sampler.results['logzerr'],
            "samples_nlive":    self.sampler.results['samples_n'],
            "ncall":            self.sampler.results['ncall'],
            "samples_it":       self.sampler.results['samples_it'],
            "samples_id":       self.sampler.results['samples_id'],
            "information":      self.sampler.results['information']
        }

        for ii in range(self.sampler.results['samples'].shape[1]):
            data[f'samples_v[{ii}]'] = self.sampler.results['samples'][:, ii]
        for ii in range(self.sampler.results['samples_u'].shape[1]):
            data[f'samples_u[{ii}]'] = self.sampler.results['samples_u'][:, ii]

        self.df = pd.DataFrame(data)
        self.df.to_csv(self.info['db']['nested result'], index=False)
        self.logger.warning(f"Results saved to {self.info['db']['nested result']}")
        summary_text = self.sampler.results.summary()
        if summary_text:
            self.logger.warning(f"Dynesty summary -> {summary_text}")
        self.lnX_from_LogLike = self._build_logl_to_logvol_interpolator()



    def plot_dynesty_results(self) -> None: 
        savepath = os.path.join(self.info['sample']['task_result_dir'], "dynesty_summary.png")
        import matplotlib.pyplot as plt
        maxid = self.df.shape
        self.logger.info(f"Dynesty plot: total rows={maxid[0]}")
        nlive = self.df.iloc[0]['samples_nlive']
        fig = plt.figure(figsize=(10, 11))
        ax = fig.add_axes([0.01, 0.99-0.5/11., 0.05, 0.5/11.])
        from jarvishep.plot import draw_logo_in_square
        draw_logo_in_square(ax)
        data = [
            {
                "y": self.df['samples_nlive'],
                "label":    "$\\text{Live points}$"
            },
            {   
                "y": np.exp(self.df['log_Like'] - max(self.df['log_Like'])),
                "label":    "Likelihood"
            },
            {
                "y": np.exp(self.df['log_weight']),
                "label":    "Importance\nweight PDF"
            },
            {
                "y": np.exp(self.df['log_Evidence']),
                "label":    "Evidence"
            },
            {
                "y": self.df['samples_it'],
                "label":    "Iters" 
            }
        ]
        # print()

        from matplotlib.ticker import AutoMinorLocator
        ax1 = fig.add_axes([0.15, 0.8/11., 0.83, 2/11.])
        ax1.scatter( - self.df['log_PriorVolume'], data[0]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax1.set_xlim(0., max(- self.df['log_PriorVolume']))
        ax1.set_ylim(0., max(data[0]['y']) * 1.1)
        ax1.set_ylabel(data[0]['label'], fontsize=18)
        ax1.yaxis.set_label_coords(-0.08, 0.5)
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax1.tick_params(which='major', length=7)
        ax1.tick_params(which='minor', length=4)
        # ax1.set_xticklabels([])
        plt.draw()


        ax2 = fig.add_axes([0.15, 2.8/11., 0.83, 2/11.])
        ax2.scatter( - self.df['log_PriorVolume'], data[1]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax2.set_xlim(0., max(- self.df['log_PriorVolume']))
        ax2.set_ylim(0., max(data[1]['y']) * 1.1)
        ax2.set_ylabel(data[1]['label'], fontsize=18)
        ax2.yaxis.set_label_coords(-0.08, 0.5)
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax2.tick_params(which='major', length=7)
        ax2.tick_params(which='minor', length=4)
        ax2.set_xticklabels([])


        plt.draw()

        from scipy.stats import gaussian_kde
        from scipy.interpolate import interp1d
        ax3 = fig.add_axes([0.15, 4.8/11., 0.83, 2/11.])
        from copy import deepcopy
        dtt = np.array(deepcopy(data[2]['y']))
        endid = np.where(dtt > dtt[-1])[-1][-1]
        self.logger.info(f"Dynesty plot: highlighted index={endid}")

        dtt = dtt / max(dtt)
        ax3.scatter( - self.df['log_PriorVolume'], dtt, marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        wt_kde = gaussian_kde(np.array(- self.df['log_PriorVolume']), weights=np.array(data[2]['y']), bw_method="silverman")
        logvol = np.linspace(min(- self.df['log_PriorVolume']), max(- self.df['log_PriorVolume']), 1000)
        wt = wt_kde(logvol)
        wt = wt / wt.max()
        ax3.plot(logvol, wt, '-', linewidth=1.8, color='#8bc34a', alpha=0.8)
        ax3.set_xlim(0., max(- self.df['log_PriorVolume']))
        ax3.set_ylim(0., 1.1)
        ax3.set_ylabel(data[2]['label'], fontsize=18)
        ax3.yaxis.set_label_coords(-0.08, 0.5)
        ax3.yaxis.set_minor_locator(AutoMinorLocator())
        ax3.xaxis.set_minor_locator(AutoMinorLocator())
        ax3.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax3.tick_params(which='major', length=7)
        ax3.tick_params(which='minor', length=4)
        ax3.set_xticklabels([])
        ax3.text(0.02, 0.96, "Normalized to the maximum value", ha='left', va='top', transform=ax3.transAxes)

        plt.draw()

        ax4 = fig.add_axes([0.15, 6.8/11., 0.83, 2/11.])
        ax4.scatter( - self.df['log_PriorVolume'], data[3]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax4.set_xlim(0., max(- self.df['log_PriorVolume']))
        ax4.set_ylim(0., max(data[3]['y']) * 1.4)
        ax4.set_ylabel(data[3]['label'], fontsize=18)
        ax4.yaxis.set_label_coords(-0.08, 0.5)
        ax4.yaxis.set_minor_locator(AutoMinorLocator())
        ax4.xaxis.set_minor_locator(AutoMinorLocator())
        ax4.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax4.tick_params(which='major', length=7)
        ax4.tick_params(which='minor', length=4)
        ax4.set_xticklabels([])
        plt.draw()
        offset_text = ax4.yaxis.get_offset_text().get_text()
        maxx = str(max(data[3]['y'])).split("e")[0][0:5]
        upp = offset_text.split("e")[-1]
        txt = r"$\times 10^{" + upp + r"}$"
        ax4.text(0.4, max(data[3]['y']) * 0.95, maxx+txt, ha='left', va='top')
        ax4.plot([0, max(-self.df['log_PriorVolume'])], [max(data[3]['y']), max(data[3]['y'])], '-', linewidth=0.8, color="grey", alpha=0.4)

        logzerr = np.array(self.df['log_Evidence_err'])  # error in ln(evidence)
        logzerr[~np.isfinite(logzerr)] = 0.
        ax4.fill_between(- self.df['log_PriorVolume'], np.exp(self.df['log_Evidence'] - logzerr  ), np.exp(self.df['log_Evidence'] + logzerr), color='#8bc34a', alpha=0.2)

        ax5 = fig.add_axes([0.15, 8.8/11., 0.83, 2/11.])
        ax5.scatter( - self.df['log_PriorVolume'], data[4]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax5.set_xlim(0., max(- self.df['log_PriorVolume']))
        ax5.set_ylim(0., max(data[4]['y']) * 1.1)
        ax5.set_ylabel(data[4]['label'], fontsize=18)
        ax5.yaxis.set_label_coords(-0.08, 0.5)
        ax5.yaxis.set_minor_locator(AutoMinorLocator())
        ax5.xaxis.set_minor_locator(AutoMinorLocator())
        ax5.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax5.tick_params(which='major', length=7)
        ax5.tick_params(which='minor', length=4)
        ax5.set_xticklabels([])

        ax1.axvline( -self.df.iloc[endid]['log_PriorVolume'], color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        ax2.axvline( -self.df.iloc[endid]['log_PriorVolume'], color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        ax3.axvline( -self.df.iloc[endid]['log_PriorVolume'], color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        ax4.axvline( -self.df.iloc[endid]['log_PriorVolume'], color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        ax5.axvline( -self.df.iloc[endid]['log_PriorVolume'], color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)

        plt.draw()

        ax1.set_xlabel(r"$-\ln(X)$", fontsize=24)
        plt.savefig(savepath, dpi=300)


    def _resolve_full_dataframe_path(self, fulldf):
        candidates = []
        if isinstance(fulldf, str) and fulldf:
            candidates.append(fulldf)
            root, ext = os.path.splitext(fulldf)
            if ext.lower() == ".hdf5":
                candidates.append(f"{root}.0.csv")
                candidates.append(f"{root}.csv")

        task_dir = self.info.get("sample", {}).get("task_result_dir")
        if isinstance(task_dir, str) and task_dir:
            db_dir = os.path.join(task_dir, "DATABASE")
            db_info = os.path.join(db_dir, "running.json")
            if os.path.exists(db_info):
                try:
                    with open(db_info, "r", encoding="utf-8") as f1:
                        meta = json.load(f1)
                    converted = meta.get("converted", [])
                    if isinstance(converted, str):
                        converted = [converted]
                    if isinstance(converted, list):
                        # Prefer latest converted csv first.
                        candidates.extend(list(reversed([p for p in converted if isinstance(p, str) and p])))
                    pathroot = meta.get("pathroot")
                    active_no = meta.get("activeNO")
                    if isinstance(pathroot, str) and isinstance(active_no, int):
                        candidates.append(f"{pathroot}.{active_no}.csv")
                except Exception:
                    pass

        seen = set()
        for path in candidates:
            if not path or path in seen:
                continue
            seen.add(path)
            if os.path.exists(path):
                return path
        return None

    def combine_data(self, fulldf):
        if self.df is None:
            return
        full_df_path = self._resolve_full_dataframe_path(fulldf)
        if not full_df_path:
            self.logger.warning(
                "Dynesty combine_data skipped -> no full sample table found (input={})".format(fulldf)
            )
            return
        df_full = pd.read_csv(full_df_path)
        merged_df = pd.merge(df_full, self.df, on="uuid", how='inner')
        merged_path = os.path.join(self.info['sample']['task_result_dir'], "DATABASE", "dynesty_full.csv")
        os.makedirs(os.path.dirname(merged_path), exist_ok=True)
        merged_df.to_csv(merged_path, index=False)
        if self.lnX_from_LogLike is not None and "LogL" in df_full.columns:
            df_full['log_PriorVolume'] = self.lnX_from_LogLike(df_full['LogL'])
        else:
            df_full['log_PriorVolume'] = np.nan
