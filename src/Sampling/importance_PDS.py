#!/usr/bin/env python3

import os
import sys
from Sampling.sampler import Sampling_method
from random import randint
import numpy as np
import pandas as pd
import logging
from sympy import sympify
from pandas.core.series import Series


class Importance_Possion_Disk(Sampling_method):
    def __init__(self) -> None:
        super().__init__()

    def set_config(self, cf):
        return super().set_config(cf)

    def set_scan_path(self, pth):
        return super().set_scan_path(pth)

    def initialize_generator(self, cf):
        self.cf = cf
        self.init_logger(cf['logging']['scanner'])
        self.logger.info(
            "Possion Disk Sampling Algorithm is used in this scan")
        self.set_pars()
        self.runing_card = os.path.join(
            self.path['save dir'], 'possion_run.json')
        self.status = "INIT"
        from Func_lib import get_time_lock
        self.timelock = get_time_lock(
            self.cf['default setting']['sampling']['TSavePoint'])
        self.pars['maxTry'] = self.cf['default setting']['sampling']['maxTry']
        # print(self.path, self.runing_card)

    def set_pars(self):
        pars = {
            "variable": self.scf.get("Sampling_Setting", "variables"),
            "likelihood":   self.scf.get("Sampling_Setting", "likelihood")
        }
        self.pars = {"minR":     float(
            self.scf.get("Sampling_Setting", "min R"))}
        self.decode_sampling_variable_setting(pars['variable'])
        self.decode_function()
        self.decode_likelihood(pars['likelihood'])
        # print(self.pars.items())

    def decode_sampling_variable_setting(self, prs):
        self.pars['vars'] = {}
        nn = 0
        from math import log10
        for line in prs.split('\n'):
            ss = line.split(",")
            if len(ss) == 4 and ss[1].strip().lower() == "flat":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "flat",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "expr":     sympify("{} + ({} - {}) * cube{}".format(float(ss[2].strip()), float(ss[3].strip()), float(ss[2].strip()), nn))
                }
            elif len(ss) == 4 and ss[1].strip().lower() == "log":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "log",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "expr":     sympify("10**({} + ({} - {}) * cube{})".format(log10(float(ss[2].strip())), log10(float(ss[3].strip())), log10(float(ss[2].strip())), nn))
                }
                if not (float(ss[2].strip()) > 0 and float(ss[3].strip())):
                    self.logger.error(
                        "Illegal Variable setting of {}\n\t\t-> The lower limit and high limit should be > 0 ".format(ss))
                    sys.exit(0)
            else:
                self.logger.error(
                    "Illegal Variable setting in: {} ".format(line))
                sys.exit(0)
            nn += 1
        self.logger.warning(self.pars['vars'])

    def decode_likelihood(self, lik):
        if lik[0:4] == "CHI2":
            self.pars['likelihood'] = {
                "method":   "chisquare",
                "fcsinc":   [],
                "expression":   lik[5:].strip()
            }
        elif lik[0:4] == "LIKI":
            self.pars['likelihood'] = {
                "method":   "likelihood",
                "fcsinc":   [],
                "expression":   lik[5:].strip()
            }
        if "&FC" in self.pars['likelihood']['expression']:
            for item in self.exprs:
                if item in self.pars['likelihood']['expression']:
                    self.pars['likelihood']['expression'] = self.pars['likelihood']['expression'].replace(
                        item, self.exprs[item]['expr'])
            for item in self.fcs:
                if item in self.pars['likelihood']['expression']:
                    self.pars['likelihood']['fcsinc'].append(item)

    def decode_function(self):
        self.fcs = {}
        self.exprs = {}
        for sc in self.scf.sections():
            if "Function" == sc[0:8] and self.scf.get(sc, "method") == "expression":
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                self.exprs[name] = {
                    "expr":   self.scf.get(sc, "expression")
                }
            elif "Function" == sc[0:8] and self.scf.get(sc, "method") == "interpolate 1d":
                from pandas import read_csv
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                from Func_lib import decode_path_from_file
                data = read_csv(decode_path_from_file(
                    self.scf.get(sc, "file")))
                self.fcs[name] = {
                    "data":      data
                }
                if self.scf.has_option(sc, "kind"):
                    method = self.scf.get(sc, "kind").strip().lower()
                    if method not in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', 'next']:
                        method = "linear"
                else:
                    method = "linear"
                fill_value = ""
                if self.scf.has_option(sc, "fill_value"):
                    fill_value = self.scf.get(sc, "fill_value").strip().lower()
                from scipy.interpolate import interp1d
                if fill_value == "extrapolate":
                    self.fcs[name]['expr'] = interp1d(
                        data['x'], data['y'], kind=method, fill_value=fill_value)
                else:
                    self.fcs[name]['expr'] = interp1d(
                        data['x'], data['y'], kind=method)
        print(self.fcs)
        print(self.exprs)

    def init_logger(self, cf):
        self.logger = logging.getLogger("Possion")
        handler = {                
            "stream":   logging.StreamHandler(),
            "ff":       logging.FileHandler(cf['logging path'])
        }
        self.logger.setLevel(logging.DEBUG)
        handler['stream'].setLevel(logging.INFO)
        handler['ff'].setLevel(logging.DEBUG)
        self.logger.addHandler(handler['stream'])
        self.logger.addHandler(handler['ff'])
        from logging import Formatter
        handler['stream'].setFormatter(
            Formatter(cf['stream_format'],  "%m/%d %H:%M:%S"))
        handler['ff'].setFormatter(Formatter(cf['file_format']))

    def prerun_generator(self):
        self.get_empty_data_row()
        self.logger.info(
            "Start poisson disk sampling in {} dimensions ...".format(self.pars['ndim']))
        self.points = pd.DataFrame(columns=self.pars['emptyData'].keys())
        self.data = []
        self.cubes = pd.DataFrame(columns=self.pars['emptyCube'].keys())
        self.grays = pd.DataFrame(columns=self.pars['emptyCube'].keys())
        self.pars['cubeids'] = list(self.pars['emptyCube'].keys())[2:]
        self.path['Samples_info'] = os.path.join(
            self.path['save dir'], "Samples_info")
        self.path['cubedir'] = os.path.join(self.path['save dir'], "cubes.csv")
        self.path['graydir'] = os.path.join(self.path['save dir'], 'grays.csv')
        self.path['pointsdir'] = os.path.join(
            self.path['save dir'], "Points.csv")
        self.path['datadir'] = os.path.join(
            self.path['save dir'], "AllData.csv")
        if not os.path.exists(self.path['Samples_info']):
            os.makedirs(self.path['Samples_info'])
        if self.cf['default setting']['sampling']['demo']:
            # import matplotlib
            # matplotlib.use('Agg')
            from matplotlib import pyplot as plt
            plt.close()
            self.fig = plt.figure(figsize=(10, 10))
            self.ax = self.fig.add_axes([0., 0., 1., 1.])
            plt.ion()
            plt.show()
        self.status = "READY"

    def generate_events(self):
        while not self.status == "FINISH":
            self.check_samples_status_number(False)
            if self.status == "READY":
                cube = np.random.rand(self.pars['ndim'])
                self.add_points_by_cube(cube)
                self.check_samples_status_number(True)
                self.status = "RUNNING"
            elif self.status == "RUNNING":
                self.check_running_sample_and_run_next_step()
                if self.info['nsample']['live'] < self.info['nsample']['dead']:
                    self.check_dead_sample_is_gray()
                if self.info['nsample']['redy'] > 0 and self.info['nsample']['runn'] < self.pack['config']['paraller number']:
                    self.ready_sample_start_run()
                if self.info['nsample']['live'] > 0:
                    # if self.info['nsample']['live'] > 0 and self.info['nsample']['dead'] < self.pack['config']['paraller number']:
                    self.get_new_sample()
                self.check_generator_status()
                self.message_out_status()

            if self.cf['default setting']['sampling']['demo']:
                self.plot_demo()

    def check_running_sample_and_run_next_step(self):
        doneids = []
        stopids = []
        finiids = []
        for kk, smp in self.samples.items():
            if smp.status not in ['Done', 'Stopped']:
                stat = smp.status
                smp.update_status()
                if stat != smp.status:
                    self.change_point_status_in_data(smp.id, smp.status)
                    if smp.status == "Done":
                        doneids.append(smp.id)
                        self.data.append(smp.vars)
                    elif smp.status == "Stopped":
                        stopids.append(smp.id)
                        self.data.append(smp.vars)
                    if smp.status == "Finish":
                        finiids.append(smp.id)
        if doneids + stopids + finiids:
            self.check_samples_status_number(True)

    def ready_sample_start_run(self):
        sids = []
        for kk, smp in self.samples.items():
            if smp.status == "Ready":
                sids.append(smp.id)
        if sids:
            for sid in sids:
                self.samples[sid].start_run()
                self.change_point_status_in_data(sid, "Running")
            self.check_samples_status_number(True)

    def afterrun_generator(self):
        self.cubes.to_csv(self.path['cubedir'], index=False)
        self.grays.to_csv(self.path['graydir'], index=False)
        self.points.to_csv(self.path['pointsdir'], index=False)
        from copy import deepcopy
        sdata = deepcopy(self.data)
        sdata = pd.DataFrame(sdata)
        sdata.to_csv(self.path['datadir'], index=False)
        
        # self.data.to_csv(self.path['datadir'], index=False)
        
        
        if self.cf['default setting']['sampling']['demo']:
            self.ax.set_xlim(-0.5, 1.5)
            self.ax.set_ylim(-0.5, 1.5)

            self.fig.savefig("END.png", dpi=300)

    def check_generator_status(self):
        # if self.info['nsample']['live'] == 0 and self.info['nsample']['dead'] == 0 and self.info['nsample']['tot'] == self.info['nsample']['done'] + self.info['nsample']['stop']:
        if self.info['nsample']['live'] == 0 and self.info['nsample']['dead'] == 0:
            self.status = "FINISH"

        # plt.savefig("possiondisk.mp4", fps=30)
    def plot_demo(self):
        self.ax.cla()

        live = self.cubes[self.cubes['Status'] == "Live"]
        dead = self.cubes[self.cubes['Status'] == "Dead"]
        # gray = self.grays
        gray = self.cubes[self.cubes['Status'] == "Gray"]

        ss = pd.concat([self.cubes, self.grays])

        self.ax.scatter(live['cube0'], live['cube1'], s=7.5, color="tomato")
        self.ax.scatter(dead['cube0'], dead['cube1'],
                        s=7.5, color="dodgerblue")
        self.ax.scatter(gray['cube0'], gray['cube1'], s=17.5, color="grey")

        from scipy.spatial import Voronoi, voronoi_plot_2d
        ss = pd.DataFrame({"x": ss['cube0'], "y": ss['cube1']}).to_numpy()
        if ss.shape[0] > 4:
            vor = Voronoi(ss[:, :])
            voronoi_plot_2d(vor, self.ax, show_points=False,
                            show_vertices=False, line_colors=None, line_width=0.5)
        self.ax.set_xlim(0., 1.)
        self.ax.set_ylim(0., 1.)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def get_new_sample(self):
        if self.info['nsample']['live'] > 0:
            sid = self.pick_sample_ID_by_status("Live")
            # print(sid)
            smp = self.samples[sid]
            smp.local, smp.grayed = self.get_local(sid)
            smp.trys = self.sampling_sphere(sid)

            for ii in range(smp.trys.shape[0]):
                smp.children.append(self.add_points_by_cube(smp.trys[ii]))
            for cid in smp.children:
                self.samples[cid].mother = sid

            smp.status = "Ready"
            smp.trys = None
            smp.init_logger(self.cf['logging']['scanner'])
            self.change_point_status_in_cube(sid, "Dead")
            self.change_point_status_in_data(sid, "Ready")
            self.check_samples_status_number(False)
            # print(smp.__dict__)
        # self.message_out_status()

    def check_dead_sample_is_gray(self):
        sid = self.pick_sample_ID_by_status("Dead")
        smp = self.samples[sid]
        if self.info['nsample']['live'] == 0:
            self.change_point_status_in_cube(sid, "Gray")
        else:
            smp.local, smp.grayed = self.get_local(sid, _max=3.1, _min=0.1)
            if smp.grayed:
                self.change_point_status_in_cube(sid, "Gray")
                self.check_samples_status_number(False)

    def get_local(self, sid, _max=2.1, _min=0.1):
        from copy import deepcopy
        cbs = np.array(deepcopy(self.cubes[self.pars['cubeids']]))
        # cds = np.where((np.linalg.norm(cbs - self.samples[sid].cube, axis=1) < 5.0 * self.pars['minR']) & (np.linalg.norm(cbs - self.samples[sid].cube, axis=1) > 0.1 * self.pars['minR']))
        cds = np.where((np.linalg.norm(cbs - self.samples[sid].cube, axis=1) < _max * self.pars['minR']) & (
            np.linalg.norm(cbs - self.samples[sid].cube, axis=1) > _min * self.pars['minR']))
        sts = np.array(deepcopy(self.cubes['Status']))[cds]
        grayed = False
        if np.count_nonzero(sts == "Live") == 0:
            grayed = True
        return cbs[cds], grayed

    def change_point_status_in_data(self, sid, status):
        idx = self.points[self.points['ID'] == sid].index
        self.points.at[idx, "Status"] = status

    def change_point_status_in_cube(self, sid, status):
        idx = self.cubes[self.cubes['ID'] == sid].index
        self.cubes.at[idx, "Status"] = status

        # if status != "Gray":
        #     self.cubes.at[idx, "Status"] = status
        # else:
        #     new = self.cubes.loc[self.cubes['ID'] == sid].iloc[0]
        #     new['Status'] = status
        #     # self.cubes = self.cubes.drop(self.cubes[self.cubes['ID'] == sid].index, inplace=True)
        #     self.cubes.drop(self.cubes[self.cubes['ID'] == sid].index, inplace=True)
        #     self.grays = self.grays.append(new, ignore_index=True)

    def add_points_by_cube(self, cube):
        from copy import deepcopy
        from sample import Sample
        new = Sample()
        new.status = "Free"
        new.grayed = False

        cub = deepcopy(self.pars['emptyCube'])
        from Func_lib import get_sample_id
        new.id = get_sample_id()
        cub['ID'] = new.id
        for ii in range(self.pars['ndim']):
            cub['cube{}'.format(ii)] = cube[ii]
        self.cubes = self.cubes.append(cub, ignore_index=True)

        raw = deepcopy(self.pars['emptyData'])
        raw['ID'] = cub['ID']
        for kk, vv in self.pars['vars'].items():
            raw[kk] = vv['expr'].subs(cub)
        self.points = self.points.append(raw, ignore_index=True)

        raw = dict(raw)
        raw.pop("ID")

        new.par = raw
        new.cube = cube
        new.local = []
        new.pack = self.pack
        new.path['info'] = os.path.join(
            self.path['Samples_info'], str(cub['ID']))
        new.path['scanner_run_info'] = self.path['run_info']

        self.samples[new.id] = new
        if not os.path.exists(self.samples[new.id].path['info']):
            os.makedirs(self.samples[new.id].path['info'])
        else:
            from shutil import rmtree
            rmtree(self.samples[new.id].path['info'])
            os.makedirs(self.samples[new.id].path['info'])
            self.samples[new.id].init_logger(self.cf['logging']['scanner'])
        self.check_samples_status_number()
        return new.id

    def pick_sample_ID_by_status(self, status):
        lives = self.cubes.loc[self.cubes['Status'] == status]
        # smp = lives.iloc[randint(0, lives.shape[0] -1)]
        if lives.shape[0] > 10:
            idx = randint(0, min(10, lives.shape[0] - 1))
            smp = lives.iloc[idx]
        else:
            smp = lives.iloc[0]
        return smp['ID']

    def sampling_volume(self, sid):
        vecs = np.random.standard_normal(
            size=(self.pars['maxTry'], self.pars['ndim']))
        vecs /= np.linalg.norm(vecs, axis=1)[:, None]
        scal = 1.0 + np.random.random(self.pars['maxTry'])
        vecs *= scal[:, None]
        vecs = self.samples[sid].cube + self.pars['minR'] * vecs
        vecs = vecs[np.where(np.min(vecs, axis=1) > 0.)]
        vecs = vecs[np.where(np.max(vecs, axis=1) < 1.)]

        from scipy.spatial.distance import cdist
        if self.samples[sid].local.size != 0:
            cds = cdist(vecs, self.samples[sid].local)
            vecs = np.delete(vecs, np.where(
                np.min(cds, axis=1) < self.pars['minR']), axis=0)

        if vecs.size != 0:
            cd = cdist(vecs, vecs) + \
                np.tril(np.ones((vecs.shape[0], vecs.shape[0])))
            vecs = np.delete(vecs, np.where(
                np.min(cd, axis=0) < self.pars['minR']), axis=0)

        return vecs

    def sampling_sphere(self, sid):
        vecs = np.random.standard_normal(
            size=(self.pars['maxTry'], self.pars['ndim']))
        vecs /= np.linalg.norm(vecs, axis=1)[:, None]
        vecs = self.samples[sid].cube + self.pars['minR'] * vecs
        vecs = vecs[np.where(np.min(vecs, axis=1) > 0.)]
        vecs = vecs[np.where(np.max(vecs, axis=1) < 1.)]

        from scipy.spatial.distance import cdist
        if self.samples[sid].local.size != 0:
            cds = cdist(vecs, self.samples[sid].local)
            vecs = np.delete(vecs, np.where(
                np.min(cds, axis=1) < self.pars['minR']), axis=0)

        if vecs.size != 0:
            cd = cdist(vecs, vecs) + \
                np.tril(np.ones((vecs.shape[0], vecs.shape[0])))
            vecs = np.delete(vecs, np.where(
                np.min(cd, axis=0) < self.pars['minR']), axis=0)

        return vecs

    def check_samples_status_number(self, output=False):
        self.info['nsample'] = {
            "tot":      self.points.shape[0],
            "live":     self.cubes.loc[self.cubes['Status'] == "Live"].shape[0],
            "dead":     self.cubes.loc[self.cubes['Status'] == "Dead"].shape[0],
            "gray":     self.cubes.loc[self.cubes['Status'] == "Gray"].shape[0],
            "free":     self.points.loc[self.points['Status'] == "Free"].shape[0],
            "redy":     self.points.loc[self.points['Status'] == "Ready"].shape[0],
            "runn":     self.points.loc[self.points['Status'] == "Running"].shape[0],
            "fini":     self.points.loc[self.points['Status'] == "Finish"].shape[0],
            "done":     self.points.loc[self.points['Status'] == "Done"].shape[0],
            "stop":     self.points.loc[self.points['Status'] == "Stopped"].shape[0]
        }
        if output:
            self.message_out_status()

    def message_out_status(self):
        self.logger.info(
            "{} samples in total: {} live, {} dead, {} gray:  \n\t   => {} Free, {} Ready, {} running, {} Finished, {} Done, {} Stopped".format(
                self.info['nsample']['tot'],
                self.info['nsample']['live'],
                self.info['nsample']['dead'],
                self.info['nsample']['gray'],
                self.info['nsample']['free'],
                self.info['nsample']['redy'],
                self.info['nsample']['runn'],
                self.info['nsample']['fini'],
                self.info['nsample']['done'],
                self.info['nsample']['stop']
            ))

    def get_empty_data_row(self):
        raw = {
            "ID":   None,
            "Status":   "Free"
        }
        cub = {
            "ID":   None,
            "Status":   "Live"
        }
        self.pars['ndim'] = 0
        ii = 0
        for var, item in self.pars['vars'].items():
            raw[var] = None
            self.pars['ndim'] += 1
            cub["cube{}".format(ii)] = None
            ii += 1
        self.pars['emptyData'] = Series(raw)
        self.pars['emptyCube'] = cub
