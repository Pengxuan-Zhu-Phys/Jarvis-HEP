#!/usr/bin/env python3 

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import os, sys 
from _pytest import config
import json 
import emoji
import math 
from copy import deepcopy
from matplotlib.colors import LogNorm
pwd = os.path.abspath(os.path.dirname(__file__))
jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
from plot import BudingPLOT
import update_colors
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, FixedLocator, NullLocator, AutoMinorLocator)
import time

from inner_func import update_funcs

class Figure(BudingPLOT):
    def __init__(self, config, name, sts):
        self.name = name 
        self.type = sts['type']
        self.yaml = config
        self.data = None
        self.info = sts
        self.start = time.time()
        self.path = os.path.join(self.yaml["Plot_Config"]['save_dir'], self.name)
        self.style = sts.get("style", "2col")
        
    def savefig(self, fig, plt):
        return super().savefig(fig, plt)
    
    def load_selection(self):
        if self.info.get("selection", ""):
            self.data = self.load_bool_df(self.data, self.info['selection'])
            if self.data.shape[0] < 1:
                print(emoji.emojize('\t:ghost::ghost::ghost: No data selected!!\n\tPlease check your Data or the selection condition \n', language='alias'))
            else: 
                print(emoji.emojize("\t:space_invader::space_invader::space_invader: Selected Data is -> {} rows".format(self.data.shape[0]), language='alias'))

    def load_bool_df(self, df, condition):
        # try:
        #     filtered_df = df.query(condition, engine="python")
        #     return filtered_df
        # except Exception as e:
        #     print("Errors happens when load data in condition ->{}".format(condition), e)
        #     return pd.DataFrame()
        try:
            allowed_globals = update_funcs({"np": np, "math": math})
            local_vars = df.to_dict("series")
            mask = eval(condition, allowed_globals, local_vars)
            if not isinstance(mask, pd.Series):
                mask = pd.Series(mask, index=df.index)
            return df[mask]
        except Exception as e:
            print("Errors when evaluating condition -> {}:\n\t".format(condition), e)
            return pd.DataFrame()

    def load_vars(self, expr, df=None):
        if df is None: 
            df = self.data
        allowed_globals = update_funcs({"np": np, "math": math})
        local_vars = df.to_dict("series")
        try:
            result = eval(expr, allowed_globals, local_vars)
            return result
        except Exception as e:
            print("Errors when evaluate expression -> {}\n\t:".format(expr), e)
            return pd.Series(dtype=float)

    def load_lim(self, lim):
        allowed_globals = update_funcs({"np": np, "math": math})
        res = []
        for expr in lim: 
            if isinstance(expr, str):
                xx = eval(expr, allowed_globals)
            else: 
                xx = expr 
            res.append(xx)
        return res 
    
    def set_ticks_format(self, ax):
        ax.tick_params(**self.para["ticks"])
        ax.tick_params(**self.para["ticks_major"])
        ax.tick_params(**self.para["ticks_minor"])
    
    def savefig(self, fig):
        from matplotlib.backends.backend_pdf import PdfPages
        support_fmt_list = ['ps', 'eps', 'pdf', 'pgf', 'png', 'raw',
                            'rgba', 'svg', 'svgz', 'jpg', 'jpeg', 'tif', 'tiff']
        if "save_format" not in self.yaml["Plot_Config"]:
            self.yaml["Plot_Config"]['save_format'] = ['pdf']
            
        unsupport = []
        support = []
        file_list = []
        for fmt in self.yaml["Plot_Config"]['save_format']:
            if fmt.strip().lower() in support_fmt_list:
                support.append(fmt.strip().lower())
                file_list.append("'{}.{}'".format(
                    self.name, fmt.strip().lower()))
            else:
                unsupport.append("'*.{}'".format(fmt.strip().lower()))    
                
        if len(support) == 0:
            support = ['pdf']
            file_list = ["'{}.{}'".format(self.name, 'pdf')]
            print("\tTimer: {:.2f} Second;  Message from '{}' -> No support picture format found in configure file! Default format '*.pdf' used in this plot".format(
                time.time()-self.start, self.name))
            
        for fmt in support: 
            if (fmt == 'ps'):
                fig.savefig("{}.pdf".format(self.path), format='pdf')
                self.compress_figure_to_PS(self.path)
            elif fmt == "pdf" and ('ps' not in self.yaml["Plot_Config"]['save_format']):
                fig.savefig("{}.{}".format(self.path, fmt))
            else:
                fig.savefig("{}.{}".format(self.path, fmt), dpi=150)

        if ('pdf' not in support) and ('ps' in support):
            os.remove("{}.pdf".format(self.path))

        print("\nFigure {} saved in the path\n\t\t-> {} \n\t\t>> {}.".format(
            self.name, os.path.dirname(self.path), ", >> ".join(file_list)))
        if unsupport:
            print(emoji.emojize('\t:ghost::ghost::ghost: Figure format unsupport -> {}. '.format(
                ", ".join(unsupport)), use_aliases=True))

    def compress_figure_to_PS(self, figpath):
        os.system('pdf2ps {}.pdf {}.ps'.format(figpath, figpath))

    def draw_text_labels(self, ax):
        """
        Draw text labels based on 'text' entries in self.info.
        Each text entry should be a dict with:
          - 'text': the string to draw
          - 'x', 'y': coordinates in data or axes space
          - optional 'transform': 'axes' to use axes coordinates, otherwise data coords
          - optional 'style': dict of text properties for ax.text()
        """
        for item in self.info.get("text", []):
            txt = item.get("text", "")
            x = item.get("x", 0)
            y = item.get("y", 0)
            style = item.get("style", {}).copy()
            transform_type = item.get("transform")
            if transform_type == "axes":
                ax.text(x, y, txt, transform=ax.transAxes, **style)
            else:
                ax.text(x, y, txt, **style)

    def draw_lines(self, ax):
        """
        Draw custom lines based on 'lines' entries in self.info.
        Each entry should be a dict with:
          - 'x': list of x coordinates
          - 'y': list of y coordinates
          - optional 'transform': 'axes' to use axes coords
          - optional 'style': dict of line properties for ax.plot()
        """
        for item in self.info.get("lines", []):
            xs = item.get("x", [])
            ys = item.get("y", [])
            style = item.get("style", {}).copy()
            if item.get("transform") == "axes":
                ax.plot(xs, ys, transform=ax.transAxes, **style)
            else:
                ax.plot(xs, ys, **style)
    
class Scatter(Figure):
    def __init__(self, config, name, sts):
        super().__init__(config=config, name=name, sts=sts) 
        with open(os.path.join(pwd, "cards/scatter.json"), 'r') as f1: 
            setting = json.loads(f1.read())
            self.para = setting[self.style]
    
    def serialize(self, fig, ax):
        if self.info.get("serialize", False):
            import pickle
            file = {"fig": fig, "ax": ax}
            pickle_file = "{}.pkl".format(self.path)
            with open(pickle_file, 'wb') as f:
                pickle.dump(file, f)
            print(f"Figure saved to {pickle_file}")
    
    def draw(self):
        plt.close("all")
        fig = plt.figure(**self.para["figure"])
        ax  = fig.add_axes(**self.para["ax"])
        for spine in ax.spines.values():
            spine.set_linewidth(self.para["spines"]['linewidth'])       
        self.load_selection()
        xx = self.info["parameters"][0]
        yy = self.info["parameters"][1]
        if self.info.get("group", ""):
            nn = 0 
            from itertools import cycle
            st = cycle(deepcopy(self.para['scatter']))
            for item in self.info['group']: 
                dfs = self.load_bool_df(self.data, item['condition'])
                xi = self.load_vars(xx['expr'], dfs)
                yi = self.load_vars(yy['expr'], dfs)
                style = next(st)
                style.update(item.get("style", {}))
                ax.scatter(xi, yi, **style)
                nn += 1
        else: 
            xx['data'] = self.load_vars(xx['expr'])   
            yy['data'] = self.load_vars(yy['expr'])
            ax.scatter(xx['data'], yy['data'], **self.para["scatter"][0])
        
        if xx['scale'].lower() == "log":
            ax.set_xscale("log")
        else:
            ax.xaxis.set_minor_locator(AutoMinorLocator())
        if yy['scale'].lower() == "log":
            ax.set_yscale("log")
        else:
            ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(self.load_lim(xx['lim']))
        ax.set_ylim(self.load_lim(yy['lim']))
        
        self.set_ticks_format(ax)
        ax.set_xlabel(r"{}".format(xx["label"]), **self.para['xlabel'])
        ax.set_ylabel(r"{}".format(yy["label"]), **self.para['ylabel'])
        # Draw any configured lines
        self.draw_lines(ax)
        # Draw any configured text labels
        self.draw_text_labels(ax)
        self.savefig(fig)
        self.serialize(fig, ax)
        if self.yaml['Plot_Config']['screen_show']: 
            plt.show()
 

                    
class ScatterC(Figure):
    def __init__(self, config, name, sts):
        super().__init__(config=config, name=name, sts=sts)
        with open(os.path.join(pwd, "cards/scatter_color.json"), "r") as f1: 
            setting = json.loads(f1.read())
            self.para = setting[self.style]

    def load_color_lim(self, cc):
        if cc.get("lim", None) == None:
            if "data" in cc.keys(): 
                cc['lim'] = [
                    min(cc['data']), 
                    max(cc['data'])
                ]
                

    def draw(self): 
        plt.close() 
        print(self.data.shape)
        self.legend = self.info.get("legend", False)

        fig = plt.figure(**self.para['figure'])
        ax  = fig.add_axes(**self.para["ax"])
        axc = fig.add_axes(**self.para["axc"])
        for spine in ax.spines.values():
            spine.set_linewidth(self.para["spines"]['linewidth'])   
        self.load_selection()
        xx = self.info["parameters"][0]
        yy = self.info["parameters"][1]
        cc = self.info["color"]
        xx['data'] = self.load_vars(xx['expr'])   
        yy['data'] = self.load_vars(yy['expr'])
        print(cc)
        cc['data'] = self.load_vars(cc['expr'])
        self.load_color_lim(cc)

        if self.info.get("group", []):
            nn = 0 
            from itertools import cycle
            st = cycle(deepcopy(self.para['scatter']))
            for item in self.info['group']: 
                dfs = self.load_bool_df(self.data, item['condition'])
                xi = self.load_vars(xx['expr'], dfs)
                yi = self.load_vars(yy['expr'], dfs)
                ci = self.load_vars(cc['expr'], dfs)
                if item.get("use_cmap", True):
                    style = self.info.get("scatter_style", {})
                    a = deepcopy(self.para['scatterc'])
                    a.update(style)
                    style = a 
                    if item.get("style", {}):
                        item['style'].pop("color", None)
                        item['style'].pop("facecolor", None)
                        style.update(item['style'])
                    if cc.get('scale', "linear").lower() == "log":
                        from matplotlib.colors import LogNorm 
                        c1 = ax.scatter(xi, yi, c=ci, **style, norm=LogNorm(cc['lim'][0], cc['lim'][1]), label=r"{}".format(item['name']))
                    else: 
                        c1 = ax.scatter(xi, yi, c=ci, **style, vmin=cc['lim'][0], vmax=cc['lim'][1], label=r"{}".format(item['name']))
                    # ax.scatter(xi, yi, **style)
                else:
                    xi = self.load_vars(xx['expr'], dfs)
                    yi = self.load_vars(yy['expr'], dfs)
                    style = next(st)
                    style.update(self.info.get("scatter_style", {}))
                    style.pop("cmap", None)
                    style.update(item.get("style", {}))
                    ax.scatter(xi, yi, **style, label=r"{}".format(item['name']))
                    nn += 1
        else: 
            style = self.info.get("scatter_style", {})
            self.para['scatterc'].update(style)
            if cc.get('scale', "linear").lower() == "log":
                from matplotlib.colors import LogNorm
                c1 = ax.scatter(xx['data'], yy['data'], c=cc['data'], 
                                **self.para['scatterc'], 
                                norm=LogNorm(cc['lim'][0], vmax=cc['lim'][1]),
                                )
            else:
                c1 = ax.scatter(xx['data'], yy['data'], c=cc['data'], 
                                **self.para['scatterc'], 
                                vmin=cc['lim'][0], vmax=cc['lim'][1]
                                )
        
        if cc.get('scale', "linear").lower() == "log":
            from matplotlib.colors import LogNorm
            plt.colorbar(c1, axc, norm=LogNorm(cc['lim'][0], cc['lim'][1]), orientation='vertical', extend='neither')
            axc.set_yscale("log")
        else: 
            plt.colorbar(c1, axc, orientation="vertical", extend="neither")
            axc.yaxis.set_minor_locator(AutoMinorLocator())
        
        axc.tick_params(**self.para["axc_ticks"])
        axc.tick_params(**self.para["axc_ticks_major"])
        axc.tick_params(**self.para["axc_ticks_minor"])
        
        # for spine in axc.spines.values():
        #     spine.set_linewidth(self.para["spines"]['linewidth'])   
                                    
        if xx['scale'].lower() == "log":
            ax.set_xscale("log")
        else:
            ax.xaxis.set_minor_locator(AutoMinorLocator())
        if yy['scale'].lower() == "log":
            ax.set_yscale("log")
        else:
            ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_xlim(self.load_lim(xx['lim']))
        ax.set_ylim(self.load_lim(yy['lim']))
        
        self.set_ticks_format(ax)
        ax.set_xlabel(r"{}".format(xx["label"]), **self.para['xlabel'])
        ax.set_ylabel(r"{}".format(yy["label"]), **self.para['ylabel'])
        axc.set_ylabel(r"{}".format(cc["label"]), **self.para["ylabel"])
        # Draw any configured lines
        self.draw_lines(ax)
        # Draw any configured text labels
        self.draw_text_labels(ax)
        if self.legend: 
            ax.legend(**self.legend)
        
        self.savefig(fig)
        
        # self.serialize(fig, ax)
        if self.yaml['Plot_Config']['screen_show']: 
            plt.show()


class VoronoiC(Figure):

    def __init__(self, config, name, sts):
        super().__init__(config=config, name=name, sts=sts)
        with open(os.path.join(pwd, "cards/voronoi_color.json"), "r") as f1: 
            setting = json.loads(f1.read())
            self.para = setting[self.style]


    def load_color_lim(self, cc):
        if cc.get("lim", None) == None:
            if "data" in cc.keys(): 
                cc['lim'] = [
                    min(cc['data']), 
                    max(cc['data'])
                ]
    
    def create_vars_data(self, xx):
        xl = self.load_lim(xx['lim'])
        if xx['scale'].lower() == "log":
            xa = np.log10(xx['data'].to_numpy())
            xa = (xa - np.log10(xl[0])) / (np.log10(xl[1]) - np.log10(xl[0]))
        else: 
            xa = xx['data'].to_numpy()
            xa = (xa - xl[0]) / (xl[1] - xl[0])
        return xa 
    
    def draw(self):
        plt.close() 
        from shapely.ops import unary_union
        from shapely.geometry import Polygon, MultiPolygon
        fig = plt.figure(**self.para['figure'])
        ax  = fig.add_axes(**self.para["ax"])
        axc = fig.add_axes(**self.para["axc"])
        for spine in ax.spines.values():
            spine.set_linewidth(self.para["spines"]['linewidth'])   

        self.load_selection()
        xx = self.info["parameters"][0]
        yy = self.info["parameters"][1]
        cc = self.info["color"]
        xx['data'] = self.load_vars(xx['expr'])   
        yy['data'] = self.load_vars(yy['expr'])
        cc['data'] = self.load_vars(cc['expr'])
        self.load_color_lim(cc)

        # --- group-selection logic; grey out unselected cells ---
        group = self.info.get("group", [])
        self.legend = self.info.get("legend", False)
        from scipy.spatial import Voronoi
        xa = self.create_vars_data(xx)
        ya = self.create_vars_data(yy)
        points = np.column_stack((xa, ya))
        self.get_cmap(cc)
        
        if group:
            # start with all False
            vor = Voronoi(points)
            regions, vertices = self.voronoi_finite_polygons_2d(vor)
            if self.use_wall:
                for region in regions:
                    polygon = vertices[region]
                    polygon_closed = np.concatenate([polygon, polygon[:1]], axis=0)
                    ax.plot(polygon_closed[:, 0], polygon_closed[:, 1], transform=ax.transAxes, **self.para['wall'])
            
            # mark rows matching any group condition
            for item in group:
                self.data['_inner_selected'] = False
                mask = self.load_bool_df(self.data, item['condition'])
                print("Loading subgroup -> ", item["name"], mask.shape)
                self.data.loc[mask.index, '_inner_selected'] = True
                point_labels = list(self.data.index)
                lbstag = True
                # collect styled polygons for this subgroup
                styled_polygons = []
                for ii, region in enumerate(regions):
                    polygon = vertices[region]
                    if self.data['_inner_selected'].iat[ii]:
                        if item['use_cmap']: 
                            facecolor = self.get_color(cc['data'].loc[point_labels[ii]])
                            # print("drawing cell -> {} with {}".format(ii, facecolor), end="\r", flush=True)
                            ax.fill(polygon[:, 0], polygon[:, 1],
                                facecolor=facecolor,
                                transform=ax.transAxes,
                                **self.cstyle['fill'])
                        if item.get("style", False):
                            
                            # accumulate this cell's polygon for later union
                            styled_polygons.append(Polygon(polygon))
                # merge all styled cells and fill once
                if item.get("style", False) and styled_polygons:
                    merged = unary_union(styled_polygons)
                    # after you compute `merged = unary_union(polygons)`
                    if isinstance(merged, Polygon):
                        polys = [merged]
                    elif isinstance(merged, MultiPolygon):
                        polys = list(merged.geoms)
                    else:
                        # fallback to anything with .geoms, or just wrap it
                        polys = list(getattr(merged, 'geoms', [merged]))
                    # polys = [merged] if isinstance(merged, Polygon) else list(merged)
                    cstyle = deepcopy(self.cstyle['fill'])
                    cstyle.update(item['style'])
                    for poly in polys:
                        xs, ys = poly.exterior.xy
                        if lbstag:
                            ax.fill(xs, ys, transform=ax.transAxes,
                                    **cstyle, label=r"{}".format(item['name']))
                            lbstag = False
                        else: 
                            ax.fill(xs, ys, transform=ax.transAxes,
                                    **cstyle)
                            
                            
                            
                            # # print(item.get("style", False))
                            # cstyle = deepcopy(self.cstyle['fill'])
                            # cstyle.update(item['style'])
                            # # print("drawing cell -> {} with {}".format(ii, item.get("style", False)['facecolor']), end="\r", flush=True)
                            # # time.sleep(0.001)
                            # if lbstag:
                            #     ax.fill(polygon[:, 0], polygon[:, 1],
                            #             transform=ax.transAxes, hatch='\\',
                            #             **cstyle, label=r"{}".format(item['name']))
                            #     lbstag = False 
                            # else: 
                            #     ax.fill(polygon[:, 0], polygon[:, 1],
                            #             transform=ax.transAxes,
                            #             **cstyle)
        else:
            vor = Voronoi(points)
            regions, vertices = self.voronoi_finite_polygons_2d(vor)
            if self.use_wall:
                for region in regions:
                    polygon = vertices[region]
                    polygon_closed = np.concatenate([polygon, polygon[:1]], axis=0)
                    ax.plot(polygon_closed[:, 0], polygon_closed[:, 1], transform=ax.transAxes, **self.para['wall'])
            for i, region in enumerate(regions):
                polygon = vertices[region]
                # grey fill if not selected, otherwise use colormap
                ax.fill(polygon[:, 0], polygon[:, 1],
                        facecolor=self.get_color(cc['data'][i]),
                        transform=ax.transAxes,
                        **self.cstyle['fill'])
        
        if self.legend: 
            ax.legend(**self.legend)


        # voronoi_plot_2d(vor, ax=ax, show_vertices=False, show_points=False, line_colors='black',
                # line_width=0.5, line_alpha=0.6, zorder=1)
        # if self.info.get("group", []):
        #     pass 
        # else: 
        if self.use_core:
            ax.scatter(xx['data'], yy['data'], **self.corestyle)
        # ax.set_xscale("log")
        import matplotlib.cm as cm
        plt.colorbar(cm.ScalarMappable(norm=self.norm, cmap=self.cmap), cax=axc)
        
        for spine in axc.spines.values():
            spine.set_linewidth(self.para["spines"]['linewidth'])  
        axc.tick_params(**self.para["axc_ticks"])
        axc.tick_params(**self.para["axc_ticks_major"])
        axc.tick_params(**self.para["axc_ticks_minor"])

        if xx['scale'].lower() == "log":
            ax.set_xscale("log")
        else:
            ax.xaxis.set_minor_locator(AutoMinorLocator())
        if yy['scale'].lower() == "log":
            ax.set_yscale("log")
        else:
            ax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.set_xlim(self.load_lim(xx['lim']))
        ax.set_ylim(self.load_lim(yy['lim']))
        
        self.set_ticks_format(ax)
        ax.set_xlabel(r"{}".format(xx["label"]), **self.para['xlabel'])
        ax.set_ylabel(r"{}".format(yy["label"]), **self.para['ylabel'])
        axc.set_ylabel(r"{}".format(cc["label"]), **self.para["ylabel"])
 
        # Draw any configured lines
        self.draw_lines(ax)
        # Draw any configured text labels
        self.draw_text_labels(ax)
        self.savefig(fig)
        if self.yaml['Plot_Config']['screen_show']: 
            plt.show()


    def voronoi_finite_polygons_2d(self, vor, radius=None):
        """
        https://gist.github.com/pv/8036995
        """
        if vor.points.shape[1] != 2:
            raise ValueError("Requires 2D input")

        new_regions = []
        new_vertices = vor.vertices.tolist()

        center = vor.points.mean(axis=0)
        if radius is None:
            radius = vor.points.ptp().max() * 2

        # 构造每个点对应的所有边（ridge）的字典
        all_ridges = {}
        for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
            all_ridges.setdefault(p1, []).append((p2, v1, v2))
            all_ridges.setdefault(p2, []).append((p1, v1, v2))

        for p1, region_index in enumerate(vor.point_region):
            vertices_idx = vor.regions[region_index]
            if all(v >= 0 for v in vertices_idx):
                new_regions.append(vertices_idx)
                continue
            
            ridges = all_ridges[p1]
            new_region = [v for v in vertices_idx if v >= 0]

            for p2, v1, v2 in ridges:
                if v2 < 0:
                    v1, v2 = v2, v1
                if v1 >= 0 and v2 >= 0:
                    continue
                t = vor.points[p2] - vor.points[p1]
                t = t / np.linalg.norm(t)
                n = np.array([-t[1], t[0]])
                midpoint = vor.points[[p1, p2]].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                far_point = vor.vertices[v2] + direction * radius
                new_region.append(len(new_vertices))
                new_vertices.append(far_point.tolist())

            # 将区域顶点按逆时针排序
            vs = np.asarray([new_vertices[v] for v in new_region])
            c = vs.mean(axis=0)
            angles = np.arctan2(vs[:,1]-c[1], vs[:,0]-c[0])
            new_region = np.array(new_region)[np.argsort(angles)]
            new_regions.append(new_region.tolist())

        return new_regions, np.asarray(new_vertices)

    def get_cmap(self, cc):
        import matplotlib.cm as cm
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        self.cstyle = deepcopy(self.para['cell'])
        vstyle = self.info.get("voronoi_style", {})
        self.cmap = cm.get_cmap(self.para['cell']['cmap'])
        self.corestyle = deepcopy(self.para['core'])
        self.use_core = vstyle.get("use_core", self.para['cell']['use_core'])
        self.use_wall = vstyle.get("use_wall", self.para['cell']['use_wall'])
        if vstyle.get("cell", None):
            if vstyle.get("cmap", False) or vstyle.get("fill", False): 
                self.cstyle.update(vstyle)
                try:
                    self.cmap = cm.get_cmap(vstyle.get("cmap", None))
                except: 
                    self.cmap = cm.get_cmap(self.para['cell']['cmap'])
            if vstyle.get("core"): 
                self.corestyle.update(vstyle['core'])

        if cc['scale'].lower() != "log": 
            self.norm = Normalize(vmin=cc['lim'][0], vmax=cc['lim'][1])
        else: 
            self.norm = LogNorm(vmin=cc['lim'][0], vmax=cc['lim'][1])
    
    def get_color(self, cvalue):
        return self.cmap(self.norm(cvalue))

def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))