import os
import fiona
from pathlib import Path
import geopandas as gpd
import pandas as pd
#from bokeh.models import GeoJSONDataSource
#from bokeh.plotting import figure, show, output_file
#from bokeh.tile_providers import get_provider, Vendors
import matplotlib.pyplot as plt

from shapely.geometry import Point, LineString, Polygon, MultiLineString
from shapely import ops


#import os
#import json
#from datetime import datetime
#from pathlib import Path
#from functools import partial
#import tkinter as tk
#from tkinter import filedialog
#import numpy as np
#import pandas as pd
#import geopandas as gpd
#import pyproj
#from shapely.ops import transform, unary_union
#from cartopy.geodesic import Geodesic


class Enc:

    def __init__(self, path):
        self.path = path
        self.name = self.path.name
        self.band = self.name[2:5]

    def list_layers(self):
        return fiona.listlayers(str(self.path))

    def __str__(self):
        return str(self.path)

    def get_objects(self, acronym, geom_type):
        try:
            objects = gpd.read_file(str(self.path), layer=acronym)
            objects = objects[objects.geom_type == geom_type]
        
            if acronym is 'SLCONS' and geom_type is 'LineString':
                #no_piers = objects['CATSLC'] != 4  # pier (jetty)
                #no_fenders = objects['CATSLC'] != 14  # fender
                #not_submerged = objects['WATLEV'] != 3  # always under water/submerged
                #objects = objects[no_piers & no_fenders & not_submerged]

                piers = objects['CATSLC'] == 4  
                fenders = objects['CATSLC'] == 14  
                submerged = objects['WATLEV'] == 3 
                objects = objects[piers | fenders | submerged]

            if acronym is 'M_COVR':
                no_coverage = objects['CATCOV'] != 2
                objects = objects[no_coverage]

            return objects

        except Exception as e:
            print(e)


class Shorex:

    wgs84 = {'init': 'epsg:4326'}
    webmerc = {'init': 'epsg:3857'}

    def __init__(self, enc_dir, ref_dir, bands_to_process):
        self.enc_dir = enc_dir
        self.ref_dir = ref_dir
        self.geojsons = []
        self.bands_to_process = bands_to_process
        self.enc_objects = self.ref_dir / 'ENC_Objects.gpkg'
        self.ref_bands_gpkg_path = self.ref_dir / 'Reference_Bands.gpkg'
        self.cusp_ref_gpkg_path = self.ref_dir / 'CUSP_Reference.gpkg'
        self.sliver_tolerance = 0.0000005

    @staticmethod
    def set_env_vars(env_name):
        user_dir = os.path.expanduser('~')
        conda_dir = Path(user_dir).joinpath('AppData', 'Local', 
                                            'Continuum', 'anaconda3')
        env_dir = conda_dir / 'envs' / env_name
        share_dir = env_dir / 'Library' / 'share'
        script_path = conda_dir / 'Scripts'
        gdal_data_path = share_dir / 'gdal'
        proj_lib_path = share_dir

        if script_path.name not in os.environ["PATH"]:
            os.environ["PATH"] += os.pathsep + str(script_path)
        os.environ["GDAL_DATA"] = str(gdal_data_path)
        os.environ["PROJ_LIB"] = str(proj_lib_path)

    def export_to_geojson(self, objects):
        for geom_type, v in objects.items():
            for acronym, objects_list in v.items():
                acronym = acronym.replace('_', '')
                geojson_name = '{}_{}.geojson'.format(acronym, geom_type)
                geojson_path = self.enc_dir / geojson_name
                object_df = pd.concat(objects_list, ignore_index=True)
                object_gdf = gpd.GeoDataFrame(
                    object_df, crs=self.wgs84).to_crs(self.webmerc)
                object_gdf.to_file(str(geojson_path), driver='GeoJSON')
                self.geojsons.append(geojson_path)

    def export_to_gpkg(self, objects, band):
        for geom_type, v in objects.items():
            for acronym, objects_list in v.items():
                if not all(o is None for o in objects_list):  # if not all None
                    df = pd.concat(objects_list)
                    df['band'] = [band] * df.shape[0]
                    gdf = gpd.GeoDataFrame(df, geometry='geometry', 
                                            crs=self.wgs84).explode()
                    layer = '{}_{}_Band{}'.format(acronym, geom_type, band)
                
                    try:
                        gdf.to_file(self.enc_objects, layer=layer, driver='GPKG')
                    except Exception as e:
                        print(e, '({} - {})'.format(acronym, geom_type))

    def plot_objects(self):

        c = {
            'Point': {'OBSTRN': 'red',},

            'LineString': {'COALNE': 'blue',
                           'SLCONS': 'green',},

            'Polygon': {'LNDARE': 'orange',
                        'MCOVR': 'lightgray',
                        'MQUAL': 'gray',
                        'DEPARE': 'blue',}
             }

        p = figure(title='ENC Objects', 
                   x_axis_type="mercator", y_axis_type="mercator",
                   plot_width=1500, plot_height=900)

        tile_provider = get_provider(Vendors.CARTODBPOSITRON)
        p.add_tile(tile_provider)

        for g in self.geojsons:
            acronym, geom_type = g.stem.split('_')
            with open(str(g)) as f:
                geojson_data = f.read()

            source = GeoJSONDataSource(geojson=geojson_data)
            
            if geom_type == 'Polygon':
                p.patches('xs', 'ys', source=source, 
                          color=c[geom_type][acronym], alpha=0.4)

            elif geom_type == 'LineString':
                p.multi_line('xs', 'ys', source=source, 
                             color=c[geom_type][acronym], line_width=2)

            elif geom_type == 'Point':
                p.circle('x', 'y', source=source, color=c[geom_type][acronym])

        output_file("tile.html")
        show(p)

    #def create_ref_shoreline(self, objects, band):
    #    coalne = objects['LineString']['COALNE'], 
    #    slcons = objects['LineString']['SLCONS'],
    #    rivers = objects['LineString']['RIVERS'],

    #    coalne_df = pd.concat(coalne, ignore_index=True)
    #    slcons_df = pd.concat(slcons, ignore_index=True)
    #    rivers_df = pd.concat(rivers, ignore_index=True)

    #    coalne_gdf = gpd.GeoDataFrame(coalne_df.geometry, crs=self.wgs84)
    #    slcons_gdf = gpd.GeoDataFrame(slcons_df.geometry, crs=self.wgs84)
    #    rivers_gdf = gpd.GeoDataFrame(rivers_df.geometry, crs=self.wgs84)

    #    shoreline = coalne_gdf.geometry.union(slcons_gdf.geometry)
    #    shoreline = shoreline.union(rivers_gdf.geometry)
    #    gdf = gpd.GeoDataFrame(geometry=shoreline, crs=self.wgs84).explode()
    #    layer = 'CUSP_Reference_Band{}'.format(band)
    #    gdf.to_file(self.ref_bands_gpkg_path, layer=layer, driver='GPKG')

    def create_ref_shoreline_2(self, objects, band):

        lndare = objects['Polygon']['LNDARE']
        mcovr = objects['Polygon']['M_COVR']
        slcons = objects['LineString']['SLCONS']

        lndare_df = pd.concat(lndare, ignore_index=True)
        mcovr_df = pd.concat(mcovr, ignore_index=True)
        slcons_df = pd.concat(slcons, ignore_index=True)

        lndare_gdf = gpd.GeoDataFrame(lndare_df.geometry, crs=self.wgs84).explode().reset_index()
        mcovr_gdf = gpd.GeoDataFrame(mcovr_df.geometry, crs=self.wgs84).explode().reset_index()
        slcons_gdf = gpd.GeoDataFrame(slcons_df.geometry, crs=self.wgs84).explode().reset_index()
        
        enc_polys = mcovr_gdf

        lakare = objects['Polygon']['LAKARE']
        if not all(o is None for o in lakare):  # if not all None
            lakare_df = pd.concat(lakare, ignore_index=True)
            lakare_gdf = gpd.GeoDataFrame(lakare_df.geometry, crs=self.wgs84).explode().reset_index()
            cols = lndare_gdf.columns
            lndare_gdf = lndare_gdf.drop([c for c in cols if c is not 'geometry'], axis=1)
            lndare_gdf = gpd.overlay(lndare_gdf, lakare_gdf, how='difference').explode().reset_index()

        rivers = objects['Polygon']['RIVERS']
        if not all(o is None for o in rivers):  # if not all None
            rivers_df = pd.concat(rivers, ignore_index=True)
            rivers_gdf = gpd.GeoDataFrame(rivers_df.geometry, crs=self.wgs84).explode().reset_index()
            cols = lndare_gdf.columns
            lndare_gdf = lndare_gdf.drop([c for c in cols if c is not 'geometry'], axis=1)
            lndare_gdf = gpd.overlay(lndare_gdf, rivers_gdf, how='difference').explode().reset_index()

        
        mcscl = objects['Polygon']['M_CSCL']
        if not all(o is None for o in mcscl):  # if not all None
            mcscl_df = pd.concat(mcscl, ignore_index=True)
            mcscl_gdf = gpd.GeoDataFrame(mcscl_df.geometry, crs=self.wgs84).explode().reset_index()
            enc_polys = gpd.overlay(enc_polys, mcscl_gdf, how='union').explode().reset_index()

        print('intersecting LNDARE with M_COVR to get pre-reference...')
        sindex_lndare = lndare_gdf.sindex
        sindex_slcons = slcons_gdf.sindex
        shoreline_bits = []
        
        num_enc_polys = enc_polys.shape[0]
        for i, poly in enumerate(enc_polys.geometry, 1):
            print('enc_poly {} of {}...'.format(i, num_enc_polys))
            possible_lndare_idx = list(sindex_lndare.intersection(poly.bounds))
            possible_lndares = lndare_gdf.iloc[possible_lndare_idx]

            precise_lndares = possible_lndares.intersection(poly)
            precise_lndares = precise_lndares[precise_lndares.geom_type == 'Polygon']
            precise_lndares = precise_lndares[~precise_lndares.is_empty]

            num_landare = precise_lndares.shape[0]
            for j, lndare in enumerate(precise_lndares.geometry, 1):
                print('lndare {} of {} (enc_poly {} of {})...'.format(j, num_landare, i, num_enc_polys))
                #try:
                possible_slcons_idx = list(sindex_slcons.intersection(lndare.bounds))
                possible_slcons = slcons_gdf.iloc[possible_slcons_idx]

                precise_slcons = possible_slcons.intersection(lndare).explode()
                precise_slcons = precise_slcons[~precise_slcons.is_empty]
                precise_slcons = precise_slcons[precise_slcons.geom_type == 'LineString']

                precise_slcons_geom = MultiLineString([g for g in precise_slcons.geometry])
                slcons = gpd.GeoDataFrame(geometry=[precise_slcons_geom], crs=self.wgs84)

                lndare_poly_gdf = gpd.GeoDataFrame(geometry=[lndare], crs=self.wgs84)
                if precise_slcons.shape[0] > 0:
                    clipped_shoreline = lndare_poly_gdf.exterior.difference(slcons)

                    if not clipped_shoreline.iloc[0].is_empty:
                        verts = []

                        for l in clipped_shoreline.explode():
                            verts.extend(list(l.coords))
                        verts.append(verts[0])

                        shoreline_geom = LineString(verts)
                        shoreline_poly_linestring = gpd.GeoDataFrame(geometry=[shoreline_geom], crs=self.wgs84)

                        shoreline_line = shoreline_poly_linestring.difference(poly.boundary)
                        shoreline_bits.append(shoreline_line.explode().reset_index())
                else:
                    shoreline_line = lndare_poly_gdf.boundary.difference(poly.boundary)
                    shoreline_bits.append(shoreline_line.explode().reset_index())
            
        print('saving results...')
        df = pd.concat(shoreline_bits, ignore_index=True)
        print(df)
        gdf = gpd.GeoDataFrame(geometry=df[0], crs=self.wgs84)
        gdf = gdf[gdf.geom_type == 'LineString']
        
        gdf.to_file(str(self.ref_bands_gpkg_path), 
                    layer='CUSP_Reference_Band{}'.format(band),
                    driver='GPKG')

    def gen_cusp_band_regions(self):

        bands = []
        for b in self.bands_to_process:
            layer= r'M_COVR_Polygon_Band{}'.format(b)
            gdf = gpd.read_file(str(self.enc_objects), layer=layer)
            bands.append(gdf)

        bands.reverse()  # [6, 5, 4, 3, 2, 1]
        current_band_coverage = bands[0]  # i.e., largest-scale band
        for b in bands[1:6]:  # start with largest-scale band and "work down"
            b_not_higher = gpd.overlay(b, current_band_coverage, how='difference').explode()
            b_not_higher = b_not_higher[b_not_higher.geometry.area > self.sliver_tolerance]
            df = pd.concat([current_band_coverage, b_not_higher], ignore_index=True)
            current_band_coverage = gpd.GeoDataFrame(df, geometry='geometry', crs=self.wgs84)
        
        layer = 'CUSP_band_regions'
        cusp_band_regions = current_band_coverage.explode()
        cusp_band_regions.to_file(self.cusp_ref_gpkg_path, 
                                  layer=layer, driver='GPKG')
        
    def create_cusp_ref(self):
        band_regions = gpd.read_file(str(self.cusp_ref_gpkg_path), 
                                     layer='CUSP_band_regions')
        
        for b in self.bands_to_process:
            print('reading Band {} CUSP reference shoreline...'.format(b))
            b_ref =  gpd.read_file(str(self.ref_bands_gpkg_path), 
                                   layer='CUSP_Reference_Band{}'.format(b))

            sindex = b_ref.sindex
            lines = []

            print('clipping Band {} reference with band {} regions...'.format(b, b))
            b_regions = band_regions[band_regions['band'] == b]

            for poly in b_regions.geometry:
                possible_lines_idx = list(sindex.intersection(poly.bounds))
                possible_lines = b_ref.iloc[possible_lines_idx]
                precise_lines = possible_lines.intersection(poly)
                lines.append(precise_lines.explode())
            
            print('saving results...')
            df = pd.concat(lines, ignore_index=True)
            gdf = gpd.GeoDataFrame(geometry=df)
            gdf = gdf[gdf.geom_type == 'LineString']
            gdf.crs = self.wgs84
            gdf.to_file(str(self.cusp_ref_gpkg_path), 
                        layer='CUSP_Reference_Band{}_CLIPPED'.format(b))


#             coords = list(polygon.exterior.geometry.iloc[0].coords)
#gdf = g.GeoDataFrame(geometry=[Polygon(coords[2:])], crs=crs)


def main():

    enc_dir = Path(r'C:\ENCs\All_ENCs\ENC_ROOT')
    ref_dir = Path(r'Z:\ENC_Shoreline_Extractor')
    bands_to_process = [1, 2, 3, 4, 5, 6]
    shorex = Shorex(enc_dir, ref_dir, bands_to_process)
    shorex.set_env_vars('enc_shorex')

    #for band in shorex.bands_to_process:
    #    print('processing band {}...'.format(band))
    #    objects_to_extract = {'Point': [],
    #                          'LineString': ['COALNE', 'SLCONS'],
    #                          'Polygon': ['LNDARE', 'LNDRGN', 
    #                                      'M_COVR', 'M_CSCL',
    #                                      'RIVERS', 'LAKARE']}

    #    objects = {'Point': {}, 'LineString': {}, 'Polygon': {}}

    #    encs = list(shorex.enc_dir.rglob('US{}*.000'.format(band)))
        
    #    num_encs = len(encs)

    #    for i, enc in enumerate(encs, 1):
    #        enc = Enc(enc)
    #        print('{} ({} of {})'.format(enc, i, num_encs))

    #        for geom_type in objects_to_extract.keys():
    #            for acronym in objects_to_extract[geom_type]:
    #                enc_objs = enc.get_objects(acronym, geom_type)

    #                if acronym not in objects[geom_type].keys():
    #                    objects[geom_type].update({acronym: []})

    #                objects[geom_type][acronym].append(enc_objs)

    #    shorex.export_to_gpkg(objects, band)
    #    #shorex.export_to_geojson(objects)
    #    #shorex.plot_objects()

    #    print(objects)

    #    shorex.create_ref_shoreline_2(objects, band)

    shorex.gen_cusp_band_regions()
    shorex.create_cusp_ref()


if __name__ == '__main__':
    main()

    #user_dir = os.path.expanduser('~')
    #conda_dir = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 'anaconda3')
    #env_dir = conda_dir / 'envs' / 'cusp'
    #share_dir = env_dir / 'Library' / 'share'
    #script_path = conda_dir / 'Scripts'
    #gdal_data_path = share_dir / 'gdal'
    #proj_lib_path = share_dir

    #if script_path.name not in os.environ["PATH"]:
    #    os.environ["PATH"] += os.pathsep + str(script_path)
    #os.environ["GDAL_DATA"] = str(gdal_data_path)
    #os.environ["PROJ_LIB"] = str(proj_lib_path)

    #band_dir = Path(r'C:\Users\Nick.Forfinski-Sarko\Documents\ArcGIS\Projects\ENC_Shoreline_Extractor\ENC_Shoreline_Extractor.gdb')

    #geod = Geodesic()

    #def distance(pt1, pt2):  # from https://pelson.github.io/2018/coast-path/
    #    result = np.array(geod.inverse(np.asanyarray(pt1), np.asanyarray(pt2)))
    #    return result[:, 0]


    #def linestring_distance(geom):  # from https://pelson.github.io/2018/coast-path/
    #    if hasattr(geom, 'geoms'):
    #        return sum(linestring_distance(subgeom) for subgeom in geom.geoms)
    #    else:
    #        points = np.array(geom.coords)
    #        return distance(points[:-1, :2], points[1:, :2]).sum()


    #length = {}
    #for b in range(1, 2):

    #    band_gdf = gpd.read_file(str(band_dir), layer='_70k')
        
    #    print('70k...'.format(b))
    #    print('-' * 50)

    #    band_length = 0
    #    for geom in band_gdf.geometry:
    #        band_length += linestring_distance(geom)

    #    length[str(b)] = band_length
    #    print(length)


    #    print()

    #length_df = pd.DataFrame(length, index=[0])
    #print(length_df.values)