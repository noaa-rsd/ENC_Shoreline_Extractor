import os
import fiona
from pathlib import Path
import geopandas as gpd
import pandas as pd
from bokeh.models import GeoJSONDataSource
from bokeh.plotting import figure, show, output_file
from bokeh.tile_providers import get_provider, Vendors


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
                no_piers = objects['CATSLC'] != 4  # pier (jetty)
                no_fenders = objects['CATSLC'] != 14  # fender
                not_submerged = objects['WATLEV'] != 3  # always under water/submerged
                objects = objects[no_piers & no_fenders & not_submerged]

            if acronym is 'M_COVR':
                no_coverage = objects['CATCOV'] != 2
                objects = objects[no_coverage]

            return objects

        except Exception as e:
            print(e)


class Shorex:

    wgs84 = {'init': 'epsg:4326'}
    webmerc = {'init': 'epsg:3857'}

    def __init__(self, enc_dir, ref_dir):
        self.enc_dir = enc_dir
        self.ref_dir = ref_dir
        self.geojsons = []
        self.gpkg_path = self.ref_dir / 'CUSP.gpkg'

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
                df = pd.concat(objects_list)
                df['band'] = [band] * df.shape[0]
                gdf = gpd.GeoDataFrame(df, geometry='geometry', 
                                       crs=self.wgs84).explode()
                layer = '{}_{}_Band{}'.format(acronym, geom_type, band)
                gdf.to_file(self.gpkg_path, layer=layer, driver='GPKG')

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

    def create_ref_shoreline(self, coalne, slcons, band):
        coalne_df  = pd.concat(coalne, ignore_index=True)
        slcons_df  = pd.concat(slcons, ignore_index=True)

        coalne_gdf = gpd.GeoDataFrame(coalne_df.geometry, crs=self.wgs84)
        slcons_gdf = gpd.GeoDataFrame(slcons_df.geometry, crs=self.wgs84)

        shoreline = coalne_gdf.geometry.union(slcons_gdf.geometry)
        gdf = gpd.GeoDataFrame(geometry=shoreline, crs=self.wgs84).explode()
        layer = 'CUSP_Reference_Band{}'.format(band)
        gdf.to_file(self.gpkg_path, layer=layer, driver='GPKG')

    def gen_cusp_band_regions(self, bands_to_process):

        bands = []
        for b in bands_to_process:
            layer= r'M_COVR_Polygon_Band{}'.format(b)
            gdf = gpd.read_file(str(self.gpkg_path), layer=layer)
            bands.append(gdf)

        gdf = gpd.GeoDataFrame(pd.concat(bands).geometry, 
                               crs=self.wgs84).explode()
        b1, b2, b3, b4, b5 = bands

        b4_not_5 = gpd.overlay(b4, b5, how='difference')
        df = pd.concat([b5, b4_not_5], ignore_index=True)
        b54 = gpd.GeoDataFrame(df, geometry='geometry', crs=self.wgs84)

        b3_not_45 = gpd.overlay(b3, b54, how='difference')
        df = pd.concat([b54, b3_not_45], ignore_index=True)
        b543 = gpd.GeoDataFrame(df, geometry='geometry', crs=self.wgs84)

        b2_not_345 = gpd.overlay(b2, b543, how='difference')
        df = pd.concat([b543, b2_not_345], ignore_index=True)
        b5432 = gpd.GeoDataFrame(df, geometry='geometry', crs=self.wgs84)

        b1_not_2345 = gpd.overlay(b1, b5432, how='difference')
        df = pd.concat([b5432, b1_not_2345], ignore_index=True)
        b54321 = gpd.GeoDataFrame(df, geometry='geometry', crs=self.wgs84)

        layer = 'CUSP_band_regions'
        cusp_band_regions = b54321.explode()
        cusp_band_regions.to_file(self.gpkg_path, layer=layer, driver='GPKG')


def main():

    enc_dir = Path(r'C:\ENCs\All_ENCs\ENC_ROOT')
    ref_dir = Path(r'Z:\ENC_Shoreline_Extractor')
    shorex = Shorex(enc_dir, ref_dir)
    shorex.set_env_vars('enc_shorex')

    bands_to_process = [1, 2, 3, 4, 5]

    for band in bands_to_process:
        print('processing band {}...'.format(band))
        objects_to_extract = {'Point': [],
                              'LineString': ['COALNE', 'SLCONS'],
                              'Polygon': ['M_COVR']}

        objects = {'Point': {}, 'LineString': {}, 'Polygon': {}}

        encs = list(shorex.enc_dir.rglob('US{}*.000'.format(band)))
        num_encs = len(encs)

        for i, enc in enumerate(encs, 1):
            enc = Enc(enc)
            print('{} ({} of {})'.format(enc, i, num_encs))

            for geom_type in objects_to_extract.keys():
                for acronym in objects_to_extract[geom_type]:
                    enc_objs = enc.get_objects(acronym, geom_type)

                    if acronym not in objects[geom_type].keys():
                        objects[geom_type].update({acronym: []})

                    objects[geom_type][acronym].append(enc_objs)

        shorex.export_to_gpkg(objects, band)
        #shorex.export_to_geojson(objects)
        #shorex.plot_objects()

        shorex.create_ref_shoreline(objects['LineString']['COALNE'], 
                                    objects['LineString']['SLCONS'],
                                    band)

    shorex.gen_cusp_band_regions(bands_to_process)


if __name__ == '__main__':
    main()

