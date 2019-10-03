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
        self.name = path.name

    def list_layers(self):
        return fiona.listlayers(str(enc))

    def __str__(self):
        return str(self.path)

    def get_objects(self, acronym, geom_type):
        objects = gpd.read_file(str(enc), layer=acronym)
        objects = objects[objects.geom_type == geom_type]
        #objects = objects.loc[objects['CATZOC'] == 6]
        objects = objects.to_crs({'init': 'epsg:3857'})
        return objects


class Shorex:

    def __init__(self, enc_dir):
        self.enc_dir = enc_dir

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

    @staticmethod
    def export_to_geojson(gdfs, s57_geojson_path):
        rdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), 
                               crs=gdfs[0].crs)
        rdf.to_file(str(s57_geojson_path), driver='GeoJSON')
        print(rdf)

    @staticmethod
    def plot_objects(s57_geojson_path):
        output_file("tile.html")
        tile_provider = get_provider(Vendors.CARTODBPOSITRON)

        with open(str(s57_geojson_path)) as f:
            geojson_data = f.read()

        source = GeoJSONDataSource(geojson=geojson_data)

        p = figure(
            title="A test map", 
            x_axis_type="mercator", 
            y_axis_type="mercator",
            plot_width=1500,
            plot_height=900)

        p.add_tile(tile_provider)
        #p.multi_line('xs', 'ys', source=source, color='gray', line_width=2)
        #p.circle('x', 'y', source=source, color='gray')
        p.patches('xs', 'ys', source=source, alpha=0.4)
        ##p.multi_polygons('xs', 'ys', source=source, alpha=0.4)

        outfp = r"C:\Users\nickf\Downloads\All_ENCs\test_map.html"
        show(p)


if __name__ == '__main__':
    enc_dir = Path(r'Z:\ENCs')
    shorex = Shorex(enc_dir)
    shorex.set_env_vars('enc_shorex')
    gdfs = []

    for enc in shorex.enc_dir.rglob('US5*.000'):
        enc = Enc(enc)
        print(enc)
        enc_objs = enc.get_objects('LNDARE', 'Polygon')
        gdfs.append(enc_objs)

    geojson_path = enc_dir / 'ENC_{}_{}.geojson'.format('LNDARE', 'Polygon')
    shorex.export_to_geojson(gdfs, geojson_path)
    shorex.plot_objects(geojson_path)

