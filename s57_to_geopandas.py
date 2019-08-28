import os
import fiona
from fiona._drivers import GDALEnv
from pathlib import Path
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import contextily as ctx
from bokeh.models import GeoJSONDataSource
from bokeh.plotting import figure, show, output_file
from bokeh.tile_providers import get_provider, Vendors


gdal_data =  Path(r'C:\Users\nickf\.conda\envs\s57\Library\share\gdal')
proj_lib =Path(r'C:\Users\nickf\.conda\envs\s57\Library\share\proj')

os.environ["GDAL_DATA"] = str(gdal_data)
os.environ["PROJ_LIB"] = str(proj_lib)

enc_dir = Path(r'C:\Users\nickf\Downloads\All_ENCs\ENC_ROOT')

gdfs = []

s57_geojson_path = Path(r'C:\Users\nickf\Downloads\All_ENCs\s57_DEPCNT.geojson')

for enc in enc_dir.rglob('US5*.000'):

    #enc = Enc()

    lyrs = fiona.listlayers(str(enc))
    print(enc)

    for lyr in lyrs:
        if lyr == 'M_QUAL':
            s57 = gpd.read_file(str(enc), layer=lyr)
            s57 = s57.to_crs({'init': 'epsg:3857'})
            s57 = s57.loc[s57['CATZOC'] == 6]
            gdfs.append(s57)

rdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs=gdfs[0].crs)
rdf.to_file(str(s57_geojson_path), driver='GeoJSON')

output_file("tile.html")
tile_provider = get_provider(Vendors.CARTODBPOSITRON)

with open(str(s57_geojson_path)) as f:
    geojson_data = f.read()

s57_source = GeoJSONDataSource(geojson=geojson_data)

p = figure(
    title="A test map", 
    x_range=(-2000000, 6000000), 
    y_range=(-1000000, 7000000), 
    x_axis_type="mercator", 
    y_axis_type="mercator",
    plot_width=1200,
    plot_height=600)
p.add_tile(tile_provider)
#p.multi_line('xs', 'ys', source=s57_source, color='gray', line_width=2)
#p.circle('x', 'y', source=s57_source, color='gray')
p.patches('xs', 'ys', source=s57_source, alpha=0.4)
##p.multi_polygons('xs', 'ys', source=s57_source, alpha=0.4)

outfp = r"C:\Users\nickf\Downloads\All_ENCs\test_map.html"
show(p)
