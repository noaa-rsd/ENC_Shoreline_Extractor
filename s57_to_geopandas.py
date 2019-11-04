import os
import numpy as np
from pathlib import Path
import pandas as pd
import geopandas as gpd
import fiona
from shapely.geometry import LineString, MultiLineString
from cartopy.geodesic import Geodesic


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
        self.sliver_tolerance = 0.0000005  # ? units
        self.slcons_tolerance = 250  # meters
        self.geod = Geodesic()

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

    def cookie_cut_inland_waters(self, lndare, lakare, rivers):
        lndare_df = pd.concat(lndare, ignore_index=True)
        lndare_gdf = gpd.GeoDataFrame(lndare_df.geometry, crs=self.wgs84).explode().reset_index()

        def difference_objects(lndare_gdf, objects):
            """differe LNDARE with other objects to get polygon the exterior
            and interior(s) of which are defined to be the shoreline"""

            if not all(o is None for o in objects):  # if not all None
                objects_df = pd.concat(objects, ignore_index=True)
                objects_gdf = gpd.GeoDataFrame(objects_df.geometry, crs=self.wgs84).explode().reset_index()
                cols = lndare_gdf.columns
                lndare_gdf = lndare_gdf.drop([c for c in cols if c is not 'geometry'], axis=1)
                lndare_gdf = gpd.overlay(lndare_gdf, objects_gdf, how='difference').explode().reset_index()            
            return lndare_gdf

        lndare_gdf = difference_objects(lndare_gdf, lakare)
        lndare_gdf = difference_objects(lndare_gdf, rivers)
        return lndare_gdf

    def union_mcscl(self, enc_polys, mcscl):
        if not all(o is None for o in mcscl):  # if not all None
            mcscl_df = pd.concat(mcscl, ignore_index=True)
            mcscl_gdf = gpd.GeoDataFrame(mcscl_df.geometry, crs=self.wgs84).explode().reset_index()
            enc_polys = gpd.overlay(enc_polys, mcscl_gdf, how='union').explode().reset_index()
        return enc_polys

    def save_shoreline(self, shoreline_bits, band):
        print('saving results...')
        df = pd.concat(shoreline_bits, ignore_index=True)
        gdf = gpd.GeoDataFrame(geometry=df[0], crs=self.wgs84)
        gdf = gdf[gdf.geom_type == 'LineString']
        gdf.to_file(str(self.ref_bands_gpkg_path), 
                    layer='CUSP_Reference_Band{}'.format(band),
                    driver='GPKG')

    def erase_slcons(self, poly, sindex_slcons, slcons_gdf, land_poly):

        def remove_short_linestrings(slcons):

            def distance(pt1, pt2):
                result = np.array(self.geod.inverse(np.asanyarray(pt1), np.asanyarray(pt2)))
                return result[:, 0]

            def linestring_distance(geom):  # from https://pelson.github.io/2018/coast-path/
                if hasattr(geom, 'geoms'):
                    return sum(linestring_distance(subgeom) for subgeom in geom.geoms)
                else:
                    points = np.array(geom.coords)
                    return distance(points[:-1, :2], points[1:, :2]).sum()

            # keep "large piers"
            lengths_m = np.asarray([linestring_distance(g) for g in list(slcons)])
            keep_indx = lengths_m < self.slcons_tolerance
            return slcons[keep_indx].explode().reset_index()

        lndare_poly_gdf = gpd.GeoDataFrame(geometry=[land_poly], crs=self.wgs84)
        
        possible_slcons_idx = list(sindex_slcons.intersection(land_poly.bounds))
        possible_slcons = slcons_gdf.iloc[possible_slcons_idx]

        precise_slcons = possible_slcons.intersection(land_poly)
        linestring_indx = precise_slcons.geom_type == 'LineString'
        multilinestring_indx = precise_slcons.geom_type == 'MultiLineString'
        precise_slcons = precise_slcons[linestring_indx | multilinestring_indx]
        precise_slcons = precise_slcons[~precise_slcons.is_empty]

        # TODO: include polygon interiors
        if precise_slcons.shape[0] > 0:
            precise_slcons = remove_short_linestrings(precise_slcons)
            precise_slcons_geom = MultiLineString([g for g in precise_slcons[0]])
            slcons = gpd.GeoDataFrame(geometry=[precise_slcons_geom], crs=self.wgs84)
        
            lndare_exterior_modified = lndare_poly_gdf.exterior.difference(slcons)

            def splice_lndare_exterior(clipped_shoreline):
                verts = []
                for cs in clipped_shoreline.explode():
                    verts.extend(list(cs.coords))
                verts.append(verts[0])
                spliced_lndare_exterior = gpd.GeoDataFrame(geometry=[LineString(verts)],
                                                            crs=self.wgs84)
                return spliced_lndare_exterior

            if not lndare_exterior_modified.iloc[0].is_empty:
                spliced_lndare_exterior = splice_lndare_exterior(lndare_exterior_modified)
                return spliced_lndare_exterior.difference(poly.boundary).explode().reset_index()

            print('INTERIORS', '-' * 50)
            print(lndare_poly_gdf.interiors)
        else:
            return lndare_poly_gdf.boundary.difference(poly.boundary).explode().reset_index()

    def create_shoreline(self, objects, band):
        m_covr = objects['Polygon']['M_COVR']
        m_cscl = objects['Polygon']['M_CSCL']
        lndare = objects['Polygon']['LNDARE']
        lakare = objects['Polygon']['LAKARE']
        rivers = objects['Polygon']['RIVERS']

        slcons = objects['LineString']['SLCONS']
        slcons_df = pd.concat(slcons, ignore_index=True)
        slcons_gdf = gpd.GeoDataFrame(slcons_df.geometry, crs=self.wgs84)#.explode().reset_index()

        lndare_gdf = self.cookie_cut_inland_waters(lndare, lakare, rivers)

        print('intersecting LNDARE with M_COVR to get pre-reference...')
        sindex_slcons = slcons_gdf.sindex
        sindex_lndare = lndare_gdf.sindex
        shoreline_bits = []
        
        mcovr_df = pd.concat(m_covr, ignore_index=True)
        mcovr_gdf = gpd.GeoDataFrame(mcovr_df.geometry, crs=self.wgs84).explode().reset_index()
        enc_polys = self.union_mcscl(mcovr_gdf, m_cscl)
        num_enc_polys = enc_polys.shape[0]
        for i, poly in enumerate(enc_polys.geometry, 1):
            print('enc_poly {} of {}...'.format(i, num_enc_polys))
            possible_lndare_idx = list(sindex_lndare.intersection(poly.bounds))
            possible_lndares = lndare_gdf.iloc[possible_lndare_idx]

            #''' some lndare features overlap; union to avoid parts 
            #of lndare exteriors ending up in the shoreline layer'''
            possible_lndares = possible_lndares.unary_union  # possible_lndares is MultiPolygon?
            precise_lndares = possible_lndares.intersection(poly)

            if precise_lndares.geom_type == 'MultiPolygon':
                precise_lndares = [l for l in precise_lndares if l.geom_type == 'Polygon']
            elif precise_lndares.geom_type == 'Polygon':
                precise_lndares = [precise_lndares]

            try:
                num_landare = len(precise_lndares)
                for j, land_poly in enumerate(precise_lndares, 1):
                    #print('lndare {} of {} (enc_poly {} of {})...'.format(j, num_landare, i, num_enc_polys))
                    shoreline_bit = self.erase_slcons(poly, sindex_slcons, slcons_gdf, land_poly)
                    shoreline_bits.append(shoreline_bit)
            except Exception as e:
                print(e)
            
        self.save_shoreline(shoreline_bits, band)

    def gen_cusp_band_regions(self):
        bands = []
        for b in self.bands_to_process:
            layer= r'M_COVR_Polygon_Band{}'.format(b)
            gdf = gpd.read_file(str(self.enc_objects), layer=layer)
            bands.append(gdf)

        bands.reverse()
        current_coverage = bands[0]  # i.e., largest-scale band
        for b in bands[1:6]:  # start with largest-scale band and "work down"
            b_not_higher = gpd.overlay(b, current_coverage, how='difference').explode()
            b_not_higher = b_not_higher[b_not_higher.geometry.area > self.sliver_tolerance]
            df = pd.concat([current_coverage, b_not_higher], ignore_index=True)
            current_coverage = gpd.GeoDataFrame(df, geometry='geometry', crs=self.wgs84)
        
        layer = 'CUSP_band_regions'
        current_coverage.explode().to_file(self.cusp_ref_gpkg_path,
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


def main():
    enc_dir = Path(r'C:\ENCs\All_ENCs\ENC_ROOT')
    ref_dir = Path(r'Z:\ENC_Shoreline_Extractor')
    bands_to_process = [5]
    shorex = Shorex(enc_dir, ref_dir, bands_to_process)
    shorex.set_env_vars('enc_shorex')

    for band in shorex.bands_to_process:
        print('processing band {}...'.format(band))
        objects_to_extract = {'Point': [],
                              'LineString': ['COALNE', 'SLCONS'],
                              'Polygon': ['LNDARE', 'LNDRGN', 
                                          'M_COVR', 'M_CSCL',
                                          'RIVERS', 'LAKARE']}

        objects = {'Point': {}, 'LineString': {}, 'Polygon': {}}

        encs = list(shorex.enc_dir.rglob('US{}CA*.000'.format(band)))
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
        shorex.create_shoreline(objects, band)

    shorex.gen_cusp_band_regions()
    shorex.create_cusp_ref()


if __name__ == '__main__':
    main()
    