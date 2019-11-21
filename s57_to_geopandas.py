import os
import numpy as np
from pathlib import Path
import pandas as pd
import geopandas as gpd
import fiona
from shapely.geometry import Polygon, LineString, MultiLineString, MultiPolygon
from shapely.ops import polygonize, polygonize_full
from cartopy.geodesic import Geodesic
import matplotlib.pyplot as plt
import cProfile


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
        self.enc_objects = self.ref_dir / 'ENC_Objects_BAND.gpkg'
        self.ref_bands_gpkg_path = self.ref_dir / 'Reference_Bands.gpkg'
        self.cusp_ref_gpkg_path = self.ref_dir / 'CUSP_Reference.gpkg'
        self.sliver_tolerance = 0.0000005
        self.poly_tolerance = 0.0000001
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

    def export_mcovr_list(self, mcovr_list, band):
        df = pd.concat(mcovr_list)
        df['band'] = [band] * df.shape[0]
        gdf = gpd.GeoDataFrame(df, geometry='geometry',
                                crs=self.wgs84).explode()
                
        try:
            gpkg = str(self.enc_objects).replace('BAND', 'BAND{}'.format(band))
            layer = 'M_COVR_Polygon_Band{}'.format(band)
            print('exporting {} to {}...'.format(layer, gpkg))
            gdf.to_file(gpkg, layer=layer, driver='GPKG')
        except Exception as e:
            print(e, '({} - {})'.format(acronym, geom_type))

    def export_to_gpkg(self, objects, band):
        for geom_type, v in objects.items():
            for acronym, objects_list in v.items():
                if acronym == 'M_COVR':
                    if not all(o is None for o in objects_list):  # if not all None
                        df = pd.concat(objects_list)
                        df['band'] = [band] * df.shape[0]
                        gdf = gpd.GeoDataFrame(df, geometry='geometry',
                                                crs=self.wgs84).explode()
                
                        try:
                            gpkg = str(self.enc_objects).replace('BAND', 'BAND{}'.format(band))
                            layer = '{}_{}_Band{}'.format(acronym, geom_type, band)
                            print('exporting {} to {}...'.format(layer, gpkg))
                            gdf.to_file(gpkg, layer=layer, driver='GPKG')
                        except Exception as e:
                            print(e, '({} - {})'.format(acronym, geom_type))

    def cookie_cut_inland_waters(self, lndare, lakare, rivers):
        lndare_df = pd.concat(lndare, ignore_index=True)
        lndare_gdf = gpd.GeoDataFrame(lndare_df.geometry, crs=self.wgs84).explode().reset_index()

        def difference_objects(lndare_gdf, objects):
            """differe LNDARE with other objects to get polygon the exterior
            and interior(s) of which are defined to be the shoreline"""

            objects_df = pd.concat(objects, ignore_index=True)
            objects_gdf = gpd.GeoDataFrame(objects_df.geometry, crs=self.wgs84).explode().reset_index()

            if not objects_gdf.empty:
                cols = lndare_gdf.columns
                lndare_gdf = lndare_gdf.drop([c for c in cols if c is not 'geometry'], axis=1)
                lndare_gdf = gpd.overlay(lndare_gdf, objects_gdf, how='difference').explode().reset_index()            
            return lndare_gdf
        
        try:
            lndare_gdf = difference_objects(lndare_gdf, lakare)
        except Exception as e:
            print(e)

        try:
            lndare_gdf = difference_objects(lndare_gdf, rivers)
        except Exception as e:
            print(e)

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
        gdf = gpd.GeoDataFrame(geometry=df.geometry, crs=self.wgs84)
        gdf = gdf[gdf.geom_type == 'LineString']
        gdf.to_file(str(self.ref_bands_gpkg_path), 
                    layer='CUSP_Reference_Band{}'.format(band),
                    driver='GPKG')

    def create_shoreline(self, objects, band):
        print('assembling ENC objects...')
        m_covr = objects['Polygon']['M_COVR']
        m_cscl = objects['Polygon']['M_CSCL']
        lndare = objects['Polygon']['LNDARE']
        lakare = objects['Polygon']['LAKARE']
        rivers = objects['Polygon']['RIVERS']
        slcons = objects['LineString']['SLCONS']
        
        slcons_df = pd.concat(slcons, ignore_index=True)
        slcons_gdf = gpd.GeoDataFrame(slcons_df.geometry, crs=self.wgs84)

        print('cookie-cutting out LAKARE and RIVERS...')
        lndare_gdf = self.cookie_cut_inland_waters(lndare, lakare, rivers)

        print('intersecting LNDARE with M_COVR to get pre-reference...')
        sindex_slcons = slcons_gdf.sindex
        sindex_lndare = lndare_gdf.sindex
        shoreline = []
        
        mcovr_df = pd.concat(m_covr, ignore_index=True)
        mcovr_gdf = gpd.GeoDataFrame(mcovr_df.geometry, crs=self.wgs84).explode().reset_index()
        enc_polys = self.union_mcscl(mcovr_gdf, m_cscl)
        num_enc_polys = enc_polys.shape[0]
        for i, poly in enumerate(enc_polys.geometry, 1):

            def get_land_polys(lndares):
                #''' some lndare features overlap; union to avoid non-shoreline parts 
                #of lndare exteriors ending up in the shoreline layer'''
                lndares = gpd.GeoSeries(lndares)
                lndares = lndares[lndares.geom_type == 'Polygon']
                #lndares = gpd.GeoSeries(lndares).unary_union
                #lndares = gpd.GeoSeries(lndares).explode()
                #lndares = lndares[lndares.geom_type == 'Polygon'].reset_index()[0]
                lndares = gpd.GeoDataFrame(geometry=lndares, crs=self.wgs84)
                return lndares.geometry

            def erase_enc_poly(land_polys, poly):
                if land_polys:
                    land_multipoly_boundary = MultiPolygon(land_polys).boundary
                    shoreline_geom = land_multipoly_boundary.difference(poly.boundary)
                else:
                    shoreline_geom = []
                return gpd.GeoSeries(shoreline_geom, crs=self.wgs84).explode()

            def erase_slcons(land_poly, sindex_slcons, slcons_gdf):

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

                    # keep "large piers" (as defined by threshold)
                    lengths_m = np.asarray([linestring_distance(g) for g in slcons])
                    keep_indx = lengths_m < self.slcons_tolerance
                    return slcons[keep_indx]

                def splice_lndare_exterior(clipped_shoreline):
                    verts = []
                    for cs in clipped_shoreline:
                        verts.extend(cs.coords)
                    verts.append(verts[0])
                    return Polygon(verts)

                # account for Lake-Tahoe--type case (with interior ring)
                lndare_rings = [land_poly.exterior] + list(land_poly.interiors)
                shoreline_ring_bits = []
                for land_ring in lndare_rings:
                    land_ring_poly = Polygon(land_ring)
                    lndare_ring_gdf = gpd.GeoDataFrame(geometry=[land_ring_poly], crs=self.wgs84)
        
                    possible_slcons_idx = list(sindex_slcons.intersection(land_ring_poly.bounds))
                    possible_slcons = slcons_gdf.iloc[possible_slcons_idx]

                    precise_slcons = possible_slcons.intersection(land_ring_poly)
                    precise_slcons = precise_slcons.unary_union
                    geo = gpd.GeoSeries(precise_slcons).explode()
                    lines_geom = geo[geo.geom_type == 'LineString'].reset_index()[0]
                    precise_slcons = gpd.GeoDataFrame(geometry=lines_geom, crs=self.wgs84)

                    if not precise_slcons.empty:
                        precise_slcons = remove_short_linestrings(precise_slcons.geometry)
                        precise_slcons_multiline = MultiLineString(list(precise_slcons))
                        slcons = gpd.GeoDataFrame(geometry=[precise_slcons_multiline], crs=self.wgs84)
                        lndare_bound_mod = lndare_ring_gdf.boundary.difference(slcons)
                        lndare_bound_mod = lndare_bound_mod.explode().reset_index()[0]

                        if not lndare_bound_mod.empty:
                            spliced_shoreline_poly = splice_lndare_exterior(lndare_bound_mod)
                            shoreline_ring_bits.append(spliced_shoreline_poly)

                    else:
                        # precise_slcons will be empty if only ends of slcons touch lndare
                        shoreline_ring_bits.append(Polygon(lndare_ring_gdf.exterior.iloc[0]))

                return shoreline_ring_bits

            def get_land_linestrings_1(land_polys):
                shoreline_bits = []
                for l in land_polys:
                    slcon_erased = erase_slcons(l, sindex_slcons, slcons_gdf)
                    shoreline_bits.extend(slcon_erased)
                return erase_enc_poly(shoreline_bits, poly)

            def get_land_linestrings_2(land_polys):
                shoreline_bits = []
                for l in land_polys:
                    slcon_erased = erase_slcons(l, sindex_slcons, slcons_gdf)
                    shoreline_bits.extend(slcon_erased)
                return gpd.GeoSeries(shoreline_bits, crs=self.wgs84)

            def get_land_linestrings_3(land_polys):
                return erase_enc_poly(list(land_polys), poly)

            def get_land_linestrings_4(land_polys):
                return land_polys.boundary.explode()

            print('band {} enc_poly {} of {}...'.format(band, i, num_enc_polys))
            possible_lndare_idx = list(sindex_lndare.intersection(poly.bounds))
            possible_lndares = lndare_gdf.iloc[possible_lndare_idx]
            prelim_land_polys = possible_lndares.intersection(poly)
                
            possible_poly_slcons_idx = list(sindex_slcons.intersection(poly.bounds))
            possible_poly_slcons = slcons_gdf.iloc[possible_poly_slcons_idx]
            poly_slcons = possible_poly_slcons.intersection(poly)

            land_polys = get_land_polys(prelim_land_polys)
            indx_poly = land_polys.intersects(poly.boundary)
            indx_slcons = land_polys.intersects(MultiLineString(list(poly_slcons)))
            
            s1_indx = indx_poly & indx_slcons
            s2_indx = ~indx_poly & indx_slcons
            s3_indx = indx_poly & ~indx_slcons
            s4_indx = ~indx_poly & ~indx_slcons

            shoreline_1 = get_land_linestrings_1(land_polys[s1_indx])
            shoreline_2 = get_land_linestrings_2(land_polys[s2_indx])
            shoreline_3 = get_land_linestrings_3(land_polys[s3_indx])
            shoreline_4 = get_land_linestrings_4(land_polys[s4_indx])

            shorelines = [shoreline_1, shoreline_2, shoreline_3, shoreline_4]

            for s in shorelines:
                if not s.empty:
                    shoreline.append(s)
            
        return shoreline

    def gen_cusp_band_regions(self):
        bands = []
        for b in self.bands_to_process:
            try:
                layer= r'M_COVR_Polygon_Band{}'.format(b)
                gpkg = str(self.enc_objects).replace('BAND', 'BAND{}'.format(b))
                gdf = gpd.read_file(gpkg, layer=layer)
                bands.append(gdf)
            except Exception as e:
                print(e)

        bands.reverse()
        current_coverage = bands[0]  # i.e., largest-scale band
        for i, b in enumerate(bands[1:], 1):  # start with largest-scale band and "work down"
            b_not_higher = gpd.overlay(b, current_coverage, how='difference').explode().reset_index()
            b_not_higher = b_not_higher[b_not_higher.geometry.area > self.sliver_tolerance]
            b_not_higher.geometry = b_not_higher.geometry.simplify(tolerance=self.poly_tolerance, 
                                                                   preserve_topology=True)
            df = pd.concat([current_coverage, b_not_higher], ignore_index=True)
            current_coverage = gpd.GeoDataFrame(df, geometry='geometry', crs=self.wgs84)
        
        layer = 'CUSP_band_regions'
        current_coverage.to_file(self.cusp_ref_gpkg_path,
                                 layer=layer, driver='GPKG')
        
    def create_cusp_ref(self):
        band_regions = gpd.read_file(str(self.cusp_ref_gpkg_path), 
                                     layer='CUSP_band_regions')
        
        for b in self.bands_to_process:
            try:
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
            except Exception as e:
                print(e)


def main():

    enc_dir = Path(r'C:\ENCs\All_ENCs\ENC_ROOT')
    ref_dir = Path(r'Z:\ENC_Shoreline_Extractor')
    bands_to_process = [4]
    shorex = Shorex(enc_dir, ref_dir, bands_to_process)
    shorex.set_env_vars('enc_shorex')

    for band in shorex.bands_to_process:
        print('processing band {}...'.format(band))
        objects_to_extract = {'LineString': ['COALNE', 'SLCONS'],
                              'Polygon': ['LNDARE', 'LNDRGN', 
                                          'M_COVR', 'M_CSCL',
                                          'RIVERS', 'LAKARE']}

        encs = list(shorex.enc_dir.rglob('US{}*.000'.format(band)))
        encs = [e for e in encs if e.stem in ['US4MN1FH']]
        num_encs = len(encs)
        
        mcovr_list = []
        band_shoreline = []

        for i, enc in enumerate(encs, 1):
            objects = {'LineString': {}, 'Polygon': {}}
            enc = Enc(enc)
            print('-' * 50)
            print('{} ({} of {})'.format(enc, i, num_encs))

            for geom_type in objects_to_extract.keys():
                for acronym in objects_to_extract[geom_type]:
                    enc_objs = enc.get_objects(acronym, geom_type)

                    if acronym not in objects[geom_type].keys():
                        objects[geom_type].update({acronym: []})

                    objects[geom_type][acronym].append(enc_objs)

            try:
                mcovr_list.extend(objects['Polygon']['M_COVR'])
                #cProfile.runctx('shorex.create_shoreline(objects, band)', globals(), locals(), sort="cumtime")
                band_shoreline.extend(shorex.create_shoreline(objects, band))
            except Exception as e:
                print(e)

        #try:
        shorex.export_mcovr_list(mcovr_list, band)
        shorex.save_shoreline(band_shoreline, band)
        #except Exception as e:
        #    print(e)

    shorex.gen_cusp_band_regions()
    shorex.create_cusp_ref()


if __name__ == '__main__':
    main()
    