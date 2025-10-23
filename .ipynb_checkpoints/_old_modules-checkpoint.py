import numpy as np
import pandas as pd
import ast
from odc.stac import stac_load


# standard values
chunks = {"x": 1024, "y": 1024, "time": 1}


def get_alias_band(product,alias):
    # create list for each row and save as new row entry
    # overview file of datasets
    df = pd.read_csv("measurements.csv")
    dfs = df.loc[df['product']==product]
    # measurements_aliases = ['QA_PIXEL', 'blue', 'green', 'red', 'nir', 'swir_1', 'swir_2']
    bands_to_keep = []
    bands_alias = []
    for r in range(len(dfs.index)):
        # r['measurement']
        _alias = dfs.iloc[r,:]['aliases']
        _alias = ast.literal_eval(_alias)
        _meas = dfs.iloc[r,:]['measurement']
        _alias.append(_meas)
        # check if any of the "measurements_aliases" is inside _alias
        _match = set(alias).intersection(_alias)
        if len(_match) > 0:
        # if any(map(lambda v: v in measurements_aliases, _alias)):
            # actual band name
            bands_to_keep.append(_meas)
            # keep also the alias name, at least return this as info
            bands_alias.append(_match.pop())
    
    print(f"Product: {product}")
    # print(f"Bands found in product: {bands_alias}")
    # print(f"Band names real       : {bands_to_keep}")
    
    
    df_out = pd.DataFrame({'measurement': bands_to_keep, 'alias': bands_alias})
    print(df_out)

    return bands_to_keep, bands_alias

def load_product_ts(catalog=None, product=None, **kwargs):

    assert 'measurements' in kwargs.keys(), \
           "\n! <measurements> is required !"
    assert 'QA_PIXEL' in kwargs['measurements'] or 'SLC' in kwargs['measurements'], \
           "\n! <measurements> must contains 'QA_PIXEL' or 'SLC' band !"
    assert 'resolution' in kwargs.keys() and 'output_crs' in kwargs.keys(), \
       "\n! <resolution> and <output_crs> are required !"
    
    # Extract optional parameters with default values if not provided
    time = kwargs.get('time')
    longitude = kwargs.get('longitude')
    latitude = kwargs.get('latitude')
    measurements = kwargs.get('measurements')
    output_crs = kwargs.get('output_crs')
    resolution = kwargs.get('resolution', [None, None])
    rename = kwargs.get('rename', False)
    alias_names = kwargs.get('alias_names')
    chunks = kwargs.get('chunks', None)

    # Build query with conditional checks for optional params
    datetime_range = f"{time[0]}/{time[1]}" if time else None
    bbox = (longitude[0], latitude[0], longitude[1], latitude[1]) if longitude and latitude else None

    query = catalog.search(
        collections=[product],
        datetime=datetime_range,
        limit=100,
        bbox=bbox
    )
    items = list(query.items())

    # Load identified items
    lazy_ds = stac_load(
        items,
        lon=longitude,
        lat=latitude,
        bands=measurements,
        crs=output_crs,
        resolution=resolution[1] if resolution else None,
        chunks=chunks,
    )

    # Rename variables if needed
    if rename and alias_names:
        _olds = list(lazy_ds.data_vars.keys())
        for r in range(len(alias_names)):
            if _olds[r] != alias_names[r]:  # only non-identical, otherwise error
                lazy_ds = lazy_ds.rename({_olds[r]: alias_names[r]})

    return lazy_ds
    
# def load_product_ts(catalog=None, product=None, time=None, longitude=None, latitude=None, measurements=None, output_crs=None,
#                     resolution=None, rename=False, alias_names = None):
#     query = catalog.search(
#         collections=[product],
#         datetime=f"{time[0]}/{time[1]}",
#         limit=100,
#         bbox=(longitude[0], latitude[0],
#               longitude[1], latitude[1])
#     )
#     items = list(query.items())
    
#     # load identified items
#     lazy_ds = stac_load(
#         items,
#         lon=longitude,
#         lat=latitude,
#         bands=measurements,
#         crs=output_crs,
#         resolution=resolution[1],
#         chunks=chunks,
#     )
#     if rename:
#         _olds = list(lazy_ds.data_vars.keys())
#         for r in range(len(alias_names)):
#             if _olds[r] != alias_names[r]:  # only non-identical, otherwise error
#                 lazy_ds = lazy_ds.rename({_olds[r]:alias_names[r]})
#     return(lazy_ds)


def load_product(catalog , product, longitude, latitude, measurements, output_crs, resolution, rename=False, alias_names = None):
    query = catalog.search(
        collections=[product],
        # datetime=f"{time[0]}/{time[1]}",
        limit=100,
        bbox=(longitude[0], latitude[0],
              longitude[1], latitude[1])
    )
    items = list(query.items())
    
    # load identified items
    lazy_ds = stac_load(
        items,
        lon=longitude,
        lat=latitude,
        bands=measurements,
        crs=output_crs,
        resolution=resolution[1],
        chunks=chunks,
    )
    
    lazy_ds['QA_PIXEL'].attrs['units'] = 'bit_index'
    lazy_ds['QA_PIXEL'].attrs['flags_definition'] = []
    lazy_ds['QA_PIXEL'] = lazy_ds['QA_PIXEL'].astype(np.int64)
    if rename:
        _olds = list(lazy_ds.data_vars.keys())
        for r in range(len(alias_names)):
            if _olds[r] != alias_names[r]:  # only non-identical, otherwise error
                lazy_ds = lazy_ds.rename({_olds[r]:alias_names[r]})
    return(lazy_ds)


def overview_data(product):
    # create list for each row and save as new row entry
    # overview file of datasets
    df = pd.read_csv("measurements.csv")
    dfs = df.loc[df['product']==product]
    # measurements_aliases = ['QA_PIXEL', 'blue', 'green', 'red', 'nir', 'swir_1', 'swir_2']
    bands_to_keep = []
    bands_alias = []
    for r in range(len(dfs.index)):
        # r['measurement']
        _alias = dfs.iloc[r,:]['aliases']
        _alias = ast.literal_eval(_alias)
        _meas = dfs.iloc[r,:]['measurement']
        _alias.append(_meas)
        # check if any of the "measurements_aliases" is inside _alias
        _match = set(measurements_aliases).intersection(_alias)
        if len(_match) > 0:
        # if any(map(lambda v: v in measurements_aliases, _alias)):
            # actual band name
            bands_to_keep.append(_meas)
            # keep also the alias name, at least return this as info
            bands_alias.append(_match.pop())
    
    print(f"Product: {product}")
    # print(f"Bands found in product: {bands_alias}")
    # print(f"Band names real       : {bands_to_keep}")
    
    
    df_out = pd.DataFrame({'measurement': bands_to_keep, 'alias': bands_alias})
    print(df_out)

    return bands_to_keep, bands_alias


def create_scl_clean_mask(scl, valid_cats = [4, 5, 6, 7, 11]):
    return xr.apply_ufunc(np.isin, scl, valid_cats, dask='allowed')



# def create_qa_pixel_clean_mask(qa_pixel, valid_cats = [4, 5, 6, 7, 11]):
#     return xr.apply_ufunc(np.isin, scl, valid_cats, dask='allowed')
    
def lsc2_qa_clean(da_qa, valid_bits = None, invalid_bits = None,
                  rm_high_conf = True, rm_mid_conf = True, force_snow = False):
    """
    Description:
      Creates a clean mask of a Landsat Collection 2 dataset using a given band
      and a list of valid bits.
      
      The function was written to supports loading the standard SDC products:
      - landsat_tm_c2_l2
      - landsat_etm_c2_l2
      - landsat_ot_c2_l2
      And tested with the following measurements:
      - QA_PIXEL
      - SR_CLOUD_QA (Landsat tm and etm only)
      
      Default <valid_bits> depends of the selected quality measurement:
      - QA_PIXEL: [5, 6, 7] (Snow, Clear, Water)
      - SR_CLOUD_QA: [0, 4, 5] (dark_dense_vegetation, snow, water)
      Your can use `show_bit_flags(<product>, <measurement>)`
      
      You can visualize the flag attributed to bit(s) using the `show_bit_flags` function
      Or have a look at official documentation:
      - https://www.usgs.gov/media/files/landsat-4-7-collection-2-level-2-science-product-guide
      - https://www.usgs.gov/media/files/landsat-8-collection-2-level-2-science-product-guide
      
    Args:
      da_qa:        quality band of a Landsat Collection 2 xarray.DataArray
      valid_bits:   (OPTIONAL) list of valid bits (see description above)
      invalid_bits: (OPTIONAL) list of bits considered an invalid (as often a pixel can
                    be flagged with several categories
      rm_high_conf: (OPTIONAL, only used with QA_PIXEL measurement and True by default)
                    Removes high confidence cloud, cloud shadow and cirrus to prioritize
                    snow
      rm_mid_conf:  (OPTIONAL, only used with QA_PIXEL measurement and True by default)
                    Removes mid confidence cloud to prioritize snow
      force_snow:   (OPTIONAL, False by default) consider any pixel with snow high confidence
                    as data even if it also has a high (or medium) confidence level for cloud,
                    cloud shadow or cirrus
      
    Output:
      clean_mask (boolean numpy array)
    """
    
    # Check submitted input
    assert isinstance(da_qa, xr.DataArray), \
            "SCRIPT INTERRUPTED: da_qa should be an xarray.DataArray"
    assert da_qa.units == 'bit_index', \
        f"SCRIPT INTERRUPTED: no flag_definition found for measurement '{da_qa.name}'"
    
    if da_qa.name != 'QA_PIXEL' and rm_high_conf is not None \
       and rm_mid_conf is not None:
        warnings.warn("<rm_high_conf> and <rm_mid_conf> not used")
    
    # List and count all da_qa unique values
    ok_qas, da_cnt = unik_count(da_qa.values)
    ok_qas = list(ok_qas)
    # Return bit encoding    
    bit_len = bit_length(da_qa)
    
    if da_qa.name == 'QA_PIXEL':
        # List snow forced qas
        if force_snow:
            snow_qas = []
            for v in sorted(ok_qas):
                b = str(bin(v))[2:].zfill(bit_len)[::-1]
                if b[5] == '1' and b[12:14] == '11':
                    snow_qas.append(v)
        
        if valid_bits is None:
            valid_bits = [5, 6, 7]
        
        # Get rid of high confidence cloud, cloud shadow and cirrus
        if rm_high_conf:
            for v in sorted(ok_qas):
                b = str(bin(v))[2:].zfill(bit_len)[::-1]
                if b[8:10] == '11' or b[10:12] == '11' or b[14:16] == '11':
                    ok_qas.remove(v)
        
        # Get rid of mid confidence cloud
        if rm_mid_conf:
            for v in sorted(ok_qas):
                b = str(bin(v))[2:].zfill(bit_len)[::-1]
                if b[8:10] == '01':
                    ok_qas.remove(v)
            
        # Append snow forced qas
        if force_snow:
            ok_qas = list(set(ok_qas + snow_qas))            
            
    elif da_qa.name == 'SR_CLOUD_QA':
        if valid_bits is None:
            valid_bits = [0, 4, 5]
    else:
        if valid_bits is None:
            sys.exit('SCRIPT INTERRUPTED: you need to define <valid_bits>, ' \
                     'use `show_bit_flags(<product>, <measurement>)` to make sure ' \
                     'you define proper bits measurement and values')
    
    # Keep only valid_bits
    data_qas = []
    for v in sorted(ok_qas):
        b = str(bin(v))[2:].zfill(bit_len)[::-1]
        for c in valid_bits:
            if b[c] == '1':
                data_qas.append(v)
    data_qas = sorted(list(set(data_qas)))
    
    # Remove invalid_bits
    if invalid_bits is not None:
        for v in sorted(data_qas):
            b = str(bin(v))[2:].zfill(bit_len)[::-1]
            for c in invalid_bits:
                if b[c] == '1':
                    data_qas.remove(v)
                    break
        data_qas = sorted(list(set(data_qas)))      
    
    bool_arr = xr.apply_ufunc(np.isin, da_qa, data_qas, dask = 'allowed').values
    return xr.DataArray(bool_arr, coords=da_qa.coords, dims=da_qa.dims, attrs=da_qa.attrs)

def decod_bit(v, da_qa):
    bits_dict = {'0': 'No',
                 '1': 'Yes',
                '00': 'No confidence level set',
                '01': 'Mid/Reserved',
                '10': 'Low',
                '11': 'High'}
    df_decoded = pd.DataFrame(columns=['bits', 'zo', 'flag', 'Value', 'tmp'])
    bit_len = bit_length(da_qa)
    b = str(bin(v))[2:].zfill(bit_len)[::-1]
    for k, itm in da_qa.flags_definition.items():        
        if isinstance(itm['bits'], list):
            zo = b[itm['bits'][0]:itm['bits'][1] + 1]
            tv = itm['bits'][0]
        else:
            zo = b[itm['bits']]
            tv = itm['bits']
        new_row = pd.DataFrame([{'bits': itm['bits'], 'zo':zo, 'flag':k, 'Value':bits_dict[zo], 'tmp': tv}])
        df_decoded = pd.concat([df_decoded, new_row], axis=0, ignore_index=True)
    return df_decoded.sort_values('tmp').drop(['tmp'], axis = 1).reset_index(drop=True)


def lsc2_ds_scale(ds_fun, p, vr = True):
    # ! Quality band needs to be 'QA_PIXEL
    # for more details refer to:
    # table 5-1 of https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/LSDS-1618_Landsat-4-7_C2-L2-ScienceProductGuide-v4.pdf
    # table 6-1 of https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/LSDS-1619_Landsat-8-9-C2-L2-ScienceProductGuide-v4.pdf
    
    data=[
          ['landsat_tm_c2_l2', f"{os.path.dirname(os.path.abspath(__file__))}/data/lsc2_scale_table_457.csv"],
          ['landsat_etm_c2_l2', f"{os.path.dirname(os.path.abspath(__file__))}/data/lsc2_scale_table_457.csv"],
          ['landsat_ot_c2_l2', f"{os.path.dirname(os.path.abspath(__file__))}/data/lsc2_scale_table_89.csv"],
         ]
    
    df = pd.DataFrame(data, columns=['product', 'table'])
    df = pd.read_csv(list(df.loc[df['product'] == p]['table'])[0])
    df['all_names'] = df['Band_Designation'] + ',' + df['Aliases']
    df['all_names'] = df.all_names.str.split(',')
    
    for v in ds_fun.var():
        if v != 'QA_PIXEL':
            v_row = df[df['all_names'].apply(lambda x: v in x)].iloc[0]

            # set nodata
            if not pd.isna(v_row.Fill_Value):
                ds_fun[v] = ds_fun[v].where(ds_fun[v] != v_row.Fill_Value)

            # restrict to valid range if requested
            if vr and not pd.isna(v_row.Valid_Range_min) and not pd.isna(v_row.Valid_Range_max):
                ds_fun[v] = ds_fun[v].where(np.all([ds_fun[v] >= v_row.Valid_Range_min,
                                                    ds_fun[v] <= v_row.Valid_Range_max], axis = 0))

            # Apply Multiplicative_Scale_Factor
            if not pd.isna(v_row.Multiplicative_Scale_Factor):
                ds_fun[v] = ds_fun[v] * v_row.Multiplicative_Scale_Factor
            # Apply Additive_Offset
            if not pd.isna(v_row.Additive_Offset):
                ds_fun[v] = ds_fun[v] + v_row.Additive_Offset
    
    return ds_fun
