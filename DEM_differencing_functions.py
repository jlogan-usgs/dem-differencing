'''
Functions for DEM differencing and estimating volumetric uncertainty.

'''
import numpy as np
import osgeo
import rasterio
from rasterio.coords import BoundingBox
from rasterio.mask import mask

import pandas as pd
import geopandas as gpd
import math

from pathlib import Path
import matplotlib.pyplot as plt

import os
import shutil

from osgeo import gdal, gdalconst
gdal.UseExceptions()

#For zonal stats
from rasterio.features import rasterize
from rasterstats import zonal_stats

# For spatial semivariogram
from pyinterpolate import read_txt
from pyinterpolate import build_experimental_variogram
from pyinterpolate import TheoreticalVariogram, build_theoretical_variogram

# For spatial semivariogram, using gstat
import skgstat as skg
from scipy.optimize import curve_fit
from skgstat import models

# For resample
from rasterio.enums import Resampling

# For raster outline to polygon shapefile
from rasterio.features import dataset_features

#For create polygon from points
from shapely.geometry import Polygon
from pyproj import CRS, Transformer

#-------------------------------------------------------
def DEM_difference(dem1, dem2, dod):
    '''
    Use rasterio to subtract dem1 from dem2 and write output to a new geotiff (dod)
    
    dem1: Path object to DEM1
    dem2: Path object to DEM2
    dod: Path object to output DoD
    '''
    #DEM difference
    #with ChatGPT assitance
    # Open the input raster files using rasterio's context manager
    with rasterio.open(dem1) as dataset1, rasterio.open(dem2) as dataset2:
        # Read the raster data as NumPy arrays
        raster_data1 = dataset1.read(1)
        raster_data2 = dataset2.read(1)
    
        # Retrieve the affine transformation matrices
        transform1 = dataset1.transform
        transform2 = dataset2.transform
    
        # Determine the common extent
        bounds1 = dataset1.bounds
        bounds2 = dataset2.bounds
        intersection_bounds = BoundingBox(
            max(bounds1.left, bounds2.left),
            max(bounds1.bottom, bounds2.bottom),
            min(bounds1.right, bounds2.right),
            min(bounds1.top, bounds2.top)
            )
    
       # Update the transformation matrices to the common extent
        window1 = dataset1.window(*intersection_bounds)
        window2 = dataset2.window(*intersection_bounds)
        transform1 = rasterio.windows.transform(window1, transform1)
        transform2 = rasterio.windows.transform(window2, transform2)
    
        # Read the subset of raster data within the common extent
        raster_data1 = dataset1.read(1, window=window1)
        raster_data2 = dataset2.read(1, window=window2)
    
        # Compute the difference between the two rasters
        difference = raster_data2 - raster_data1
        
        # Set nodata values in the difference array
        nodata_mask = np.logical_or(raster_data1 == dataset1.nodata, raster_data2 == dataset2.nodata)
        difference[nodata_mask] = dataset1.nodata
    
        # Retrieve the metadata from one of the input rasters
        metadata = dataset1.meta
    
        # Update the metadata to reflect the difference raster
        metadata.update(
            count=1,
            height=difference.shape[0],
            width=difference.shape[1],
            transform=transform1,
            nodata=dataset1.nodata,
            compress='lzw'
        )
    
        # Create a new raster file for the difference
        output_file = dod
        with rasterio.open(str(output_file), 'w', **metadata) as output_dataset:
            output_dataset.write(difference, 1)
    
    print(f"****************\nRaster difference complete.\n    Performed the following operation:\n        {dem2.name} - {dem1.name} = {dod.name}\n\nOutput written to {dod}\n****************\n")

def weighted_average(dataframe, value, weight):
    '''
    Calculate weighted average of a dataframe column, with vals and weights in columns
    ''' 
    val = dataframe[value]
    wt = dataframe[weight]
    return (val * wt).sum() / wt.sum()

def clip_raster_with_polygon(unclipped_raster, clip_shapefile, out_raster):
    '''
    Use rasterio to clip a raster using a shapefile
    
    unclipped_raster: Path object to raster
    clip_shapefile: Path object to clipping shapefile
    out_raster: Path object to output clipped raster (geotiff, .tif)
    '''
    
    # Read the shapefile using geopandas
    shapefile = gpd.read_file(clip_shapefile)
    
    # Read the raster using rasterio
    with rasterio.open(unclipped_raster) as src:
        # Clip the raster using the shapefile geometry
        clipped_raster, clipped_transform = mask(src, shapefile.geometry, crop=True)
        clipped_meta = src.meta
    
    # Update the metadata of the clipped raster
    clipped_meta.update({
        'height': clipped_raster.shape[1],
        'width': clipped_raster.shape[2],
        'transform': clipped_transform,
        'compress':'lzw'
    })
    
    # Write the clipped raster to the output path
    with rasterio.open(out_raster, 'w', **clipped_meta) as dst:
        dst.write(clipped_raster)
        
    print(f'****************\nRaster clipping complete.\n    Clipped raster saved to: {out_raster}\n****************\n')
    
def mask_raster_with_polygon(unmasked_raster, mask_shapefile, out_raster):
    '''
    Use rasterio to mask a raster using a shapefile
    
    unmasked_raster: Path object to raster
    mask_shapefile: Path object to masking shapefile
    out_raster: Path object to output clipped raster (geotiff, .tif)
    '''
    
    # Read the shapefile using geopandas
    shapefile = gpd.read_file(mask_shapefile)
    
    # Read the raster using rasterio
    with rasterio.open(unmasked_raster) as src:
        # Mask the raster using the shapefile geometry (invert=True)
        clipped_raster, clipped_transform = mask(src, shapefile.geometry, invert=True)
        clipped_meta = src.meta
    
    # Update the metadata of the clipped raster
    clipped_meta.update({
        'height': clipped_raster.shape[1],
        'width': clipped_raster.shape[2],
        'transform': clipped_transform,
        'compress':'lzw'
    })
    
    # Write the clipped raster to the output path
    with rasterio.open(out_raster, 'w', **clipped_meta) as dst:
        dst.write(clipped_raster)
        
    print(f'****************\nRaster masking complete.\n    Masked raster saved to: {out_raster}\n****************\n')
    
def mask_raster_with_raster(inras, maskras, outras):
    '''
    Mask a raster (remove cells) using a mask raster.  Must have same crs and cells should be aligned.
    
    inras: Path object to raster to be masked
    maskras: Path object to mask raster
    outras: Path object to output masked raster (geotiff, .tif)
    '''
    #ensure that inputs are path objects
    inras = Path(inras)
    maskras = Path(maskras)
    
    #make temporary directory for shapefile
    tempdir = Path(maskras.parent,'tempshp')
    tempdir.mkdir(exist_ok=True)

    #make temporary masking polygon
    maskpolypath = Path(maskras.parent,'tempshp',maskras.stem + '_boundary_temp.shp')
    raster_outline_to_polyshp(maskras, maskpolypath)
    print('... masking shapefile will be removed after operation\n')

    #mask raster with temporary polygon
    mask_raster_with_polygon(inras, maskpolypath, outras)

    #delete temporary polygon and dir
    shutil.rmtree(tempdir)

def create_polygon_from_points(points_file, output_shapefile, x_col='x', y_col='y', in_crs_epsg='6339', out_crs_epsg='6339'):
    '''
    Creates a polygon shapefile from a csv of points with x and y columns
    
    points_file: Path object to csv with points and header with col labels
    output_shapefile: Path object to clipping shapefile
    x_col: str name of x column in header
    y_col: str name of y column in header
    in_crs_epsg: EPSG of CRS for input
    out_crs_epsg: EPSG of CRS for output shp
    
    '''
    # Read the points from the CSV file using pandas
    df = pd.read_csv(points_file)

    # Extract x and y coordinates
    points = list(zip(df[x_col], df[y_col]))

    # Create a Polygon object from the points
    polygon = Polygon(points)

    # Create a GeoDataFrame with the polygon, set crs
    data = gpd.GeoDataFrame(geometry=[polygon], crs=in_crs_epsg)
    
    # Transform the polygon to output crs (default=EPSG:6339)
    data = data.to_crs(out_crs_epsg)

    # Save the GeoDataFrame as a shapefile 
    data.to_file(output_shapefile, driver='ESRI Shapefile')
    

def raster_outline_to_polyshp(inras, outpolyshp):
    '''
    Creates a polygon shapefile of outline of raster, excluding no_data.
    from https://gis.stackexchange.com/questions/429662/extracting-raster-outline-to-a-vector-geometry

    inras: Path object to raster
    outpolyshp: Path object to output shapefile
    '''
    
    with rasterio.open(inras) as ds:
        gdf = gpd.GeoDataFrame.from_features(dataset_features(ds, bidx=1, as_mask=True, geographic=False, band=False))
        gdf = gdf.set_crs(ds.crs)
        gdf.to_file(outpolyshp)

    print(f'****************\nPolygon shapefile saved to {outpolyshp}\n****************\n')
    

def calculate_volume(dod):
    '''
    Report volume, gets cell size from raster metadata. Returns dataframe with values.

    dod: Path object to DoD
    '''
    #Calc sum
    with rasterio.open(dod) as src:
        # Read the raster data as NumPy arrays
        masked_raster = src.read(1, masked=True)
        #get cell size x, y
        cellsize_x, cellsize_y  = src.res
        
        #calc net vol, area    
        vol_net = np.sum(masked_raster) * cellsize_x * cellsize_y
        area_net = masked_raster.count() * cellsize_x * cellsize_y

        #calc deposition (where dod is positive) and area
        vol_pos = np.sum(masked_raster[masked_raster > 0]) * cellsize_x * cellsize_y
        area_pos = masked_raster[masked_raster > 0].count() * cellsize_x * cellsize_y

        #calc erosion (where dod is negative) and area
        vol_neg = np.sum(masked_raster[masked_raster < 0]) * cellsize_x * cellsize_y
        area_neg = masked_raster[masked_raster < 0].count() * cellsize_x * cellsize_y
       
        #report
        # print(f'****************\nNet volume change\n****************\n    Net volume: {vol_net} m^3\n    Area: {area_net} m^2\n    Number of cells: {masked_raster.count()}\n    Cell size X: {cellsize_x}\n    Cell size Y: {cellsize_y}\n')
        # print(f'****************\nDeposition (DoD positive)\n****************\n    Deposition volume: {vol_pos} m^3\n    Area: {area_pos} m^2\n    Number of cells: {masked_raster[masked_raster > 0].count()}\n    Cell size X: {cellsize_x}\n    Cell size Y: {cellsize_y}\n')
        # print(f'****************\nErosion (DoD negative)\n****************\n    Erosion volume: {vol_neg} m^3\n    Area: {area_neg} m^2\n    Number of cells: {masked_raster[masked_raster < 0].count()}\n    Cell size X: {cellsize_x}\n    Cell size Y: {cellsize_y}\n')
    
        results_dict = {'dod': str(Path(dod).name), 
                        'cellsize_x': cellsize_x,
                        'cellsize_y': cellsize_y,
                        'net_vol': vol_net,
                        'net_area': area_net,
                        'deposition_vol': vol_pos,
                        'deposition_area': area_pos,
                        'erosion_vol': vol_neg,
                        'erosion_area': area_neg
                       }
        
        df = pd.DataFrame([results_dict])

    return df

def integer_align_raster(inras, outras):
    '''
    Align raster to integer bounds.
    from gis.stackexchange.com/questions/296770/aligning-many-rasters-using-pyqgis-or-python
    
    Use gdalwarp targetAlignedPixels (-tap) to ensure cell alignment.
    -tap requires -tr argument.  Make sure -tr is set to whole cm by rounding to 2 deminal places.

    inras: Path object to input raster
    outras: Path object to output raster
    '''
    src = gdal.Open(str(inras), gdalconst.GA_ReadOnly)
    srcProj = src.GetProjection()
    ulx, xres_in, xskew, uly, yskew, yres_in  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres_in)
    lry = uly + (src.RasterYSize * yres_in)
    #get target bounds (expanded to nearest integer
    t_ulx = math.floor(ulx)
    t_uly = math.ceil(uly)
    t_lrx = math.ceil(lrx)
    t_lry = math.floor(lry)
    #make sure cell size is rounded to cm
    xres = np.round(xres_in,2)
    yres = np.round(yres_in,2)
    #gdalwarp
    ds = gdal.Warp(str(outras), src, format='GTiff', outputBounds=[t_ulx, t_lry, t_lrx, t_uly], xRes=xres, yRes=yres, targetAlignedPixels=True, resampleAlg=gdal.GRA_Bilinear, creationOptions=['COMPRESS=LZW'])

    print(f'****************\nAlignment complete.\n    Input cell size: {xres}\n    Input origin: {ulx}, {uly}\n\n    Output cell size: {new_cell_size}\n    Output origin: {t_ulx}, {t_uly}\n\nInteger aligned raster written to: {outras}\n****************\n')

def integer_align_resample_raster_bilinear(inras, new_cell_size, outras):
    '''
    Align raster to integer bounds, and resamples to new cell size.
    from gis.stackexchange.com/questions/296770/aligning-many-rasters-using-pyqgis-or-python
    
    Use gdalwarp targetAlignedPixels (-tap) to ensure cell alignment.

    inras: Path object to input raster
    new_cell_size: float value for new cell size
    outras: Path object to output raster
    '''
    src = gdal.Open(str(inras), gdalconst.GA_ReadOnly)
    dest = gdal.Open(str(inras), gdalconst.GA_ReadOnly)
    srcProj = src.GetProjection()
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    #get target bounds (expanded to nearest integer
    t_ulx = math.floor(ulx)
    t_uly = math.ceil(uly)
    t_lrx = math.ceil(lrx)
    t_lry = math.floor(lry)
    #gdalwarp
    ds = gdal.Warp(str(outras), src, format='GTiff', outputBounds=[t_ulx, t_lry, t_lrx, t_uly], xRes=new_cell_size, yRes=new_cell_size, targetAlignedPixels=True, resampleAlg=gdal.GRA_Bilinear, creationOptions=['COMPRESS=LZW'])

    print(f'****************\nAlignment complete.\n    Input cell size: {xres}\n    Input origin: {ulx}, {uly}\n\n    Output cell size: {new_cell_size}\n    Output origin: {t_ulx}, {t_uly}\n\nInteger aligned raster written to: {outras}\n****************\n')


def resample_raster_bilinear(input_raster, new_cell_size, output_raster):
    '''
    Use this function to resample a raster using bilinear interpolation. Value provided in argument 
    will be new cell size for output raster.  Rasters must be geotiff.

    input_raster: Path object to input raster
    new_cell_size: floating point value for cell size of new raster
    output_raster: Path object to resampled output raster
    '''
    # Read the raster using rasterio
    with rasterio.open(input_raster) as src:
        # Get the raster values as a masked numpy array
        # raster_values = src.read(1,masked=True)
        profile = src.profile.copy()
        #get cell size x, y
        cellsize_x, cellsize_y  = src.res

        #determine scale factor
        scale_factor = cellsize_x/new_cell_size

        # resample data to target shape
        data = src.read(
                            out_shape=(
                                src.count,
                                int(src.height * scale_factor),
                                int(src.width * scale_factor)
                                ),
                                resampling=Resampling.bilinear
                            )
    
        # scale image transform
        transform = src.transform * src.transform.scale(
                                                        (src.width / data.shape[-1]),
                                                        (src.height / data.shape[-2])
                                                    )

        profile.update({"height": data.shape[-2],
                    "width": data.shape[-1],
                   "transform": transform,
                   "compress":"lzw"}
                   )

    with rasterio.open(output_raster, "w", **profile) as dataset:
        dataset.write(data)
        
    

def vertical_adjust_raster_uniform(input_raster, adjustment_value, output_raster):
    '''
    Use this function to vertically adjust a raster with a uniform adjustment value. Value provided in argument will be added to input raster
    and saved to output raster.  Rasters must be geotiff.

    input_raster: Path object to input raster
    adjustment_value: floating point value that will be added to input_raster to create output_raster
    output_raster: Path object to adjusted output raster
    '''

    # Read the raster using rasterio
    with rasterio.open(input_raster) as src:
        # Get the raster values as a masked numpy array
        raster_values = src.read(1,masked=True)
        src_meta = src.meta

    # Add value
    out_raster_values = raster_values + adjustment_value

    #Write new raster with same meta
    with rasterio.open(output_raster, 'w', **src_meta) as dst:
        # Get the raster values as a numpy array
        dst.write(out_raster_values, indexes=1)



def vertical_adjust_raster_with_raster(input_raster, adjustment_raster, output_raster):
    '''
    Use this function to apply a vertical adjustment raster (such as an error trend surface) to an input raster. Rasters must be geotiff.

    input_raster: Path object to input raster
    adjustment_raster: adjustment raster that will be added to input_raster to create output_raster
    output_raster: Path object to adjusted output raster
    '''

    with rasterio.open(input_raster) as dataset1, rasterio.open(adjustment_raster) as dataset2:   
        # Retrieve the affine transformation matrices
        transform1 = dataset1.transform
        transform2 = dataset2.transform
    
        # Determine the common extent
        bounds1 = dataset1.bounds
        bounds2 = dataset2.bounds
        intersection_bounds = BoundingBox(
            max(bounds1.left, bounds2.left),
            max(bounds1.bottom, bounds2.bottom),
            min(bounds1.right, bounds2.right),
            min(bounds1.top, bounds2.top)
            )
    
       # Update the transformation matrices to the common extent
        window1 = dataset1.window(*intersection_bounds)
        window2 = dataset2.window(*intersection_bounds)
        transform1 = rasterio.windows.transform(window1, transform1)
        transform2 = rasterio.windows.transform(window2, transform2)
    
        # Read the subset of raster data within the common extent
        input_raster_values = dataset1.read(1, window=window1)
        adjustment_raster_values = dataset2.read(1, window=window2)
    
        # Add value
        out_raster_values = input_raster_values + adjustment_raster_values
        
        # Set nodata values in the difference array
        nodata_mask = np.logical_or(input_raster_values == dataset1.nodata, adjustment_raster_values == dataset2.nodata)
        out_raster_values[nodata_mask] = dataset1.nodata
    
        # Retrieve the metadata from one of the input rasters
        metadata = dataset1.meta
    
        # Update the metadata to reflect the difference raster
        metadata.update(
            count=1,
            height=out_raster_values.shape[0],
            width=out_raster_values.shape[1],
            transform=transform1,
            nodata=dataset1.nodata,
            compress='lzw'
        )
    
        # Create a new raster file for the difference
        output_file = output_raster
        with rasterio.open(str(output_file), 'w', **metadata) as output_dataset:
            output_dataset.write(out_raster_values, 1)
    
    print(f"****************\nRaster adjustment complete.\n    Performed the following operation:\n        {input_raster.name} + {adjustment_raster.name} = {output_raster.name}\n\nOutput written to {output_raster}\n****************\n")


def stable_area_stats(raster, polyshp, add_to_attribute_table: bool=False):
    '''
    Calculate zonal statstics for a polygon shapefile,
    and optionally add them as attributes to attribute table.
    Use this to evaluate systematic and random error for stable areas 
    in DEM of difference.

    raster: Path object to raster
    polyshp: Path object to polygon shapefile
    add_to_attribute_table: boolean to add statistics to attribute table (default False)

    df: dataframe with zonal stats
    '''

    #define MAE function to be added to rasterstats
    def mae(x):
        return np.ma.mean(np.abs(x))

    # Read the shapefile using geopandas
    shapefile = gpd.read_file(polyshp)
    
    # Read the raster using rasterio
    with rasterio.open(raster) as src:
        # Get the raster values as a numpy array
        raster_values = src.read(1)
        # Get the affine transformation from the raster
        transform = src.transform
        # Perform the zonal statistics using rasterstats
        stats = zonal_stats(shapefile, raster_values, affine=transform, nodata=src.nodata, stats=['mean', 'median', 'max', 'min', 'std', 'sum'], add_stats={'mae':mae})
    
    # Extract the statistics and add them as new attributes in the shapefile
    shapefile['raster'] = raster.name
    shapefile['mean'] = [stat['mean'] for stat in stats]
    shapefile['median'] = [stat['median'] for stat in stats]
    shapefile['max'] = [stat['max'] for stat in stats]
    shapefile['min'] = [stat['min'] for stat in stats]
    shapefile['std'] = [stat['std'] for stat in stats]
    shapefile['sum'] = [stat['sum'] for stat in stats]
    shapefile['mae'] = [stat['mae'] for stat in stats]
    shapefile['area_m_sq'] = shapefile['geometry'].area

    df = shapefile[['id','mean','median','max','min','std','sum','mae','area_m_sq']]
       
    if add_to_attribute_table:
        # Save attributes to the updated shapefile
        shapefile.to_file(polyshp)

    #print stats
    print('Zonal statistics calculation complete.\n')

    return df

def SpatiallyCorrelatedRandomErrorAnalysis_DataPrep(raster, polyshp, poly_id_field='id', poly_id=1):
    '''
    Use this function to select a portion of the raster on which to perform the analysis
    to determine spatially correlated random error. 
    This function will clip a raster using a shapefile and return a dataframe with x,y,z values 
    for subsequent semi-variogram analysis.

    raster: Path object to input raster
    polyshp: Path object to polygon shapefile with polygon for stable area
    poly_id_field: name of field in attribute table used to select polygon (str)
    poly_id: polygon ID to select single polygon for analysis (to select from poly_id_field).

    '''

    # Read the polygon shapefile using geopandas
    polygon_gdf = gpd.read_file(polyshp)
    # select only one polygon
    polygon_gdf = polygon_gdf.loc[polygon_gdf[poly_id_field] == poly_id]
    
    # Read the raster using rasterio
    with rasterio.open(raster) as src:
        # Clip the raster using the shapefile geometry
        clipped_raster, clipped_transform = rasterio.mask.mask(src, polygon_gdf.geometry, crop=True)
        clipped_meta = src.meta
    
    # Update the metadata of the clipped raster
    clipped_meta.update({
        'height': clipped_raster.shape[1],
        'width': clipped_raster.shape[2],
        'transform': clipped_transform
    })
    
    # Extract the x, y coordinates from the transform
    x_coords = []
    y_coords = []
    z_values = []
    
    for row in range(clipped_raster.shape[1]):
        for col in range(clipped_raster.shape[2]):
            x, y = rasterio.transform.xy(clipped_transform, row, col)
            value = clipped_raster[0, row, col]  # Assume single-band raster
    
            if value != clipped_meta['nodata']:
                x_coords.append(x)
                y_coords.append(y)
                z_values.append(value)
    
    # Create a pandas DataFrame with x, y, z values
    df = pd.DataFrame({'x': x_coords, 'y': y_coords, 'z': z_values})

    return df

def SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram(df, cellsize, max_range=80, plot: bool=True):
    '''
    Use this function to create and plot a preliminary semivariogram to be used to help estimate 
    the correct range subsequent semivariogram plotting.  Run SpatiallyCorrelatedRandomErrorAnalysis_DataPrep
    function first, to create dataframe of DoD raster cell x,y,z from clipped stable area in DoD.

    df: dataframe with x,y,z columns
    cellsize: cell size of DEM of difference
    max_range: maximum range for prelim semivariogram (something larger than predicted range, default 80 meters)
    plot: boolean to plot semivariogram (default=True)

    '''
    #Get the semi-variogram values through kriging 
    #(using pyinterpolate library, with help from docs here: 
    # https://pyinterpolate.readthedocs.io/en/latest/usage/tutorials/Semivariogram%20Estimation%20(Basic).html
    
    #sample list
    data = df[['x','y','z']].to_numpy()
    
    # Create experimental semivariogram
    experimental_variogram = build_experimental_variogram(input_array=data, step_size=cellsize, max_range=max_range)
    
    if plot:
    #show plot of semivariogram
        experimental_variogram.plot(plot_semivariance=True, plot_covariance=True, plot_variance=True)
        
    return experimental_variogram

def SpatiallyCorrelatedRandomErrorAnalysis_FitSphericalModel(experimental_variogram, nugget=0, plot: bool=True):
    '''
    Use this function to plot a spherical semivariogram.  Must run SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram function first.

    experimental_variogram: variogram object from SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram function
    var_range: range estimated from semivariogram plotted with SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram function
    nugget: nugget, should leave at 0 for this analysis (default: 0 meters)
    plot: boolean to plot semivariogram (default=True)
    
    '''
    #Use plot from "SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram" function to estimate range
    # sill = experimental_variogram.variance
    
    semivariogram_model = TheoreticalVariogram()
    fitted = semivariogram_model.autofit(
        experimental_variogram=experimental_variogram,
        model_types='spherical',
        nugget=0)
    print(f"\n\nModel type: {fitted['model_type']}\nNugget: {fitted['nugget']}\nOptimized Sill: {fitted['sill']} (USE THIS VALUE FOR SPATIALLY CORRELATED RANDOM ERROR VOLUMETRIC UNCERTAINTY CALCULATIONS)\nOptimized Range: {fitted['range']} (USE THIS VALUE FOR SPATIALLY CORRELATED RANDOM ERROR VOLUMETRIC UNCERTAINTY CALCULATIONS)\nRMSE: {fitted['rmse']}")
    
    if plot:
    #show plot of semivariogram
        semivariogram_model.plot()
        
    return fitted

def SpatiallyCorrelatedRandomErrorAnalysis_OptimizedModel(experimental_variogram, nugget=0, plot: bool=True):
    '''
    Use this function to plot determine an optimized semivariogram model for the stable areas.  Must run SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram function first.

    experimental_variogram: variogram object from SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram function
    var_range: range estimated from semivariogram plotted with SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram function
    nugget: nugget, should leave at 0 for this analysis (default: 0 meters)
    plot: boolean to plot semivariogram (default=True)
    
    '''
    #Use plot from "SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram" function to estimate range
    # sill = experimental_variogram.variance
    
    semivariogram_model = TheoreticalVariogram()
    # fitted = semivariogram_model.autofit(
    #     experimental_variogram=experimental_variogram,
    #     model_types='all',
    #     nugget=0,
    #     rang=var_range,
    #     sill=sill)
    # print(f"Chosen model type: {fitted['model_type']}\nNugget: {fitted['nugget']}\nSill: {fitted['sill']}\nRange: {fitted['range']}\nRMSE: {fitted['rmse']}")
    
    fitted = semivariogram_model.autofit(experimental_variogram=experimental_variogram, nugget=0)
    print(f"\n\nOptimized model type: {fitted['model_type']}\nNugget: {fitted['nugget']}\nOptimized Sill: {fitted['sill']} (USE THIS VALUE FOR SPATIALLY CORRELATED RANDOM ERROR VOLUMETRIC UNCERTAINTY CALCULATIONS)\nOptimized Range: {fitted['range']} (USE THIS VALUE FOR SPATIALLY CORRELATED RANDOM ERROR VOLUMETRIC UNCERTAINTY CALCULATIONS)\nRMSE: {fitted['rmse']}")
    
    if plot:
    #show plot of semivariogram
        semivariogram_model.plot()
        
    return fitted

def SpatiallyCorrelatedRandomErrorAnalysis_FitSphericalModel_gstat(df, n_lags=50, use_nugget: bool=False, plot: bool=True):
    '''
    Use this function to plot a spherical semivariogram, using the scikit-gstat module.  Run SpatiallyCorrelatedRandomErrorAnalysis_DataPrep
    function first, to create dataframe of DoD raster cell x,y,z from clipped stable area in DoD.
    from:   https://scikit-gstat.readthedocs.io/en/latest/auto_examples/tutorial_01_getting_started.html
            https://scikit-gstat.readthedocs.io/en/latest/userguide/variogram.html

    df: dataframe with x,y,z columns
    n_lags: number of lags (don't set this too high (>80?) or too low (too small will reduce resolution of output range).  Default 50?
    use_nugget = default False
    plot: boolean to plot semivariogram (default=True)
        
    '''

    coords = df[["x", "y"]].to_numpy()
    val = df[["z"]].to_numpy()

    V = skg.Variogram(coords, val.flatten(), maxlag='median', n_lags=n_lags, use_nugget=use_nugget, normalize=False)
    # fig1 = V.plot(show=False)
    
    V.estimator = 'matheron'
    V.model = 'spherical'

    #get data from variogram and use curve_fit to find range, sill
    xdata = V.bins
    ydata = V.experimental
    p0 = [np.mean(xdata), np.mean(ydata), 0]
    cof, cov =curve_fit(models.spherical, xdata, ydata, p0=p0)
    range, sill, nugget = (cof[0], cof[1], cof[2])
    rmse = V.rmse
    print(f"    range: {range}\n    sill: {sill}\n    rmse:    {rmse}\n    nugget: {nugget}")
    
    if plot:
        #print fit model
        xi =np.linspace(xdata[0], xdata[-1], 100)
        yi = [models.spherical(h, *cof) for h in xi]
        fig2 = plt.plot(xdata, ydata, 'og')
        plt.plot(xi, yi, '-b');
        plt.xlabel('Distance (m)')
        plt.ylabel('Variance')
        plt.show()

    return range, sill, nugget, rmse, V, xdata, ydata

def aggregate_stable_stats(dod, polyshp, poly_id_field='id', poly_id=['all'], max_range=50, dod_for_spatially_correlated_random_error='same', downsamplefrac=1.0, gstat: bool=True, pyinterp: bool=True):
    '''
    Aggregate spatial stats (including spatially correlated random error) for stable area polygons into a table.

    dod: dod path
    polyshp: shapefile with stable area polygons
    poly_id_field: name of field with polygon ids
    poly_id: list of polygon ids for which spatially correlated random error stats will be calculated.  If 'all', all will be used.
    max_range: max. range for spatial variogram estimates
    dod_for_spatially_correlated_random_error:  dod to be used for spatially correlated random error analysis. 
                                                If blank ('same'), the dod referenced in 'dod' argument will be used.
                                                Different dod (downsampled) can be supplied if memory requirements limit resolution.
    downsamplefrac: fraction to downsample input point data before spatially correlated random error calcs, if needed. 
                                                If 1.0, then all points will be used.  If less than 1.0 (0-1), then input points 
                                                will be randomly sampled to fraction of input.
    gstat: Run spatially correlated random error analysis using Scikit gstat. bool default True
    pyinterp: Run spatially correlated random error analysis using pyinterp. bool default True
    '''

    #make sure objects are paths
    dod = Path(dod)
    polyshp = Path(polyshp)

    #set dod params
    if dod_for_spatially_correlated_random_error=='same':
        dod_for_spatially_correlated_random_error = dod

    # Get cell sizes using rasterio
    with rasterio.open(dod) as src:
        dodcellsize_x, dodcellsize_y  = src.res
    if dod_for_spatially_correlated_random_error != dod:
        with rasterio.open(dod_for_spatially_correlated_random_error) as src:
            scedodcellsize_x, scedodcellsize_y  = src.res
    else:
        scedodcellsize_x = dodcellsize_x

    #run stable area stats
    df = stable_area_stats(dod, polyshp, add_to_attribute_table=False)

    #set new col
    df['dod'] = dod.name
    df['dod_cellsize'] = dodcellsize_x
    df['sce_dod'] = dod_for_spatially_correlated_random_error.name
    df['sce_dod_cellsize'] = scedodcellsize_x
    df['sce_max_range'] = max_range
    df['sce_num_pts'] = np.nan
    df['gstat_sill'] = np.nan
    df['gstat_range'] = np.nan
    df['gstat_rmse'] = np.nan
    df['pyinterp_sill'] = np.nan
    df['pyinterp_range'] = np.nan
    df['pyinterp_rmse'] = np.nan

    #iterate through rows (polygons) in df to run spatially correlated error analysis
    for index, row in df.iterrows():
        if 'all' in poly_id or row[poly_id_field] in poly_id:
            #run spatially correlated random error analysis on this polygon
            scedf = SpatiallyCorrelatedRandomErrorAnalysis_DataPrep(dod_for_spatially_correlated_random_error,
                                                                     polyshp,
                                                                     poly_id_field=poly_id_field, 
                                                                     poly_id=row['id'])
                                                                     
            #sample input points if needed
            if downsamplefrac != 1:
                scedf = scedf.sample(frac=downsamplefrac, random_state=1)
            df.loc[index,'sce_num_pts'] = len(scedf)
            
            if gstat:
                try:
                    range, sill, nugget, rmse, _, _, _ = SpatiallyCorrelatedRandomErrorAnalysis_FitSphericalModel_gstat(scedf, n_lags=max_range, use_nugget=False)
                    df.loc[index,'gstat_sill'] = sill
                    df.loc[index,'gstat_range'] = range
                    df.loc[index,'gstat_rmse'] = rmse
                except:
                    print('Problem encountered while executing function SpatiallyCorrelatedRandomErrorAnalysis_FitSphericalModel_gstat.\n    Setting values to -9999 and continuing to next polygon.')
                    df.loc[index,'gstat_sill'] = -9999
                    df.loc[index,'gstat_range'] = -9999
                    df.loc[index,'gstat_rmse'] = -9999

            if pyinterp:
                try:
                    experimental_variogram = SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram(scedf, scedodcellsize_x, max_range=max_range, plot=False);
                    fitted = SpatiallyCorrelatedRandomErrorAnalysis_FitSphericalModel(experimental_variogram)
                    df.loc[index,'pyinterp_sill'] = fitted['sill']
                    df.loc[index,'pyinterp_range'] = fitted['range']
                    df.loc[index,'pyinterp_rmse'] = fitted['rmse']
                except: 
                    print('Problem encountered while executing function SpatiallyCorrelatedRandomErrorAnalysis_CreateSemivariogram or SpatiallyCorrelatedRandomErrorAnalysis_FitSphericalModel.\n    Setting values to -9999 and continuing to next polygon.')
                    df.loc[index,'pyinterp_sill'] = -9999
                    df.loc[index,'pyinterp_range'] = -9999
                    df.loc[index,'pyinterp_rmse'] = -9999

    return df
    
def calculateVolumetricUncertainty(dod, sigma_sys, scre_sill, scre_range, sigma_re, confidence_level=95):
    '''
    Calculate volumetric uncertainty for DoD using spatial stats split out into net change, deposition, and erosion.

    dod: dod path
    sigma_sys: residual systematic error in adjusted (final) DoD. (systematic error)
    scre_sill: sill value in the semivariogram model of vertical residuals in stable areas of the 
                adjusted (final) DoD (spatially correlated random error)
    scre_range: range of the semivariogram model of vertical residuals in stable areas of the 
                adjusted (final) DoD (spatially correlated random error)
    sigma_re: standard deviation of vertical residuals in stable areas of the adjusted (final) DoD (uncorrelated random error)
    confidence_level: confidence level for output uncertainty estimates, default=95 (2 sigma)
    
    '''
    
    #get cell size, volumes, and areas using calculate_volume function
    voldf = calculate_volume(dod)
    #drop, rename cellsize cols
    voldf.drop(['cellsize_y'], axis=1, inplace=True)
    voldf.rename(columns={'cellsize_x': 'cellsize'}, inplace=True)
    #insert total uncertainty columns at left for ease of use
    voldf.insert(4,'net_total_vol_error',None)
    voldf.insert(7,'deposition_total_vol_error',None)
    voldf.insert(10,'erosion_total_vol_error',None)    
    #add columns for additional values
    voldf.loc[:, ['confidence_level','sigma_sys','scre_sill','scre_range','sigma_re',
                 'net_mean_re', 'net_vol_re', 'net_mean_scre', 'net_vol_scre', 'net_mean_sys', 'net_vol_sys',
                 'deposition_mean_re', 'deposition_vol_re', 'deposition_mean_scre', 'deposition_vol_scre','deposition_mean_sys', 'deposition_vol_sys',
                 'erosion_mean_re', 'erosion_vol_re', 'erosion_mean_scre', 'erosion_vol_scre', 'erosion_mean_sys', 'erosion_vol_sys'                    
             ]] = None
             
    #populate
    voldf['confidence_level']=confidence_level
    voldf['sigma_sys']=sigma_sys
    voldf['scre_sill']=scre_sill
    voldf['scre_range']=scre_range
    voldf['sigma_re']=sigma_re

    # Uncorrelated random error (equation 12)
    if confidence_level == 95:
        sigma_re = sigma_re * 1.96
    for prefix in ['net', 'deposition', 'erosion']:
        n = voldf[prefix + '_area'].values[0] / voldf['cellsize'].values[0] / voldf['cellsize'].values[0]
        voldf[prefix + '_mean_re'] = sigma_re / np.sqrt(n)
        voldf[prefix + '_vol_re'] = np.sqrt(n) * np.square(voldf['cellsize'].values[0]) * sigma_re

    # Spatially correlated random error (equations 14 and 16)
    sigma_sc = np.sqrt(scre_sill)

    if confidence_level == 95:
        sigma_sc = sigma_sc * 1.96
    for prefix in ['net', 'deposition', 'erosion']:
        n = voldf[prefix + '_area'].values[0] / voldf['cellsize'].values[0] / voldf['cellsize'].values[0]
        voldf[prefix + '_mean_scre'] = (sigma_sc / np.sqrt(n)) * np.sqrt((np.pi * np.square(scre_range)) / 5 * np.square(voldf['cellsize'].values[0]))
        voldf[prefix + '_vol_scre'] = 0.79 * scre_range * np.sqrt(n) * voldf['cellsize'].values[0] * sigma_sc

    # Systematic error (equation 20)
    # sigma_sys = Residual systematic error in adjusted (final) DoD.
    #    Here the mean of residuals in stable areas of the adjusted (final) DoD is interpreted as the mean systematic uncertainty.
    #    Unsure if this represents 1 sigma or 2 sigma?  Does it need to be multiplied by 1.96 if desired conf_level is 95%?  (I don't know)
    for prefix in ['net', 'deposition', 'erosion']:
        n = voldf[prefix + '_area'].values[0] / voldf['cellsize'].values[0] / voldf['cellsize'].values[0]    
        voldf[prefix + '_mean_sys'] = sigma_sys
        voldf[prefix + '_vol_sys'] = n * np.square(voldf['cellsize'].values[0]) * sigma_sys

    #Calc total vol uncertainties
    for prefix in ['net', 'deposition', 'erosion']:
        voldf[prefix + '_total_vol_error'] = np.sqrt(np.square(voldf[prefix + '_vol_re'].values[0]) 
                                                     + np.square(voldf[prefix + '_vol_scre'].values[0]) 
                                                     + np.square(voldf[prefix + '_vol_sys'].values[0])
                                                    )
    cols = ['dod', 'cellsize', 'net_vol', 'net_total_vol_error', 'net_area', 'net_vol_re', 'net_vol_scre','net_vol_sys']
    netdf = voldf[cols]
    
    cols = ['dod', 'cellsize', 'deposition_vol', 'deposition_total_vol_error', 'deposition_area', 'deposition_vol_re', 'deposition_vol_scre','deposition_vol_sys']
    depositiondf = voldf[cols]
    
    cols = ['dod', 'cellsize', 'erosion_vol', 'erosion_total_vol_error', 'erosion_area', 'erosion_vol_re', 'erosion_vol_scre','erosion_vol_sys']
    erosiondf = voldf[cols]
    
    return voldf, netdf, depositiondf, erosiondf