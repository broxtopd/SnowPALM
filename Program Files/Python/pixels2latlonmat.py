from osgeo import osr
from osgeo import gdal
import sys
from osgeo.gdalconst import *
import time
import numpy as np
import copy
# USAGE: python pixels2latlonmat.py data.tif lat.tif lon.tif

# The following method translates given latitude/longitude pairs into pixel locations on a given GEOTIF
# INPUTS: geotifAddr - The file location of the GEOTIF
#      latLonPairs - The decimal lat/lon pairings to be translated in the form [[lat1,lon1],[lat2,lon2]]
# OUTPUT: The pixel translation of the lat/lon pairings in the form [[x1,y1],[x2,y2]]
# NOTE:   This method does not take into account pixel size and assumes a high enough 
#	  image resolution for pixel size to be insignificant
def latLonToPixel(ds, latLonPairs):
	# Get a geo-transform of the dataset
	gt = ds.GetGeoTransform()
	# Create a spatial reference object for the dataset
	srs = osr.SpatialReference()
	srs.ImportFromWkt(ds.GetProjection())
	# Set up the coordinate transformation object
	dsrs = osr.SpatialReference()
	dsrs.ImportFromEPSG(4326)
	srsLatLong = dsrs.CloneGeogCS()
	ct = osr.CoordinateTransformation(srsLatLong,srs)
	# Go through all the point pairs and translate them to latitude/longitude pairings
	pixelPairs = []
	for point in latLonPairs:
		# Change the point locations into the GeoTransform space
		(point[1],point[0],holder) = ct.TransformPoint(point[1],point[0])
		# Translate the x and y coordinates into pixel values
		x = (point[1]-gt[0])/gt[1]
		y = (point[0]-gt[3])/gt[5]
		# Add the point to our return array
		pixelPairs.append([int(x),int(y)])
	return pixelPairs
# The following method translates given pixel locations into latitude/longitude locations on a given GEOTIF
# INPUTS: geotifAddr - The file location of the GEOTIF
#      pixelPairs - The pixel pairings to be translated in the form [[x1,y1],[x2,y2]]
# OUTPUT: The lat/lon translation of the pixel pairings in the form [[lat1,lon1],[lat2,lon2]]
# NOTE:   This method does not take into account pixel size and assumes a high enough 
#	  image resolution for pixel size to be insignificant
def pixelToLatLon(ds,pixelPairs):
	# Get a geo-transform of the dataset
	gt = ds.GetGeoTransform()
	# Create a spatial reference object for the dataset
	srs = osr.SpatialReference()
	srs.ImportFromWkt(ds.GetProjection())
	# Set up the coordinate transformation object
	dsrs = osr.SpatialReference()
	dsrs.ImportFromEPSG(4326)
	srsLatLong = dsrs.CloneGeogCS()
	ct = osr.CoordinateTransformation(srs,srsLatLong)
	# Go through all the point pairs and translate them to pixel pairings
	latLonPairs = []
	for point in pixelPairs:
		# Translate the pixel pairs into untranslated points
		ulon = point[0]*gt[1]+gt[0]
		ulat = point[1]*gt[5]+gt[3]
		# Transform the points to the space
		(lon,lat,holder) = ct.TransformPoint(ulon,ulat)
		# Add the point to our return array
		latLonPairs.append([lat,lon])
 
	return latLonPairs
    
if __name__=='__main__':
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    
    # open the image
    inDs = gdal.Open(argv[1], gdal.GA_ReadOnly)
    if inDs is None:
        print('Could not open ' + argv[1])
        sys.exit(1)
        
    driver = inDs.GetDriver()
    
    # get image size
    rows = inDs.RasterYSize
    cols = inDs.RasterXSize
    
    inBand = inDs.GetRasterBand(1)
    inData = inBand.ReadAsArray(0, 0, cols, rows).astype(np.float32)
    
    starttime = time.time()
    pixelPairs = []
    for y in range(rows):
        for x in range(cols):
            pixelPairs.append([x,y])
   
    latLonPairs = pixelToLatLon(inDs,pixelPairs)

    outData1 = copy.deepcopy(inData)
    outData2 = copy.deepcopy(inData)
    c = 0
    for y in range(rows):
        for x in range(cols):
            outData1[y,x] = latLonPairs[c][0]
            outData2[y,x] = latLonPairs[c][1]
            c = c+1
    
    outDs = driver.Create(argv[2], cols, rows, 1, GDT_Float32)
    if outDs is None:
      print('Could not create ' + argv[2])
      sys.exit(1)
    outBand = outDs.GetRasterBand(1)

    # write the output data
    outBand.WriteArray(outData1, 0, 0)
    outBand.FlushCache()
    stats = outBand.GetStatistics(0, 1)

    # georeference the image and set the projection
    outDs.SetGeoTransform(inDs.GetGeoTransform())
    outDs.SetProjection(inDs.GetProjection())
    outDs = None

    outDs = driver.Create(argv[3], cols, rows, 1, GDT_Float32)
    if outDs is None:
      print('Could not create ' + argv[3])
      sys.exit(1)
    outBand = outDs.GetRasterBand(1)

    # write the output data
    outBand.WriteArray(outData2, 0, 0)
    outBand.FlushCache()
    stats = outBand.GetStatistics(0, 1)

    # georeference the image and set the projection
    outDs.SetGeoTransform(inDs.GetGeoTransform())
    outDs.SetProjection(inDs.GetProjection())
    outDs = None

    inDs = None
    
    print('Done creating matrix of latitudes: ' + str(time.time()-starttime) + ' seconds')   
    