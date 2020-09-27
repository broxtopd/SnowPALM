import sys,os
import gdal
from gdalconst import *

# Get georeferencing information from a raster file and print to text file

src = sys.argv[1]
fname_out = sys.argv[2]

ds = gdal.Open(src, GA_ReadOnly)
if ds is None:
    print('Content-Type: text/html\n')
    print('Could not open ' + src)
    sys.exit(1)

# Get the geotransform, the georeferencing, and the dimensions of the raster to match
transform = ds.GetGeoTransform()
wkt = ds.GetProjection()
rows = ds.RasterYSize
cols = ds.RasterXSize
ulx = transform[0]
uly = transform[3]
pixelWidth = transform[1]
pixelHeight = transform[5]
lrx = ulx + (cols * pixelWidth)
lry = uly + (rows * pixelHeight)

f_out = open(fname_out, 'w')
f_out.write(str(pixelWidth) + '\n')
f_out.write(str(pixelHeight) + '\n')
f_out.write(str(ulx) + '\n')
f_out.write(str(uly) + '\n')
f_out.write(str(lrx) + '\n')
f_out.write(str(lry))
f_out.close()
