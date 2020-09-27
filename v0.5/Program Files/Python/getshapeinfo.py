import sys,os
import ogr

# Get georeferencing information, as well as value of 'Name' field (if 
# exists) from a sahpefile and print to text file

src = sys.argv[1]
fname_out = sys.argv[2]

driver = ogr.GetDriverByName('ESRI Shapefile')
ds = driver.Open(src, 0)

if ds is None:
    print('Content-Type: text/html\n')
    print('Could not open ' + src)
    sys.exit(1)

ulx = 1E30
uly = -1E30
lrx = -1E30
lry = 1E30
names = []
layer = ds.GetLayer()
ldefn = layer.GetLayerDefn()
schema = [ldefn.GetFieldDefn(n).name  for n in range(ldefn.GetFieldCount())]
for feature in layer:
    geom = feature.GetGeometryRef()
    (minX, maxX, minY, maxY) = geom.GetEnvelope()
    ulx = min(ulx,minX)
    uly = max(uly,maxY)
    lrx = max(lrx,maxX)
    lry = min(lry,minY)
    if "Name" in schema:
        names.append(feature.GetField("Name"))

f_out = open(fname_out, 'w')
f_out.write(str(ulx) + '\n')
f_out.write(str(uly) + '\n')
f_out.write(str(lrx) + '\n')
f_out.write(str(lry))
for name in names:
    f_out.write('\n' + str(name))
f_out.close()
