import tifffile
import matplotlib.pyplot as plt
import datetime
from glob import glob
tifFiles = glob('X:/DY2402/instantEchoOutput/*.tif')

curFile = tifFiles[0]
with tifffile.TiffFile(curFile) as tiff:
    image_data = tiff.asarray()
    extra_tags = tiff.pages[0].tags
    
x = [datetime.datetime.utcfromtimestamp(k) for k in extra_tags[0].value]
y = extra_tags[1].value
f = extra_tags[2].value

for i in range(len(image_data)):
    fig = plt.figure()
    c = plt.pcolormesh(x, y, image_data[i].T, vmin=-70, vmax=-30)
    c.cmap.set_under('w')
    plt.gca().invert_yaxis()
    plt.gcf().autofmt_xdate()
    plt.title(f[i])
    plt.colorbar(c)
    plt.show()

    
