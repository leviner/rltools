"""
This is where I'll compile various maps that utilize the AK etopo1 file.

Usage:
import akMaps
maps = akMaps.akMaps()
maps.chukchiScatter(lat = np.array([68,68,69]),lon = np.array([-170,-171,-170]),values = np.array([20,35,27]),max=40,extent = [-155,-175,64,73])
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import rcParams
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

class akMaps:
    
    def __init__(self):
        etopoFile = os.path.dirname(__file__)+'/etopo1_bedrock.asc'
        self.setOlevels(etopoFile)
    
    def setOlevels(self,etopoFile):
        #Initialize the topography
        topo_file = open(etopoFile, 'r')
        ncols = int(topo_file.readline().split()[1])
        nrows = int(topo_file.readline().split()[1])
        xllcorner = float(topo_file.readline().split()[1])
        yllcorner = float(topo_file.readline().split()[1])
        cellsize = float(topo_file.readline().split()[1])
        topo_file.close()
        dres = 1
        etopo = np.loadtxt(etopoFile, skiprows=5)
        etopo[:nrows+1, :] = etopo[nrows+1::-1, :]
        self.etopo = etopo[::dres, ::dres]

        # Create longitude and latitude vectors for etopo
        lons = np.arange(xllcorner, xllcorner+cellsize*ncols, cellsize)[::dres]
        lats = np.arange(yllcorner, yllcorner+cellsize*nrows, cellsize)[::dres]
        self.rlons, self.rlats = (np.meshgrid(lons[:-1],lats[:-1]))

        # Draw etopo1, first for land and then for the ocean, with different colormaps
        llevels = np.arange(-500,2251,100) # check etopo.ravel().max()
        #lcs = m.contourf(rlons, rlats, etopo, llevels, cmap=cm.terrain)
        #olevels = np.arange(-3500,1,100) # check etopo.ravel().min()
        self.olevels1 = [-200,-150,-100,-50,0]
        self.olevels2 = [-2000,-1500,-1000,-500]
        self.olevels3 = [0,10000]
        self.olevels4 = [0]

    
    def chukchiScatter(self,lat,lon,values,max,extent=[-155,-171,59,73.5]):
        figure = plt.figure(figsize=(30,20))
        rcParams['contour.negative_linestyle'] = 'solid'
        ax1=plt.subplot(141,projection=ccrs.Mercator())
        rcParams['lines.linewidth'] = .5
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linestyle='--',zorder=3)
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlocator = mticker.FixedLocator([-172, -168, -164, -160, -156, -154])
        gl.ylocator = mticker.FixedLocator([57,60,62,64,66,68,70,72,74])
        gl.xlabel_style = {'size':16}
        gl.ylabel_style = {'size':16}
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        c = plt.cm.jet(values/int(max))
        a = ax1.scatter(lon, lat,s=40+5**np.log10(values),facecolors='None',edgecolors = c,transform=ccrs.Geodetic(),cmap='jet',zorder=6, vmin = 0,vmax=4)
        s = a.set_clim([0,4])

        lfill = ax1.contourf(self.rlons, self.rlats, self.etopo, self.olevels3, colors ='grey',transform=ccrs.PlateCarree(),zorder=4)#cmap=cm.ocean)
        cso1 = ax1.contour(self.rlons, self.rlats, self.etopo, self.olevels1, colors ='grey',transform=ccrs.PlateCarree(),zorder=1)#cmap=cm.ocean)
        cso2 = ax1.contour(self.rlons, self.rlats, self.etopo, self.olevels2, colors ='lightgrey',transform=ccrs.PlateCarree(),zorder=2)#cmap=cm.ocean)
        cso4 = ax1.contour(self.rlons, self.rlats, self.etopo, self.olevels4, colors ='k',transform=ccrs.PlateCarree(),zorder=5)#cmap=cm.ocean)
        ax1.set_extent(extent)