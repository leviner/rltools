t = xr.open_dataset('C:/Users/rober/Downloads/hawaii_soest_794e_6df2_6381_6464_6c66_07de.nc') # global netcdf file
# split up around -180
dsEast = t.sel(longitude=slice(-180,-150))
dsWest = t.sel(longitude=slice(150,180))

# revise longitude labels
dsWest['longitude2'] = dsWest.longitude-360
dsWest = dsWest.swap_dims({'longitude':'longitude2'})
dsWest = dsWest.drop('longitude')

dsEast['longitude2'] = dsEast.longitude
dsEast = dsEast.swap_dims({'longitude':'longitude2'})
dsEast = dsEast.drop('longitude')

# merge the datasets, removing the redundant data at -180
dsAll = dsWest.merge(dsEast.where(dsEast.longitude2>-180,drop=True))

#revert to original coordinate labels and save
dsAll = dsAll.rename({'longitude2':'longitude'})
dsAll.to_netcdf('C:/Users/rober/Downloads/PAR.nc')