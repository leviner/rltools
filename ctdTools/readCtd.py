import pandas as pd
from glob import glob
import ctd
import seawater as sw
import re
import numpy as np
import os
from geopy.distance import distance
from datetime import timedelta
import xarray as xr

def castMeta(file): # get the lat/lon data
    with open(file) as fd:
        try:
            for line in fd:
                match = re.search(r'NMEA', line)
                if match:
                    if re.search('(?<=Latitude)(.*)', line):
                        lat = re.search('(?<=Latitude)(.*)', line)
                        lat = lat.group()
                        lat = float(lat.split()[1])+(float(lat.split()[2])/60)
                    if re.search('(?<=Longitude)(.*)', line):
                        lon = re.search('(?<=Longitude)(.*)', line)
                        lon = lon.group()
                        lon = -1*float(lon.split()[1])+(float(lon.split()[2])/60)
                    if re.search('(?<=UTC)(.*)', line):
                        dt = re.search('(?<=UTC)(.*)', line)
                        dt = dt.group()
                        dt = pd.to_datetime(dt.split('=')[1])
        except:
            lat, lon, dt = [np.nan for i in range(3)]
        return lat, lon, dt

def readCnv(file): # get the T/S Data by Cast
    df = ctd.from_cnv(file)
    df= df[['t090C','sal00']].rename(columns={"t090C": "temp", "sal00": "sal"})
    df['cast'] = float(os.path.basename(file)[3:-4])
    return df

def formatDf(file): # convert pressure to depth
    lat, lon, dt = castMeta(file)
    df = readCnv(file)
    df['depth'] = sw.dpth(df.index,lat)
    df = df.set_index(['cast'])
    df['lat'] = lat
    df['lon'] = lon
    return df, lat, lon, dt, df.index.unique()[0]

def cnv2table(dir): #produce ctd and key dataframes from SB .cnv files
    if os.path.isdir(dir) is True:
        files=glob(dir+'*.cnv')
    else:
        files = [dir]
    dfCtd = pd.DataFrame()
    lats, lons, casts,time = [[] for i in range(4)]
    for file in files:
        df, lat, lon, dt, cast = formatDf(file)
        dfCtd = dfCtd.append(df)
        lats.append(lat)
        lons.append(lon)
        casts.append(cast)
        time.append(dt)
    dfCtd.loc[dfCtd.temp == -9.990e-29] = np.nan # replace seabird nan with numpy nan
    dfCtdKey = pd.DataFrame({'cast':casts,'lat':lats,'lon':lons,'time':time})
    return dfCtd, dfCtdKey

def csv2table(dir): # produce ctd and key dataframes from QC'd ctd data tables
    if os.path.isdir(dir) is True:
        files=glob(dir+'*.csv')
    else:
        files = [dir]
    dfCtd,dfCtdKey = [pd.DataFrame() for i in range(2)]
    for file in files:
        dfCur = pd.read_csv(file,skiprows=[1])
        cast = []
        for name in dfCur.profile_id.values:
            cast.append(float(name.split('_')[0][-3:]))
        dfCur['cast'] = cast
        dfCur['depth'] = sw.dpth(dfCur.pressure,dfCur.latitude)
        dfCur = dfCur.rename(columns={"T_28": "temp", "S_41": "sal","latitude":"lat","longitude":"lon"})
        dfCtdCur = dfCur[["cast","temp","sal","depth","lat","lon"]]
        dfCtdCur = dfCtdCur.set_index('cast')
        dfCtd = dfCtd.append(dfCtdCur)

        dfKeyCur = dfCur[["cast","lat","lon","time"]]
        dfKeyCur.lon = dfKeyCur.lon-360
        dfKeyCur.time = pd.to_datetime(dfKeyCur.time)
        dfCtdKey = dfCtdKey.append(dfKeyCur.drop_duplicates())
    return dfCtd, dfCtdKey

def nearestCtd(dfTrawls, dfCtdKey): # find nearest ctd to trawl event
    nearDist, nearCast,tEvent = [[] for i in range(3)]
    for i in range(len(dfTrawls)):
        dist = 100
        cast = 0
        for ii in range(len(dfCtdKey)):
            curDist = distance((dfTrawls.EQ_LATITUDE.values[i],dfTrawls.EQ_LONGITUDE.values[i]),(dfCtdKey.lat.values[ii],dfCtdKey.lon.values[ii])).nm
            curTime = pd.to_datetime(dfTrawls.EQ_TIME.values[i], format="%d-%b-%y %H.%M.%S.%f %p")
            if (curDist < dist) & ((curTime - dfCtdKey.time.values[ii]) < timedelta(1)):
                dist = curDist
                cast = dfCtdKey.cast.values[ii]
        if cast == 0:
            print('Trawl object '+str(i)+' had no match')
        nearDist.append(dist)
        nearCast.append(cast)
        tEvent.append(dfTrawls.EVENT_ID.values[i])
    dfMatched = pd.DataFrame({'Event_id':tEvent,'ctdCast':nearCast,'ctdDist':nearDist})
    return dfMatched

def eventTemps(dfCtd, dfCtdKey, dfTrawls):
    dfMatched = nearestCtd(dfTrawls, dfCtdKey)
    tempOpen, tempBot, tempSurf, tempCol, tempRange = [[] for i in range(5)]
    for index, trawl in dfMatched.iterrows():
        curCast = dfCtd[dfCtd.index == trawl.ctdCast]
        tempOpen.append(np.nanmean(curCast.temp[(curCast.depth > (dfTrawls.AVG_HEAD_ROPE_DEPTH[dfTrawls.EVENT_ID == trawl.Event_id].values[0])) & \
            (curCast.depth < (dfTrawls.AVG_HEAD_ROPE_DEPTH[dfTrawls.EVENT_ID == trawl.Event_id].values[0]+ \
            dfTrawls.AVG_NET_HORI_OPENING[dfTrawls.EVENT_ID == trawl.Event_id].values[0]))]))
        tempRange.append(np.nanmean(curCast.temp[(curCast.depth > (dfTrawls.MIN_HEAD_ROPE_DEPTH[dfTrawls.EVENT_ID == trawl.Event_id].values[0])) & \
            (curCast.depth < (dfTrawls.MAX_HEAD_ROPE_DEPTH[dfTrawls.EVENT_ID == trawl.Event_id].values[0]+ \
            dfTrawls.AVG_NET_HORI_OPENING[dfTrawls.EVENT_ID == trawl.Event_id].values[0]))]))
        tempBot.append(np.nanmean(curCast.temp[curCast.depth > (curCast.depth.max()-5)]))
        tempSurf.append(np.nanmean(curCast.temp[curCast.depth < 5]))
        tempCol.append(np.nanmean(curCast.temp))
    dfMatched['tempCol'] = tempCol
    dfMatched['tempBot'] = tempBot
    dfMatched['tempSurf'] = tempSurf
    dfMatched['tempOpen'] = tempOpen
    dfMatched['tempRange'] = tempRange
    return dfMatched

def nc2csv(file,year): # save PMEL standard netCDF files to csv files
    ds = xr.open_dataset(file,decode_times=False)
    df = ds.to_dataframe().reset_index()
    df['StationID'] = file.split('_')[2]
    df['Year'] = year
    df['StationNumber'] = int(file.split('_')[2][-3:])
    df = df[['StationID','Year','StationNumber','lat','lon','dep','PAR_905','T_28','S_41','OST_62','ST_70']]
    df = df.rename(columns={'lat':'LatitudeStart DD','lon':'LongitudeStart DD','dep':'Depth (m)','S_41':'PrimarySalinity PSU','T_28':'PrimaryTemperature deg. C',
                   'PAR_905':'PARIrradiance uE','OST_62':'PrimaryOxygenSat %','ST_70':'SigmaTheta kg/m3'})
    df.to_csv(file.split('.')[0]+'.csv', index=False)
    return True