# Read in .cnv file to get data
import ctd
import os
from glob import glob
import pandas as pd
import numpy as np
import re
import gsw
import tsCalc
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import argparse

import sys
#sys.path.insert(0, "R:/applications/MaceFunctions")
sys.path.insert(0, "G:/WPy64-31150/applications/MaceFunctions")

from MaceFunctions.echolab2.instruments import EK80
from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap
from MaceFunctions.echolab2.plotting.matplotlib import echogram
from MaceFunctions.echolab2.processing import  line, grid, integration
    
def readCnv(file): # get the T/S Data by Cast
    df = ctd.from_cnv(file)
    df= df[['t090C','sal00']].rename(columns={"t090C": "temp", "sal00": "sal"})
    df['cast'] = float(re.search(r'\d+', os.path.basename(file)).group())
    return df

def cnv2table(file): # convert pressure to depth
    df = readCnv(file)
    df['depth'] = gsw.z_from_p(df.index,50)
    return df.reset_index()

def castaway2table(ctdFile):
    dfCtd = pd.read_csv(ctdFile)
    dfCtd.rename(columns={'Depth (Meter)':'depth','Temperature (Celsius)':'temp','Salinity (Practical Salinity Scale)':'sal','Sound velocity (Meters per Second)':'sound_speed'},inplace=True)
    return dfCtd[['temp','sal','depth','sound_speed']]

def targetEnv(ctdFile,d,sphere,sphereRange):
    lat = float(d.nmea_data.get_datagrams('GGA')['GGA']['data'][0].lat[:2])+\
          (float(d.nmea_data.get_datagrams('GGA')['GGA']['data'][0].lat[2:])/60)

    tsbw = {18000:{512:1750,1024:1570},38000:{512:3280,1024:2430},70000:{512:4630,1024:2830},120000:{512:5490,1024:2990},200000:{512:590,1024:3050},333000:{512:590,1024:3050}}
    if ctdFile.split('.')[-1]== 'cnv': # For ship CTD data
        dfCTD= cnv2table(ctdFile)
    elif ctdFile.split('.')[-1]== 'csv': # For castaway data
        dfCTD = castaway2table(ctdFile)
    dfCTD = dfCTD.reset_index()
    materialPoperties = tsCalc.material_properties()['Tungsten carbide']
    
    sphereEnv = dfCTD.iloc[(dfCTD['depth'] - np.floor(sphereRange)-10).abs().argsort()[:1]]
    c, rho = tsCalc.water_properties(sphereEnv.sal.values, sphereEnv.temp.values, sphereEnv.depth.values, lon=0.0, lat=lat)
    
    if d.is_cw():
        f = d.frequency[0]
        pl = int(np.round(d.pulse_duration[0]*1000000))
        fr, ts = tsCalc.freq_response(f-tsbw[f][pl]/2, f+tsbw[f][pl]/2, sphere/1000/2, c, materialPoperties['c1'], materialPoperties['c2'], rho, materialPoperties['rho1'], fstep=100)
    refTS = 10*np.log10(np.mean(10**(ts/10)))

    return fr, ts, refTS

def getSV(f, files=None,fm=True):
    if isinstance(files, list):
        files = sorted(files)
    else:
        files = [files]
    print('N files:',len(files))

    print('Reading raw files...')
    curDate = files[0].split('-')[1]
    #change this to echosounder package
    ek80 = EK80.EK80()
    ek80.read_raw(files,frequencies=[f])
    d = ek80.get_channel_data(frequencies=f)
    if fm:
        d = d[f][[i for i, val in enumerate([dc.is_fm() for dc in d[f]]) if val][0]]
    else:
        d = d[f][[i for i, val in enumerate([dc.is_cw() for dc in d[f]]) if val][0]]

    cal = d.get_calibration()

    d_sv = d.get_Sv()

    return d,cal,d_sv

def calEchogram(d_sv,f,sphereRange =None,sphereRangeTol=1):
    calview = {18000:2000,38000:2000,70000:2700,120000:5000,200000:8000,333000:8000}
    if sphereRange :
        sphereRange  = sphereRange 
    else:
        sphereRange  = sphereDict[survey][leg][sphere]
    d_sv1 = d_sv.view((0,-1,1),(0,calview[f],1))
    fig_1 = figure(figsize=(12,9))
    eg = echogram.Echogram(fig_1, d_sv1,threshold=[-90,-30])
    eg.add_colorbar(fig_1)
    show()

    sphereRange = [i for i, x in enumerate(d_sv.range) if np.abs(x-sphereRange )<sphereRangeTol]
    d_sv2 = d_sv.view((0,-1,1),(sphereRange[0],sphereRange[-1],1))
    fig_1 = figure(figsize=(12,3))
    eg = echogram.Echogram(fig_1, d_sv2,threshold=[-90,-30])
    eg.add_colorbar(fig_1)
    show()

def detectSingleTargets(test_data,cal,detectParms,singleTargets):
    Sp = test_data.get_Sp(calibration=cal)
    along,athwart = test_data.get_physical_angles(calibration=cal)

    for ping in range(Sp.n_pings):
        cpv = 40 * np.log10(Sp.range) + 2 * cal.absorption_coefficient[ping] * Sp.range
        calPower = Sp.data[ping] - cpv
        maxima = argrelextrema(calPower, np.greater)
        

        for l in maxima[0]:
            PLDLval = calPower[l]-detectParms.PLDL
            if np.where(calPower[l:]< PLDLval)[0].size >0:
                right = l+ np.where(calPower[l:]< PLDLval)[0][0]-1
            else:
                continue

            if np.where(calPower[:l]< PLDLval)[0].size >0:
                left = np.where(calPower[:l]< PLDLval)[0][-1]+1
            else:
                continue
            
            xLeft = left+(PLDLval-calPower[left])/(calPower[left+1]-calPower[left])
            xRight = right+(PLDLval-calPower[right])/(calPower[right+1]-calPower[right])
            normWidth = (xRight - xLeft) / 4
            if (normWidth > detectParms.maxNormPulseLen) | (normWidth < detectParms.minNormPulseLen):
                continue
            
            eStartIdx = left + 1
            eEndIdx = right - 1
            if (eEndIdx - eStartIdx) < 1:
                continue
            
            peakAlong = along.data[ping][l]
            peakAthwart = athwart.data[ping][l]

            al = 2 * peakAlong / cal.beam_width_alongship[ping]
            at = 2 * peakAthwart / cal.beam_width_athwartship[ping]
            beamComp = 6.0206 * (al**2 + at**2 - (0.18 * al**2 * at**2))

            if (beamComp > detectParms.maxBeamComp):
                continue
            
            alongTarget = along.data[ping][eStartIdx:eEndIdx]
            athwartTarget = athwart.data[ping][eStartIdx:eEndIdx]

            sdAlng = np.std(alongTarget)
            if (sdAlng > detectParms.maxSDalong):
                continue
            sdAthw = np.std(athwartTarget)
            if (sdAthw > detectParms.maxSDathwart):
                continue
            r = (sum(Sp.range[eStartIdx:eEndIdx] * calPower[eStartIdx:eEndIdx])) /  sum(calPower[eStartIdx:eEndIdx]) -  (cal.sound_speed * cal.pulse_duration) / 4
            if (r > detectParms.excludeBelow) | (r < detectParms.excludeAbove):
                continue

            uTS = calPower[l] + (40 * np.log10(r)) +  (2 * cal.absorption_coefficient[ping] * r)
            cTS = uTS + beamComp
            if (cTS < detectParms.threshold):
                continue

            singleTargets.ping = np.append(singleTargets.ping, ping)
            singleTargets.r = np.append(singleTargets.r, r)
            singleTargets.uTS = np.append(singleTargets.uTS, uTS)
            singleTargets.cTS = np.append(singleTargets.cTS, cTS)
            singleTargets.peakAthwart = np.append(singleTargets.peakAthwart, peakAthwart)
            singleTargets.peakAlong = np.append(singleTargets.peakAlong, peakAlong)
            singleTargets.sdAlng = np.append(singleTargets.sdAlng, sdAlng)
            singleTargets.sdAthw = np.append(singleTargets.sdAthw, sdAthw)
            singleTargets.normWidth = np.append(singleTargets.normWidth, normWidth)
    
    return singleTargets

def intCalSingleTarget(d_sv,singleTargetsResults,cal,refTS=None,sphereRange=20,sphereRangeTol=1,plot=True):
    
    sphereHits = np.where(np.abs(singleTargetsResults.r-sphereRange)<(sphereRangeTol/2))
    
    plt.plot(singleTargetsResults.r[sphereHits])
    show()

    d_svOnAxis = d_sv.copy()
    pingsNoSphere = np.arange(d_sv.n_pings)[~np.in1d(np.arange(d_sv.n_pings),np.unique(singleTargetsResults.ping[sphereHits]).astype(int))]
    d_svOnAxis.delete(index_array=pingsNoSphere)
    d_svOnAxis.range = d_sv.range

    observedTS = 10*np.log10(np.mean(10**(singleTargetsResults.cTS[sphereHits]/10))) # Get my observed TS by averaging the on-axis hits
    if plot:
        plt.figure(figsize=(8,4))
        plt.subplot(1,2,1)
        plt.hist(singleTargetsResults.cTS[sphereHits],bins=100);
        plt.subplot(1,2,2)
        plt.plot(singleTargetsResults.peakAlong[sphereHits],singleTargetsResults.peakAthwart[sphereHits],'.');
        show()

    meanRange = np.mean(singleTargetsResults.r[sphereHits])
    
    # set the sphere lines for the upper/lower integration 
    upperLine = line.line(ping_time=d_svOnAxis.ping_time, data=meanRange-(sphereRangeTol*.5)) 
    lowerLine = line.line(ping_time=d_svOnAxis.ping_time, data=meanRange+(sphereRangeTol*1.5)) 

    
    fig_1 = figure(figsize=(12,3))
    eg = echogram.Echogram(fig_1, d_svOnAxis,threshold=[-90,-30])
    eg.add_colorbar(fig_1)
    show()
    
    i = integration.integrator(min_threshold_applied=False) # set up the integrator
    g = grid.grid(interval_length=10000, interval_axis='ping_number',data=d_svOnAxis, layer_axis='range',layer_thickness=100) # thors a grid thats bigger than all the data
    i_2 = i.integrate(d_svOnAxis,g,exclude_above_line=upperLine,exclude_below_line=lowerLine) # integrate with our lines 

    # Calculate the reference nasc from the reference TS, EBA during data colection, and range
    refNasc = (10**(refTS/10)*(1852**2)*4*np.pi)/((10**(np.unique(cal.equivalent_beam_angle)[0]/10))*((meanRange)**2))
    # Now print it all out
    print('***Cal Results***')
    print('Observed TS: ',observedTS) # Observed TS
    print('EBA: ',cal.equivalent_beam_angle[0]) # Observed TS
    print('Target range: ',meanRange) # Mean Range
    print('Observed sA: ',i_2.nasc[0][0]) # Integrated sA
    print('Reference sA: ',refNasc) # reference sA
    print('Number of good pings: ',len(d_svOnAxis.ping_time)) # How many hits were included
    print('Gain during collection: ',np.unique(cal.gain+cal.sa_correction)[0]) # what the gain used during data collection was
    # get our 'Sv gain' as the difference between the used gain (+sa corr, which is 0 for FM) and the ratio of the sA
    print('New Sv gain: ',np.unique((cal.gain+cal.sa_correction)-(10*np.log10(refNasc/i_2.nasc))/2))
    # Break that down into our 'Calc gain', which is the TS gain and the 'Sv corr', which is basically the equivalent of the sa correction
    print('Calc gain: ',((observedTS-refTS)/2)+np.unique(cal.gain), 'Sa corr: ',(np.unique((cal.gain+cal.sa_correction)-(10*np.log10(refNasc/i_2.nasc))/2)) -(((observedTS-refTS)/2)+np.unique(cal.gain)))
    return np.unique((cal.gain+cal.sa_correction)-(10*np.log10(refNasc/i_2.nasc))/2) # This returns the Sv gain

class detectParmsInit():
    # Set something up to mimic TS range from lobes
    PLDL = 6
    maxNormPulseLen = 20
    minNormPulseLen = .1
    maxBeamComp = .1
    maxSDalong = .6
    maxSDathwart = .6
    excludeBelow = 1e10
    excludeAbove = 0
    threshold = -50

class singleTargetsInit():
    ping = np.array([])
    r = np.array([])
    uTS = np.array([])
    cTS = np.array([])
    peakAthwart = np.array([])
    peakAlong = np.array([])
    sdAlng = np.array([])
    sdAthw = np.array([])
    normWidth = np.array([])

def main(args):
    if os.path.isdir(args.file):
        args.file = glob(args.file+'/*.raw')
    d, cal,d_sv =  getSV(args.freq,files=args.file,fm=args.FM)
    fr, ts, refTS = targetEnv(args.ctd_file,d,args.sphere_diameter,args.sphere_range)
    print('Reference TS: ',refTS)
    calEchogram(d_sv,args.freq,sphereRange =args.sphere_range,sphereRangeTol=args.sphere_range_tol)
    detectParms = detectParmsInit()
    singleTargets = singleTargetsInit()
    print('Detecting single targets...')
    singleTargetsResults = detectSingleTargets(d,cal,detectParms,singleTargets)
    svgain = intCalSingleTarget(d_sv,singleTargetsResults,cal,sphereRange=args.sphere_range,plot=True,sphereRangeTol=args.sphere_range_tol,refTS=refTS)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Lets Cal!')
    parser.add_argument('file', type=str, help='File to calibrate')
    parser.add_argument('--ctd_file', type=str, help='CTD file from calibration, either seabird .cnv or castaway .csv')
    parser.add_argument('--FM', type=bool, default=False, help='FM?')
    parser.add_argument('--freq', type=int, default=38000, help='Frequency to calibrate')
    parser.add_argument('--sphere_diameter', type=float, default=38.1, help='diameter of calibration sphere')
    parser.add_argument('--sphere_range', type=float, default=30, help='range to the sphere')
    parser.add_argument('--sphere_range_tol', type=float, default=1, help='tolerance for sphere range')
    args = parser.parse_args()
    main(args)

'''
Example Call
python DysonOnAxis.py G:/Robert/QuickCal/data/38 
--ctd_file=G:/DY2408/calibration/CalibrationCTD/000_processed.cnv 
--sphere_range=21 --freq=38000
'''