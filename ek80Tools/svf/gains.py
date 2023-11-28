import pandas as pd
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import os

class gains():
    def __init__(self,calDataDir = 'E:/BB/WorkingProcessing/data/'):
        self.calDataDir = calDataDir
        self.df = self.buildCalTable(['202200','202204','202207'],['L1','L2','L3'],[38,70,120,200])
        self.df = self.addIntCal(self.df,'E:/BB/WorkingProcessing/calValsIntWorking.csv')

    def getGains(self,f,pulse):
        if pulse == 1:
            files = self.getFMFiles(f)
            sphere, cF, cG = self.calFM(files,f,full=True)
        else:
            cG = 10*np.log10(np.nanmean(10**(self.df[(self.df.f == f)].G_cw/10)))
            cF = f*1000
        return cF, cG

    def getFMFiles(self,f):
        files = []
        cs = {38:'r', 70:'b', 120:'g', 200:'c'}
        cruises, legs, freqs = ['202200','202204','202207'],['L1','L2','L3'],[f]
        for f in freqs:
            for cruise in cruises:
                for leg in legs:
                    for pulse in ['FM']:
                        if glob(self.calDataDir + cruise +'/'+ leg+'/'+str(f)+'*'+pulse+'*.xml'):
                            [files.append(file) for file in glob(self.calDataDir + cruise +'/'+ leg+'/'+str(f)+'*'+pulse+'*.xml')]
        return files

       
    def calFM(self,files,freq,plot=False,full=False):
        if plot:
            plt.figure(figsize=(5,5))
        ranges={38:[34000,45000],70:[50000,90000],120:[95000.0, 155000.0],200:[165000.0, 260000.0]}
        calFfit = np.arange(ranges[freq][0],ranges[freq][1]+1,100)
        if len(files)>1:
            calFs,calGs = [],[]
            allFs,allGs = np.array([]),np.array([])
            ff = True
            for file in files:
                calResults = pd.read_xml(file,xpath=".//CalibrationResults")
                calF, calG = np.array([float(f) for f in calResults.Frequency.values[0].split(';')]), np.array([float(f) for f in calResults.Gain.values[0].split(';')])
                if plot:
                    plt.plot(calF,calG,'.')
                calFs.append(calF)
                calGs.append(calG)
                if ff:
                    allFs = calF
                else:
                    allFs = np.unique(np.concatenate((allFs,calF),0))
                ff = False
            for f in allFs:
                curGains = []
                for gl,fl in zip(calGs,calFs):
                    if gl[fl == f].size > 0:
                        curGains.append((gl[fl == f][0]))
                    else:
                        curGains.append(np.nan)
                allGs = np.append(allGs,10*np.log10(np.nanmean(10**(np.array(curGains)/10))))
            calFfit = np.arange(ranges[freq][0],ranges[freq][1]+1,100)
            calGfit = np.interp(calFfit, allFs,allGs)
            if plot:
                plt.plot(allFs,allGs)
                plt.plot(calFfit,calGfit,':k')
            gFC = calGfit[calFfit == ((ranges[freq][0]+ranges[freq][1])/2)]
            gFM = 10*np.log10(np.nanmean((10**(np.array(calGfit)/10)),axis=0))
            sphere=0
        else:
            calResults = pd.read_xml(files[0],xpath=".//CalibrationResults")
            calF, calG = np.array([float(f) for f in calResults.Frequency.values[0].split(';')]), np.array([float(f) for f in calResults.Gain.values[0].split(';')])
            calFfit = np.arange(ranges[freq][0],ranges[freq][1]+1,100)
            calGfit = np.interp(calFfit, calF,calG)
            gFC = calGfit[calFfit == ((ranges[freq][0]+ranges[freq][1])/2)]
            gFM = 10*np.log10(np.nanmean((10**(np.array(calG)/10)),axis=0))
            sphere=38
        if full:
            return sphere, calFfit, calGfit
        else:
            return sphere, gFC, gFM

    def calCW(self,files):
        if len(files)>1:
            sphere=0
            gCIs = np.array([])
            for file in files:
                calResults = pd.read_xml(file,xpath=".//CalibrationResults")
                gCIs = np.append(gCIs,calResults.Gain+calResults.SaCorrection)
            gCI = 10*np.log10(np.nanmean(10**(np.array(gCIs)/10)))
        else:
            calResults = pd.read_xml(files[0],xpath=".//CalibrationResults")
            gCI = calResults.Gain+calResults.SaCorrection
            sphere=38
        return sphere, gCI

    def buildCalTable(self,cruises, legs, freqs):
        dfsFM,dfsCW = [],[]
        for cruise in cruises:
            for leg in legs:
                for f in freqs:
                    for pulse in ['CW','FM']:
                        files = glob(self.calDataDir+ cruise +'/'+ leg+'/'+str(f)+'*'+pulse+'*.xml')
                        df = pd.DataFrame({'cruise':cruise, 'leg':leg, 'f':f,'nfiles':len(files)},index=[0])
                        if not files:
                            continue
                        if pulse == 'CW':
                            sphere, gCI = self.calCW(files)
                            df['sphere'] = sphere
                            df['G_cw'] = gCI
                            dfsCW.append(df)
                        elif pulse == 'FM':
                            sphere, gFC, gFM = self.calFM(files,f)
                            df['sphere'] = sphere
                            df['G_fc_fm'] = gFC
                            df['G_fave_fm'] = gFM
                            dfsFM.append(df)
        df = (pd.concat(dfsFM)).merge(pd.concat(dfsCW),on=['cruise', 'leg', 'f'],how='outer', suffixes=('_FM', '_CW'))
        return df

    def addIntCal(self,df, intFile):
        dfInt = pd.read_csv(intFile)
        cruises,legs,freqs = ['202200','202204','202207'],['L1','L2','L3'],[38,70,120,200] 
        for cruise in cruises:
                for leg in legs:
                    for f in freqs:
                        curG = dfInt[(dfInt.Cruise == int(cruise))&(dfInt.Leg == leg)&(dfInt.f == f)].G_int_fm
                        df.loc[(df.cruise==cruise)&(df.leg==leg)&(df.f==f),'G_int_fm'] = 10*np.log10(np.nanmean(10**(curG/10)))
        return df


