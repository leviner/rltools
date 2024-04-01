from glob import glob
import os
import numpy as np
from echolab2.instruments import EK80
from tqdm.notebook import tqdm
from scipy import interpolate
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import pandas as pd

from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
from echolab2.processing import mask, line, grid, afsc_bot_detector, integration

# This function will build a list of files to process based on the dates in the passive lookup table
def buildFileLists(dataDir):
    files = sorted(glob(dataDir+'*.raw'))
    goodFiles = []
    # This is to sort through the raw files and find the ones with xyz
    for file in files:
        if glob(dataDir+os.path.basename(file)[:-4]+'-ES38*-1.xyz'):
            goodFiles.append(file)
    # This is to sort through the actual data
    fileDates = [(file.split('-')[1]) for file in goodFiles]

    uniqueDates = sorted(list(set(fileDates)))
    # We're only going to use days where there is passive data
    # uniqueDates = sorted(list(set(passiveLookup.index.get_level_values('Date'))))# np.unique(fileDates)

    # To make processing a little easier, I'm gonna make the file groups into shorter lists
    n=5 # for now a max of 10 files
    fileLists = []
    for curDate in uniqueDates: 
        curFileInd= [i for i,x in enumerate(fileDates) if x==curDate]
        if len(curFileInd) > n:
            curFileInd = [curFileInd[i:i + min(n, int(round((len(curFileInd)/2))))] for i in range(0, len(curFileInd), min(n, int(round((len(curFileInd)/2)))))]
            for l in curFileInd:
                fileLists.append([goodFiles[k] for k in l])
        else:
            fileLists.append([goodFiles[k] for k in curFileInd])
    print('made', len(fileLists),'groups of files')
    return fileLists # This is a list of lists of files to process

# Integrator function, input is a list of raw files, the passive lookup table, the frequency,
# and whether or not to save the data. If you choose to save, there will be a csv for each file
# group, each frequency, and each pulse type. The csvs will be saved in the analysisFiles folder.
# We also save a csv of the calculated noise curve for each frequency and pulse type.

def integrate(rawfiles,passiveLookup,dfCal, calType='G_fc_fm',f=38000,save=False,fNo=0,outDir='output'): # add bin size as an argument
    intData = {}
    curDate = rawfiles[0].split('-')[1]

    if calType=='G_fave_fm':
        psiOffset = {38000:0.085,70000:.37,120000:.258,200000:.223}
    elif calType=='G_fc_fm':
        psiOffset = {38000:0.0,70000:0.0,120000:0.0,200000:0.0}
    elif calType=='G_int_fm':
        psiOffset = {38000:0.085,70000:.37,120000:.258,200000:.223}
    
    ek80 = EK80.EK80()
    ek80.read_raw(rawfiles,frequencies=[f])
    channels = ek80.frequency_map[f]

    for c in channels:
        curChannel = ek80.get_channel_data()[c][0]
        for file in rawfiles:
            try:
                lineFile = glob(os.path.dirname(file)+'/'+os.path.basename(file)[:-4]+'*'+c.split(' ')[-1][:-2]+'*'+c.split(' ')[-1][-1:]+'.xyz')
                l = line.read_xyz(lineFile[0])
                l.add_object_attribute('data_type','line')
                l.add_object_attribute('frequency',0)
                if file == rawfiles[0]:
                    lAll = l
                else:
                    lAll.append(l)
                xyzExist = True
            except:
                xyzExist = False

        cal = curChannel.get_calibration()
        if curChannel.is_cw():
            cal.sa_correction = np.full(len(cal.sa_correction),0)
            cal.gain = np.full(len(cal.gain),10*np.log10(np.nanmean(10**(dfCal[dfCal.f==f/1000]['G_cw']/10))))
            try:
                noise = passiveLookup.loc[(curDate,f,True)].noiseMean
            except:
                noise = np.percentile(passiveLookup.xs(f,level=1).xs(True,level=1).noiseMean,.90)
        else:
            cal.sa_correction = np.full(len(cal.sa_correction),0)
            cal.gain = np.full(len(cal.gain),10*np.log10(np.nanmean(10**(dfCal[dfCal.f==f/1000][calType]/10))))
            cal.equivalent_beam_angle = cal.equivalent_beam_angle+psiOffset[f]
            try:
                noise = passiveLookup.loc[(curDate,f,False)].noiseMean
            except:
                noise = np.percentile(passiveLookup.xs(f,level=1).xs(False,level=1).noiseMean,.90)
        
        Sv = curChannel.get_Sv(calibration=cal)
        SvNoise = noise+(20*np.log10(Sv.range))+(2*np.unique(cal.absorption_coefficient)*Sv.range)
        SvRange = Sv.range

        Sv.to_depth()
        Sv.delete(start_ping=Sv.n_pings-((Sv.n_pings - lAll.n_pings)-1), end_ping=Sv.n_pings)
        
        g = grid.grid(interval_length=50, interval_axis='ping_number',data=Sv, layer_thickness=5)
        i = integration.integrator(min_threshold_applied=False)
        
        if xyzExist:
            lAll.ping_time = Sv.ping_time
            lAll.data = lAll.data-5
            lAll.data[lAll.data>Sv.depth.max()] = Sv.depth.max()-5
            lAll.data[lAll.data<0]=lAll.data.max()
            if len(lAll.data) > len(lAll.ping_time):
                lAll.data = lAll.data[:len(lAll.ping_time)]
            surf_line = line.line(ping_time=Sv.ping_time, data=15) # generate a line at 15 m depth 
            i_2 = i.integrate(Sv,g,exclude_above_line=surf_line,exclude_below_line=lAll)
        else:
            surf_line = line.line(ping_time=Sv.ping_time, data=15) # generate a line at 15 m depth 
            i_2 = i.integrate(Sv,g,exclude_above_line=surf_line)
        
        i_2.mean_Sv[i_2.mean_Sv==-999]=np.nan
        i_2.nasc[i_2.mean_Sv==-999]=np.nan
        i_2.noise = SvNoise
        i_2.Sv_range = SvRange


        if save:
            i_2.export_to_csv(outDir+'/'+calType+'/'+curDate+'_'+str(fNo)+'_'+str(np.unique(curChannel.pulse_form)[0])+'_intCal_'+str(int(f))+'.csv')
            if fNo == 0:
                noise_df = pd.DataFrame({'range':SvRange,'noise':SvNoise})
                noise_df.to_csv(outDir+'/'+calType+'/'+curDate+'_'+str(np.unique(curChannel.pulse_form)[0])+'_'+str(int(f))+'__intCal_noise.csv',index=False)
        
        if curChannel.is_fm():
            intData['FM'] = i_2
        if curChannel.is_cw():
            intData['CW'] = i_2
    return intData

# This reads in the csv files for a given day, frequency, and pulse type and returns the integration data by cell along with the SNR. 
# This is for the original datatype and is no longer used.
def csv2SNR(date,f,pulseform,outDir='analysisFiles/'):
    
    def layerNoise(startDepth,endDepth):
        return noise[(noise['range'] >= startDepth) & (noise['range'] <= endDepth)]['noise'].mean()

    allIntData = []
    for file in glob(outDir+date+'*'+str(pulseform)+'_intCal_'+str(f)+'.csv'):
        curdf = pd.read_csv(file)
        curdf['fNo'] = file.split('_')[1]
        allIntData.append(curdf)
    intData = pd.concat(allIntData)
    
    noise = pd.read_csv(outDir+date+'_'+str(pulseform)+'_'+str(f)+'__intCal_noise.csv')
    intData['noise'] = intData.apply(lambda row : layerNoise(row['layer_start'],row['layer_end']), axis = 1)
    intData['SNR'] = intData['mean_Sv']-intData['noise']
    intData['time_middle'] = pd.to_datetime(intData.time_middle)
    intData['time_start'] = pd.to_datetime(intData.time_start)
    intData['time_end'] = pd.to_datetime(intData.time_end)

    return intData

