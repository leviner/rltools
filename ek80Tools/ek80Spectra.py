import os
from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap
from echolab2.instruments import EK80
from echolab2.plotting.matplotlib import echogram
import numpy as np
from bs4 import BeautifulSoup as bs
from scipy.signal.windows import tukey
from scipy.fft import fft
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class ek80Spectra():
    def __init__(self,windowX, windowZ, rawFiles, calFiles=None, freqs=None):
        
        if windowX == 0:
            self.ek80 = EK80.EK80()
            print('Reading raw files...')
            self.ek80.read_raw(rawFiles[0])
            self.freqs = [freqs] if freqs is not None else list(self.ek80.frequency_map.keys())
            self.rawFiles = rawFiles
            self.windowX = 10**10 #make it something absurd
        else:
            self.ek80 = EK80.EK80()
            print('Reading raw files...')
            self.ek80.read_raw(rawFiles)
            self.freqs = [freqs] if freqs is not None else list(self.ek80.frequency_map.keys())
            self.windowX = int(windowX)
        print('Found ',len(self.freqs), 'frequencies:',self.freqs)
        self.windowZ = windowZ
        self.calFiles = calFiles
        self.spectraDict = {}
        
    def calcSpectra(self):
        for self.freq in self.freqs:
            try:
                self.calFile = [cf for cf in self.calFiles if str(int(self.freq)) in cf][0]
                print('Calculating spectra for ', self.freq)
            except:
                print('Skipping ', str(self.freq),', no cal found.')
                continue
            self.__setCal__()
            
            self.spectraDict[self.freq] = {}
            self.data = self.ek80.get_channel_data(frequencies=self.freq)[self.freq][0]
            self.__setWindow__()
            
            self.spectraDict[self.freq]['rangeBinCenters'] = self.rangeBinCenters
            self.spectraDict[self.freq]['rangeBinSize'] = self.windowZ
            self.spectraDict[self.freq]['windowSampleIndex'] = self.winIndsZ
            self.spectraDict[self.freq]['windowPingStart'] = []
            self.spectraDict[self.freq]['windowPingStartTime'] = []
            self.spectraDict[self.freq]['Sv'] = []
            

            print('Calculating spectra for',len(self.winIndsZ), 'vertical windows in', len(self.nwinX),' horizontal windows')
            
            if self.windowX ==10**10:
                npings=0
                for file in self.rawFiles:
                    self.ek80.read_raw(file)
                    self.data = self.ek80.get_channel_data(frequencies=self.freq)[self.freq][0]
                    winSpec = []
                    for i in range(len(self.winIndsZ)):
                        ftmp, svtmp = self.__specMath__(self.winIndsZ[i],0,self.rangeBinCenters[i])
                        winSpec.append(svtmp)
                    self.spectraDict[self.freq]['Sv'].append(winSpec)
                    self.spectraDict[self.freq]['windowPingStart'].append(npings)
                    self.spectraDict[self.freq]['windowPingStartTime'].append(self.data.ping_time[0])
                    npings = npings+self.data.n_pings
            else:
                for winIndX in self.nwinX:
                    winSpec = []
                    for i in range(len(self.winIndsZ)):
                        ftmp, svtmp = self.__specMath__(self.winIndsZ[i],winIndX,self.rangeBinCenters[i])
                        winSpec.append(svtmp)
                    self.spectraDict[self.freq]['Sv'].append(winSpec)
                    self.spectraDict[self.freq]['windowPingStart'].append(self.windowX*winIndX)
                    self.spectraDict[self.freq]['windowPingStartTime'].append(self.data.ping_time[self.windowX*winIndX])
            
            self.spectraDict[self.freq]['frequency'] = ftmp
            

    def __setWindow__(self):
        self.range = self.data.get_power().range # we need to do this since range isn't carried in the 'raw'
        if self.windowX == 10**10:
            self.nwinX = [0]
        else:
            self.nwinX = np.arange(int(np.ceil(self.data.n_pings/self.windowX)))
        step = self.windowZ # this is in meters
        self.nfft = 2**9 # make this variable, take a window length across the channels, whats the max number of datapoints, thats the minimum fft for all channels
        deltaRange = self.range[1]-self.range[0]
        maxRange = self.range.max()
        maxRangeClipped = maxRange-75 # these shouldn't be hardcoded but it's just for clipping
        minRange = 3 # these shouldn't be hardcoded but it's just for clipping
        self.rangeBins = np.arange(minRange,maxRangeClipped+step,step)
        self.rangeBinCenters = self.rangeBins + step/2
        self.nIndsZ = int(np.ceil(step/deltaRange))
        self.winIndsZ = []
        for h in self.rangeBins:
            self.winIndsZ.append(np.where(abs(self.range-h)==np.min(abs(self.range-h)))[0][0])
        
    def __setCal__(self):
        # Currently the read_ecs in echolab doesn't work so I'm going to do it myself from the xml. I'm not going to spend time on this since
        # when working, echolab will have the cal values stored for each channel already. Currently all of this can only handle one channel because the cals
        # are all independent (i.e., 1 xml file for each channel). The reac_ecs implementation would fix this.
        if self.calFile is None:
            # This is just a place holder since there are some issue in the calibration object
            self.data.get_calibration().equivalent_beam_angle[0]
            self.data.get_calibration().frequency # this is empty? 
            self.data.get_calibration().gain[0] # this is a single value for every ping but needs to be a table
        else:
            name, extension = os.path.splitext(self.calFile)
            print('Grabbing cal data from ',extension, 'file')
            if extension == '.xml':
                content = []
                with open(self.calFile, "r") as file:
                    # Read each line in the file, readlines() returns a list of lines
                    content = file.readlines()
                    # Combine the lines in the list into a string
                    content = "".join(content)
                    bs_content = bs(content, "lxml")

                self.calF = [float(f) for f in bs_content.find('calibrationresults').find('frequency').get_text().split(';')]
                self.calG = [float(f) for f in bs_content.find('calibrationresults').find('gain').get_text().split(';')]
                self.calPsi = float(bs_content.find('equivalentbeamangle').get_text())

            # ECS can go here, but see comments above. 
            
            
    def __specMath__(self,winIndZ,winIndX,rangeBinCenter):
        complexVoltage = np.mean(self.data.complex,axis=2) # get the complext voltage for each sample from the 4 channels
        #curComplex = 
        cvWindow = (np.mean(complexVoltage[self.windowX*(winIndX):self.windowX*(winIndX+1)],axis=0))[winIndZ:(winIndZ+self.nIndsZ)] # grab the complext voltages in the first window
        self.range = self.data.get_power().range
        rWindow = self.data.get_power().range[winIndZ:(winIndZ+self.nIndsZ)] # grab the ranges of that first window
        sVector = cvWindow*rWindow #scale the complex voltage by range to account for spreading 
        b = tukey(self.nIndsZ,0.1)/(np.linalg.norm(tukey(self.nIndsZ,0.1))/np.sqrt(self.nIndsZ)) # build the tukey window with a 10% taper
        sVector = sVector*b # apply the tukey window
        sVector = fft(sVector,self.nfft) # run the fft on the now windowed scaled volltages
        FFTvecin = sVector
        fsdec = 1/self.data.sample_interval[0] # sampling rate, parameter value in the config
        FFTvec_tmp, ftmp = self.freqtransf(FFTvecin,fsdec,self.freq) # this is a frequency transformation   
        env = self.data.environment[self.data.environment != np.array(None)][0]
        self.alpha = [self.alphaFG(env['sound_speed'],env['acidity'],env['temperature'],env['depth'],env['salinity'],nomf/1000)/1000 for nomf in ftmp]
        
        f = interp1d(self.calF, self.calG,fill_value=np.nan) # 1-d interpolation of the calibration gains so that they fit the size of our frequency vector
        ftmp[ftmp< min(self.calF)] = np.nan
        ftmp[ftmp> max(self.calF)] = np.nan
        calPsi = self.calPsi + 20*np.log10(self.freq/ftmp)
        G = f(ftmp) # Gain interpolated to match the frequency values
        dt = 2*(self.range[winIndZ+self.nIndsZ] - self.range[winIndZ])/env['sound_speed']
        
        pr = np.abs(FFTvec_tmp)**2 # note sure I get what this is but I think it's the absolute value of the square of the complex data
        zet = self.data.ZTRANSDUCER # transducer impedance
        zer = self.data.ZTRANSCEIVER # transciever impedance
        pTr = self.data.transmit_power[0] #Transmit power

        svtmp = 10*np.log10(pr) +\
           [(2*rangeBinCenter*a) for a in self.alpha] - 2*G - calPsi - \
           10*np.log10(dt) +\
           10*np.log10(4/zet/pTr/(2*np.sqrt(2))**2) +\
           10*np.log10((zer+zet)/zer) - \
           10*np.log10(env['sound_speed']**3/(32*np.pi**2*ftmp**2))
        return ftmp, svtmp
      
        
    def alphaFG(self,c, pH, T, D, S,f): # requires sound speed (m/s), pH, temp(C), depth(m), salinity(ppt), and nominal frequency(kHz)
        # Attenuation Coefficient is based on Francois and Garrison, 1982 - "Sound absorption based on ocean measurements.
        # Boric Acid Contribution, P1 = 1
        A1=((8.86/c)*(10**(0.78*pH-5)))
        f1=((2.8*((S/35)**0.5))*(10**(4-(1245/(T+273)))))
        # MgSO4 Contribution
        A2=((21.44*(S/c))*(1+(0.025*T)))
        P2=(1-(1.37*(10**-4)*D)+(6.2*(10**-9)*(D**2)))
        f2=((8.17*(10**(8-(1990/(T+273)))))/(1+.0018*(S-35)))
        # Pure water contribution, where A3 is temperature dependent
        if T > 20:
            A3=((3.964*(10**-4))-(1.146*(10**-5)*T)+(1.45*(10**-7)*(T**2))-(6.5*(10**-10)*(T**3)))
        else:
            A3=((4.937*(10**-4))-(2.59*(10**-5)*T)+(9.11*(10**-7)*(T**2))-(1.5*(10**-8)*(T**3)))
        P3=((1-(3.83*(10**-5)*D)) + (4.9*(10**-10)*(D**2)))
        # Calculate and return Alpha
        alpha = (((f**2)*A1*f1)/(((f1**2)) + (f**2)))+ ((A2*P2*f2*(f**2))/((f2**2) + (f**2))) + (A3*P3*(f**2))
        return alpha

    def freqtransf(self,FFTvecin,fsdec,fc): #sampling decimation of band-limited signal
        nfft  = len(FFTvecin)
        FFTvec  = FFTvec = np.ravel([FFTvecin, FFTvecin, FFTvecin])
        fvec  = fsdec*np.linspace(0,1-(1/nfft),nfft)
        fvec   = np.ravel([fvec,fsdec+fvec,(2*fsdec)+fvec])
        if fc > (fsdec/2):
            idxmin      = round((fc-fsdec/2)/fsdec*nfft)-1
        else:
            idxmin      = 0
        FFTvec = FFTvec[idxmin:idxmin+nfft]
        fvec = fvec[idxmin:idxmin+nfft]
        return FFTvec, fvec
    
    
    def plotSpec(self,windowX,windowZ):
        fig = plt.figure(figsize=(10,5))
        for freq in self.spectraDict.keys():
            plt.plot(self.spectraDict[freq]['frequency'],self.spectraDict[freq]['Sv'][windowX][windowZ],color='dimgrey')
            plt.title(self.spectraDict[freq]['windowPingStartTime'][windowX].astype(str)+'\n'+\
                      str(self.spectraDict[freq]['rangeBinCenters'][windowZ]-self.windowZ/2)+\
                      ' to '+str(self.spectraDict[freq]['rangeBinCenters'][windowZ]+self.windowZ/2)+'m range')
        plt.xlabel('Frequency')
        plt.ylabel('Sv')
        plt.show()