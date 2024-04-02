import Calculation as calc
import numpy as np
import gains
from scipy.interpolate import interp1d
from echolab2.processing import  grid
import csv
from echolab2.processing import  line, integration
from glob import glob
import os


class inputs(object):
    '''
    This sets up an 'inputs' object that takes care of all the various variables needed for 
    the sv(f) calculations. It's a bit of a catch all, but it's a good way to keep everything
    together for a single frequency.
    '''

    def __init__(self, ek80,freq,add_cw=False, n_f_points=500):

        
        self.n_f_points = n_f_points # the number of frequency steps (Lars uses 1000), this is what will be used for the fft

        gain = gains.gains() # This should be replaced with a pyEcholab grab of xml files(s)
        self.gaindf = gain.df # This is also a future delme

        fm_channel_ids = self.get_fm_channel_ids(ek80) # get the fm channels
        self.data = ek80.get_channel_data()[fm_channel_ids[freq]][0] # get the data object for the frequency we're working on
        self.raw_files = np.unique([k['file_path']+'\\'+k['file_name'] for k in self.data.configuration])

        self.fnom = freq
        self.cal = self.data.get_calibration() # grab pyEcholab calibration params

        # Use the little gains class built for this project to get the frequency and gain vectors 
        try:
            self.f,self.gf = gain.getGains(int(self.fnom/1000),1)
        except:
            print('No gain table, using parameters as stored in header')
            self.f,self.gf = self.cal
            

        # If we're adding the CW data, we need to grab the data and the calibration and just process it up front
        if add_cw:
            f_cw,self.gf_cw = gain.getGains(int(self.fnom/1000),0)
            cw_channel_ids = self.get_cw_channel_ids(ek80)
            self.data_cw = ek80.get_channel_data()[cw_channel_ids[freq]][0]
            self.cal_cw = self.data_cw.get_calibration()
            self.cal_cw.gain = np.array([self.gf_cw] *len(self.cal_cw.gain))
            self.cal_cw.sa_correction = np.array([0]*len(self.cal_cw.sa_correction))
            self.Sv_cw = self.data_cw.get_Sv(calibration = self.cal_cw)
        
        self.get_vars() # get all the input params needed for the calulculation functions
        self.build_transmit_signal() # make our transmit signal for match filtering
        self.Sv_t = self.data.get_Sv(calibration = self.cal) # get the Sv(t) data 
        
        # If we have gps, grab it 
        try:
            self.nmea = ek80.nmea_data.interpolate(self.data, 'RMC')
        except:
            print('No NMEA data found')

    def get_fm_channel_ids(self,ek80):
        '''
        Makes a list of the channel ids where the pulse is FM.
        !!!This should be updated now that there is an 'is_cw()' property in the pyEcholab data object!!!
        '''
        fnom_dict = {34000.0:38000,50000.0:70000,95000.0:120000,165000.0:200000}
        fm_channel_ids = {}
        for channel in ek80.channel_ids:
            if ek80.get_channel_data()[channel][0].pulse_form[0]:
                fm_channel_ids[fnom_dict[ek80.get_channel_data()[channel][0].frequency_start[0]]]=channel
        return fm_channel_ids
    
    def get_cw_channel_ids(self,ek80):
        '''
        Makes a list of the channel ids where the pulse is CW.
        !!!This should be updated now that there is an 'is_cw()' property in the pyEcholab data object!!!
        '''
        cw_channel_ids = {}
        for channel in ek80.channel_ids:
            if ek80.get_channel_data()[channel][0].pulse_form[0]==0:
                cw_channel_ids[ek80.get_channel_data()[channel][0].frequency[0]]=channel
        return cw_channel_ids
    
    def get_vars(self):
        '''
        Grab all the required variables from the data to pass to the calculations library functions. Everything
        needed for Lars' functions exists in pyEcholab, it's just organized a bit differently.
        '''
      
        # transducer and receiver impedance
        self.z_td_e, self.z_rx_e = self.data.ZTRANSDUCER, self.data.ZTRANSCEIVER

        # We can go straight to the final sampling frequency after all the decimation by just taking it from the samples directly
        self.f_s_dec = 1/self.data.sample_interval[0]

        # Pulse duration from the header
        self.tau = self.data.pulse_duration[0]

        # Pulse slope from the header
        self.slope = self.data.slope[0]

        # Grab and restructure the filter stages
        self.filter_v = []
        for fk in [1,2]:
            self.filter_v.append({'h_fl_i':self.data.filters[0][fk]['coefficients'],'D':self.data.filters[0][fk]['decimation_factor']})

        # Decimated sampling interval from the header
        self.f_s_dec = np.round(1/self.data.sample_interval[0])

        # Let's back out the undecimated sampling frequency
        self.f_s = self.f_s_dec*self.filter_v[1]['D']*self.filter_v[0]['D']

        # Number of channels
        self.N_u = len(self.data.complex[0][0])

    # Grab the frequency range
        self.f_0, self.f_1 = self.data.frequency_start[0], self.data.frequency_end[0]
        f_n = (self.data.frequency_start+self.data.frequency_end)[0]/2
        self.f_m = np.ravel(np.linspace(self.f_0, self.f_1, self.n_f_points))

        psi_f_c = 10**(self.cal.equivalent_beam_angle[0]/10)
        

        # Get the range for every sample
        self.r_n = self.get_range_vector(self.data)

        # Get the delta range
        self.dr = self.r_n[1]-self.r_n[0]

        # grab c and the transmit power
        self.c = self.data.sound_velocity[0]
        self.p_tx_e = self.data.transmit_power[0]

        self.f_m = np.ravel(np.linspace(self.f_0, self.f_1, self.n_f_points))

        cal_f = np.append(self.f_m[np.where(self.f_m < self.f[0])],self.f)
        cal_g = np.append(np.array([self.gf[0]] * len(self.f_m[np.where(self.f_m < self.f[0])])), self.gf)
        cal_f = np.append(cal_f,self.f_m[np.where(self.f_m > cal_f[-1])])
        cal_g = np.append(cal_g, np.array([self.gf[-1]] * len(self.f_m[np.where(self.f_m > self.f[-1])])))
        cal_g = 10**(cal_g/10) # make it linear

        f_interp = interp1d(cal_f, cal_g,fill_value=np.nan)

        self.g_0_m = f_interp(self.f_m)

        self.alpha_m =  calc.Calculation.calcalpha(self.data,self.f_m)
        self.psi_m = calc.Calculation.calcpsi(psi_f_c, f_n, self.f_m)
        self.lambda_m = self.c / self.f_m

    def build_transmit_signal(self):
        # build the transmit signal
        y_tx_n, t = calc.Calculation.generateIdealWindowedTransmitSignal(self.f_0, self.f_1, self.tau, self.f_s, self.slope)
        # Normalize it (it should technically already be correct)
        y_tilde_tx_n = calc.Calculation.calcNormalizedTransmitSignal(y_tx_n)
        # Use the filters to decimate the signal
        y_tilde_tx_nv = calc.Calculation.calcFilteredAndDecimatedSignal(y_tilde_tx_n, self.filter_v)
        # Use only the fully decimated (final) signal for our match filter
        self.y_mf_n = y_tilde_tx_nv[-1]
        # Determine the autocorrelation of the match filter
        self.y_mf_auto_n, tau_eff =  calc.Calculation.calcAutoCorrelation(self.y_mf_n, self.f_s_dec)

    def get_bottom_xyz(self):
        '''
        Make a pyEcholab line object from an xyz file
        '''
        c_id = self.data.channel_id.split(' ')[-1]
        self.bottom_xyz = self.__get_bottom_xyz(c_id, Sv=self.Sv_t)
        try:
            self.data_cw
            c_id = self.data_cw.channel_id.split(' ')[-1]
            self.bottom_xyz_cw = self.__get_bottom_xyz(c_id, Sv=self.Sv_cw)
        except:
            pass
        
        
    def __get_bottom_xyz(self, c_id, Sv):
        raw_files_datetime = [os.path.basename(r).split('_')[1][5:-4] for r in self.raw_files]
        line_files = []
        for dt in raw_files_datetime:
            #print(os.path.dirname(self.raw_files[0])+'/*'+dt+'*'+c_id.split('-')[0]+'*'+c_id.split('_')[1]+'.xyz')
            line_files.append(glob(os.path.dirname(self.raw_files[0])+'/*'+dt+'*'+c_id.split('-')[0]+'*'+c_id.split('_')[1]+'.xyz')[0])
        if len(line_files)!=len(self.raw_files):
            print('Number of line files does not match number of raw files')

        for lf in line_files:
            l = line.read_xyz(lf)
            l.add_object_attribute('data_type','line')
            l.add_object_attribute('frequency',0)
            if lf == line_files[0]:
                lAll = l
            else:
                lAll.append(l)

        lAll.ping_time = Sv.ping_time
        lAll.data = lAll.data-lAll.transducer_draft-0.5 # added half meter
        lAll.data[lAll.data>Sv.range.max()] = Sv.range.max() - 5
        lAll.data[lAll.data<0]=lAll.data.max()
        if len(lAll.data) > len(lAll.ping_time):
            lAll.data = lAll.data[:len(lAll.ping_time)]
        if lAll.data.shape[0] != Sv.n_pings:
            print('XYZ bottom line does not match length of Sv data, proceed with caution')
            if lAll.data.shape[0] > Sv.n_pings:
                print('Trimming XYZ bottom line to match Sv data')
                lAll.data = lAll.data[:Sv.n_pings]
            elif lAll.data.shape[0] < Sv.n_pings:
                print('Padding XYZ bottom line to match Sv data')
                lAll.data = np.append(lAll.data,np.full(Sv.n_pings-len(lAll.data),lAll.data[-1]))
        
        lAll.data = lAll.data - 5
        return lAll
    
    @staticmethod
    def get_range_vector(data):
        """
        get_range_vector returns a non-corrected range vector.
        """
        # Calculate the thickness of samples with this sound speed.
        thickness = data.sample_interval[0] * data.sound_velocity[0] / 2.0
        # Calculate the range vector.
        range = (np.arange(0, data.n_samples) + data.sample_offset[0]) * thickness
        range[0] = 1e-20

        return range
    

    
class svf(object):

    '''
    Sv(f) class for FM only
    
    '''
        
    def __init__(self):
        self.Sv_sample = {}
        self.frequency = {}
        self.Sv_grid = {}
        self.Sv_noise = {}
        self.ping_time = {}
        self.svf_range = {}


    def calc_sample_Svf(self,inputs,exclude_above_line=None, exclude_below_line=None,Nfft=None, step=None,ping_start=None,ping_end=None):
        '''
        This is just a class wrapper for the calc_sample_Svf function below
        '''
        Sv_sample, frequency, svf_range, ping_time = calc_sample_Svf(inputs,exclude_above_line=exclude_above_line, exclude_below_line=exclude_below_line,Nfft=Nfft,step=step,ping_start=ping_start,ping_end=ping_end)
        self.Sv_sample[inputs.fnom] = Sv_sample
        self.frequency[inputs.fnom] = frequency
        self.svf_range[inputs.fnom] = svf_range
        self.ping_time[inputs.fnom] = ping_time


    def grid_Svf(self,inputs,interval_length=50, layer_thickness=5):
        '''
        Grid up the Sv(f) by sample data according to the input interval length and layer thickness.
        '''

        svf_grid = []

        self.g = grid.grid(interval_length=interval_length, interval_axis='ping_number',layer_axis='range',data=inputs.Sv_t, layer_thickness=layer_thickness)

        for iter_interval in range(self.g.n_intervals):
            #self.calc_noise(iter_interval)
            hold_layer = []
            for iter_layer in range(self.g.n_layers):
                layer_i = np.where([(self.svf_range[inputs.fnom]>=self.g.layer_edges[iter_layer])&(self.svf_range[inputs.fnom]<self.g.layer_edges[iter_layer]+self.g.layer_thickness)])[1]
                cur_cell_full = self.Sv_sample[inputs.fnom][self.g.ping_start[iter_interval]:self.g.ping_end[iter_interval]+1,layer_i[0]:layer_i[-1]+1]
                hold_layer.append(pMean(pMean(cur_cell_full,axis=1),axis=0))
            svf_grid.append(np.array(hold_layer))
        self.Sv_grid[inputs.fnom] = np.array(svf_grid)
        
        if inputs.nmea:
            self.g.interval_latitude_edges = inputs.nmea[1]['latitude'][self.g.interval_edges.astype(int)[:-1]]
            self.g.interval_longitude_edges = inputs.nmea[1]['longitude'][self.g.interval_edges.astype(int)[:-1]]

        self.first_raw_datetime = [os.path.basename(r).split('_')[1][5:-4] for r in inputs.raw_files][0]
            


    def get_noise(self,inputs):
        '''
        Follow De Robertis and Higgenbottom 2007 noise estimation, will work on CW or FM data
        '''
        Sv_sample_noise, frequency, svf_range, ping_time = calc_sample_Svf(inputs,ping_start=self.g.ping_start[0],ping_end=self.g.ping_start[0]+1)
        hold_layer_noise = []
        for iter_layer in range(self.g.n_layers):
            layer_i = np.where([(self.svf_range>=self.g.layer_edges[iter_layer])&(self.svf_range<self.g.layer_edges[iter_layer]+self.g.layer_thickness)])[1]
            cur_cell_full = Sv_sample_noise[:,layer_i[0]:layer_i[-1]+1]
            hold_layer_noise.append(pMean(pMean(cur_cell_full,axis=1),axis=0)\
                                    - (20*np.log10(self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))) \
                                    - (2*inputs.alpha_m*(self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))))
        Sv_noise = np.nanmin(hold_layer_noise,axis=0)

        Sv_noise_bylayer = []
        for iter_layer in range(self.g.n_layers):
            Sv_noise_bylayer.append(Sv_noise+(20*np.log10((self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))))\
             +(2*np.unique(inputs.alpha_m)*(self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))))
        
        self.Sv_noise[inputs.fnom] = np.array(Sv_noise_bylayer)



class svt(object):
    '''
    Sv(t) class for both FM and CW
    '''
    def __init__(self,pulse='FM'):
        self.Sv_sample = {}
        self.frequency = {}
        self.Sv_grid = {}
        self.sA_grid = {}
        self.Sv_noise = {}
        self.pulse = pulse

    def grid_Svt(self,inputs,interval_length=50, layer_thickness=5, exclude_below_line=None, exclude_above_line=None,new_gain=None):

        self.frequency[inputs.fnom] = inputs.fnom

        if self.pulse=='FM':
            if exclude_below_line: # something exists, either it's true, an 'xyz' or a value
                if isinstance(exclude_below_line,int) | isinstance(exclude_below_line,float): # if it's a float or an int, just make it a vector of that value
                    exclude_below_line = line.line(ping_time=inputs.data.ping_time, data=exclude_below_line) 
                elif exclude_below_line == 'xyz': # if it's the xyz, then lets grab the data
                        exclude_below_line = inputs.bottom_xyz
            if exclude_above_line: # something exists, either it's true, an 'xyz' or a value
                if isinstance(exclude_above_line,int) | isinstance(exclude_above_line,float): # if it's a float or an int, just make it a vector of that value
                        exclude_above_line = line.line(ping_time=inputs.data.ping_time, data=exclude_above_line)
                elif exclude_above_line == 'xyz': # if it's the xyz, then lets grab the data
                        exclude_above_line = inputs.top_xyz    
        elif self.pulse=='CW':
            if exclude_below_line: # something exists, either it's true, an 'xyz' or a value
                if isinstance(exclude_below_line,int) | isinstance(exclude_below_line,float): # if it's a float or an int, just make it a vector of that value
                    exclude_below_line = line.line(ping_time=inputs.data_cw.ping_time, data=exclude_below_line) 
                elif exclude_below_line == 'xyz': # if it's the xyz, then lets grab the data
                        exclude_below_line = inputs.bottom_xyz_cw
            if exclude_above_line: # something exists, either it's true, an 'xyz' or a value
                if isinstance(exclude_above_line,int) | isinstance(exclude_above_line,float): # if it's a float or an int, just make it a vector of that value
                        exclude_above_line = line.line(ping_time=inputs.data_cw.ping_time, data=exclude_above_line)
                elif exclude_above_line == 'xyz': # if it's the xyz, then lets grab the data
                        exclude_above_line = inputs.top_xyz_cw    
        
        if self.pulse=='FM':
            self.gain_type = 'default'
            if new_gain: # This is special just for this project, basically if you want to use a different gain from the gain table...    

                psiOffset = {'G_fave_fm':{38000:0.085,70000:.37,120000:.258,200000:.223},
                'G_fc_fm':{38000:0.0,70000:0.0,120000:0.0,200000:0.0},
                'G_int_fm':{38000:0.085,70000:.37,120000:.258,200000:.223}}

                self.gain_type = new_gain
                new_cal = inputs.cal
                new_cal.gain = np.array([pMean(inputs.gaindf[inputs.gaindf.f == inputs.fnom/1000][new_gain])]*len(new_cal.gain))
                new_cal.equivalent_beam_angle = new_cal.equivalent_beam_angle+psiOffset[new_gain][inputs.fnom]
                Sv = inputs.data.get_Sv(calibration=new_cal)
                self.__run_grid_Svt(Sv,inputs.fnom,interval_length,layer_thickness,exclude_above_line, exclude_below_line)
            else:    
                self.__run_grid_Svt(inputs.Sv_t,inputs.fnom,interval_length,layer_thickness,exclude_above_line, exclude_below_line)

        elif self.pulse=='CW':
            self.__run_grid_Svt(inputs.Sv_cw,inputs.fnom,interval_length,layer_thickness,exclude_above_line, exclude_below_line)

        if inputs.nmea:
            self.g.interval_latitude_edges = inputs.nmea[1]['latitude'][self.g.interval_edges.astype(int)[:-1]]
            self.g.interval_longitude_edges = inputs.nmea[1]['longitude'][self.g.interval_edges.astype(int)[:-1]]

        self.first_raw_datetime = [os.path.basename(r).split('_')[1][5:-4] for r in inputs.raw_files][0]
    
    def __run_grid_Svt(self, Sv,f,interval_length,layer_thickness,exclude_above_line, exclude_below_line):

        i = integration.integrator(min_threshold_applied=False)
        self.g = grid.grid(interval_length=interval_length, interval_axis='ping_number',layer_axis='range',data=Sv, layer_thickness=layer_thickness)
        i_2 = i.integrate(Sv,self.g,exclude_above_line=exclude_above_line,exclude_below_line=exclude_below_line)
        i_2.mean_Sv[i_2.mean_Sv==-999]=np.nan
        i_2.nasc[i_2.mean_Sv==-999]=np.nan
        self.Sv_grid[f] = i_2.mean_Sv
        self.sA_grid[f] = i_2.nasc


    def get_noise(self, inputs, passive_lookup=None):
        '''
        If a lookup table is provided, use the lookup estimate of noise, otherwise use the De Robertis and Higgenbottom 2007 noise estimation
        '''
        curDate = np.unique([os.path.basename(r).split('_')[1][5:-4].split('-')[0] for r in inputs.raw_files])[0]
        if passive_lookup is not None:
            if self.pulse== 'FM':
                Sv_noise = self.__get_noise_lookup(inputs.fnom, curDate, passive_lookup, False)
            if self.pulse== 'CW':
                Sv_noise = self.__get_noise_lookup(inputs.fnom, curDate, passive_lookup, True)
        else:
            hold_layer_noise = []
            for iter_layer in range(self.g.n_layers):
                layer_i = np.where([(inputs.Sv_t.range>=self.g.layer_edges[iter_layer])&(inputs.Sv_t.range<self.g.layer_edges[iter_layer]+self.g.layer_thickness)])[1]
                cur_cell_full = inputs.Sv_t[:,layer_i[0]:layer_i[-1]+1]
                hold_layer_noise.append(pMean(pMean(cur_cell_full,axis=1),axis=0)\
                                        - (20*np.log10(self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))) \
                                        - (2*inputs.cal.absorption_coefficient[0]*(self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))))
            Sv_noise = np.nanmin(hold_layer_noise,axis=0)
        Sv_noise_bylayer = []
        for iter_layer in range(self.g.n_layers):
            Sv_noise_bylayer.append(Sv_noise+(20*np.log10((self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))))\
            +(2*(inputs.cal.absorption_coefficient[0])*(self.g.layer_edges[iter_layer]+(self.g.layer_thickness/2))))
        
        self.Sv_noise[inputs.fnom] = np.array(Sv_noise_bylayer)


    def __get_noise_lookup(self, f, curDate, passive_lookup, CW):
        try:
            noise = passive_lookup.loc[(curDate,f,CW)].noiseMean
        except:
            noise = np.percentile(passive_lookup.xs(f,level=1).xs(CW,level=1).noiseMean,.90)
        return noise



'''
Functions for general processing
'''

# Write the gridded sv(f) data to a csv file matching the echoview exports
def write_grid_to_csv(results,freq,output_dir=''):

    with open(output_dir+results.first_raw_datetime+'_'+str(int(freq/1000))+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')

        writer.writerow(['channel','interval','layer','time_start','time_end','ping_start','ping_end',
                            'lat_start','lon_start','layer_start','layer_end','layer_thickness','frequency','mean_Sv','SNR'])
        
        for i_int in range(len(results.Sv_grid[freq])):
            for i_layer in range(len(results.Sv_grid[freq][i_int])):
                # put noise by layer in here
                if isinstance(results.frequency[freq],int): # CW Data
                    writer.writerow([freq,i_int+1,i_layer+1,
                                results.g.time_start[i_int],
                                results.g.time_end[i_int],
                                results.g.ping_start[i_int],
                                results.g.ping_end[i_int],
                                results.g.interval_latitude_edges[i_int],
                                results.g.interval_longitude_edges[i_int],
                                results.g.layer_edges[i_layer], results.g.layer_edges[i_layer]+results.g.layer_thickness,
                                results.g.layer_thickness,
                                results.frequency[freq],
                                results.Sv_grid[freq][i_int][i_layer],
                                (results.Sv_grid[freq][i_int][i_layer]-results.Sv_noise[freq][i_layer])])
                else:
                    writer.writerow([freq,i_int+1,i_layer+1, # FM Data
                                results.g.time_start[i_int],
                                results.g.time_end[i_int],
                                results.g.ping_start[i_int],
                                results.g.ping_end[i_int],
                                results.g.interval_latitude_edges[i_int],
                                results.g.interval_longitude_edges[i_int],
                                results.g.layer_edges[i_layer], results.g.layer_edges[i_layer]+results.g.layer_thickness,
                                results.g.layer_thickness,
                                ';'.join([str(j) for j in results.frequency[freq]]),
                                ';'.join([str(j) for j in (results.Sv_grid[freq][i_int][i_layer])]),
                                ';'.join([str(j) for j in (results.Sv_grid[freq][i_int][i_layer]-results.Sv_noise[freq][i_layer])])])
            

# power mean
def pMean(data,axis=0):
    return 10*np.log10(np.nanmean(10**(data/10),axis=axis))

# power standard deviation
def pStd(data,axis=0):
    return 10*np.log10(np.nanmean(10**(data/10),axis=axis)-np.std(10**(data/10),axis=axis)),10*np.log10(np.nanmean(10**(data/10),axis=axis)+np.std(10**(data/10),axis=axis))


# This is the main function that does the work for calculating Svf for each sample
def calc_sample_Svf(svf_inputs,exclude_above_line=None, exclude_below_line=None,Nfft=None,step=None,ping_start=None,ping_end=None):
    if step:
        step = step
    else:
        # This can be coded up so that we can do it in m instead of in samples based on the pulse length in the file
        step_size = {38000:10,70000:43,120000:64,200000:128} #roughly correspond to ~.5 m using a 1ms pulse
        step = step_size[svf_inputs.fnom]
    
    if Nfft:
        Nfft = Nfft
    else:
        # This can be coded up so that we can do it in m instead of in samples based on the pulse length in the file
        Nfft_size = {38000:42,70000:170,120000:256,200000:512} #roughly correspond to ~2 m using a 1ms pulse
        Nfft = Nfft_size[svf_inputs.fnom]
    
    if not ping_start:
        ping_start = 0
    if not ping_end:
        ping_end = svf_inputs.data.n_pings

    # Grab the complex data, but because of how Lars' code handle it, we have to reformat it such that instead of 
    # n-ping x 4, we have 4 x n-ping where each channel is a vector of len(n-ping)
    #y_rx_nu = data.complex[0]
    #y_rx_nu = np.array([np.array([k[0] for k in y_rx_nu]), np.array([k[1] for k in y_rx_nu]),np.array([k[2] for k in y_rx_nu]),np.array([k[3] for k in y_rx_nu])])
    Sv_full = []
    pTime = []

    # This should be match/case but I'm running 3.9
    if exclude_below_line: # something exists, either it's true, an 'xyz' or a value
        if isinstance(exclude_below_line,int) | isinstance(exclude_below_line,float): # if it's a float or an int, just make it a vector of that value
            exclude_below_line = np.full(len(range(ping_start,ping_end)),fill_value=exclude_below_line)
        elif exclude_below_line == 'xyz': # if it's the xyz, then lets grab the data
                exclude_below_line = svf_inputs.bottom_xyz.data
    if exclude_above_line: # something exists, either it's true, an 'xyz' or a value
        if isinstance(exclude_above_line,int) | isinstance(exclude_above_line,float): # if it's a float or an int, just make it a vector of that value
                exclude_above_line = np.full(len(range(ping_start,ping_end)),fill_value=exclude_above_line)
        elif exclude_above_line == 'xyz': # if it's the xyz, then lets grab the data
                exclude_above_line = svf_inputs.top_xyz.data    
                            

    for ping_no in range(ping_start,ping_end):
        #if you want to an individual ping
        y_rx_nu = svf_inputs.data.complex[ping_no]
        y_rx_nu = np.array([np.array([k[0] for k in y_rx_nu]), np.array([k[1] for k in y_rx_nu]),np.array([k[2] for k in y_rx_nu]),np.array([k[3] for k in y_rx_nu])])

        if exclude_below_line is not None:
            for c in range(len(y_rx_nu)):
                y_rx_nu[c] = np.where((exclude_below_line[ping_no]<svf_inputs.Sv_t.range), np.nan,y_rx_nu[c])
        if exclude_above_line is not None:
            for c in range(len(y_rx_nu)):
                y_rx_nu[c] = np.where((exclude_above_line[ping_no]>svf_inputs.Sv_t.range), np.nan,y_rx_nu[c])

        # Do the pulse compression 
        y_pc_nu = calc.Calculation.calcPulseCompressedSignals(y_rx_nu, svf_inputs.y_mf_n)
        # Take the avearage of the four channels
        y_pc_n = calc.Calculation.calcAverageSignal(y_pc_nu)
        # Convert the complex to received power
        p_rx_e_n = calc.Calculation.calcPower(y_pc_n, svf_inputs.z_td_e, svf_inputs.z_rx_e, svf_inputs.N_u)
        # Compensate for the spreading loss
        y_pc_s_n = calc.Calculation.calcPulseCompSphericalSpread(y_pc_n, svf_inputs.r_n)

        # Here we build the hanning window, outputting the weights, length of window, time interval, and time vector 
        w_tilde_i, N_w, t_w, t_w_n = calc.Calculation.defHanningWindow(svf_inputs.c, svf_inputs.tau, svf_inputs.dr, svf_inputs.f_s_dec, N_w=Nfft)


        # Get the discrete fourier transform of the compressed data using the hanning window
        # This ouputs the dft of the pulse compressed data, the dft of the match filter autocorrelation,
        # and the normalized dft of the pulse compressed data
        # The step is the size of the window (vertical) in samples
        Y_pc_v_m_n, Y_mf_auto_m, Y_tilde_pc_v_m_n, svf_range = calc.Calculation.calcDFTforSv(
            y_pc_s_n, w_tilde_i, svf_inputs.y_mf_auto_n, N_w, svf_inputs.f_m, svf_inputs.f_s_dec, svf_inputs.r_n, step=step)

        # From the DFT, we calculate the power
        P_rx_e_t_m_n = calc.Calculation.calcPowerFreqSv(Y_tilde_pc_v_m_n, svf_inputs.N_u, svf_inputs.z_rx_e, svf_inputs.z_td_e)

        # Last step is to calculate Sv(f) in each of the steps
        Sv_m_n = calc.Calculation.calcSvf(
                P_rx_e_t_m_n, svf_inputs.alpha_m, svf_inputs.p_tx_e, svf_inputs.lambda_m, t_w, svf_inputs.psi_m, svf_inputs.g_0_m, svf_inputs.c, svf_range)
        
        Sv_full.append(Sv_m_n)
        pTime.append(svf_inputs.data.ping_time[ping_no])

    Sv_sample= np.array(Sv_full)  
    frequency= svf_inputs.f_m  
    
    return Sv_sample, frequency, svf_range, np.array(pTime)
