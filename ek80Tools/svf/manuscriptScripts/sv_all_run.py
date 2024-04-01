from glob import glob
import svf_buildAll as svf
from echolab2.instruments import EK80
import pandas as pd
import FullWalk as fw
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

dataDir = 'E:/BB/202204_Shelikof/EK80/' #The data directory 
fl = fw.buildFileLists(dataDir) #Build the file lists from each day with a max file count of 5

frequencies = [38000,70000,120000,200000] # set the frequencies to read in
exclude_above = 5 # this can be a pyEcholab line object or an int/flt, I'm giving it 5 m as a cutoff for the whole file
exclude_below = 'xyz' # I'm going to grab the xyz files and use the line object in the file for the bottom

output_dir = 'D:/FMCW/svf/data/output/' # Set the output directory
passiveLookup = pd.read_pickle('D:/FMCW/svf/passiveLookup.pkl') # To make things easier, we've already built the passive lookup table

for files in tqdm(fl):
    ek80 = EK80.EK80()
    ek80.read_raw(files)

    svf_results = svf.svf()
    svt_results_fm = svf.svt()
    svt_results_cw = svf.svt(pulse='CW')

    # Define the frequency and channel id
    for freq in frequencies:
        inputs = svf.inputs(ek80,freq,add_cw=True,frequency_resolution=None)
        inputs.get_bottom_xyz() # Get the bottom line from the xyz files and add it to our inputs
        svf_results.calc_sample_Svf(inputs,exclude_above_line=5,exclude_below_line='xyz') 
        svf_results.grid_Svf(inputs,interval_length=50, layer_thickness=5)
        svf_results.get_noise(inputs)
        svf.write_grid_to_csv(svf_results,freq,output_dir=output_dir+'svf/')
        
        svt_results_fm.grid_Svt(inputs, interval_length=50, layer_thickness=5,exclude_below_line='xyz', exclude_above_line=5,new_gain='G_fc_fm')
        svt_results_fm.get_noise(inputs, passive_lookup=passiveLookup)
        svf.write_grid_to_csv(svt_results_fm,freq,output_dir=output_dir+'svt/fc/')

        svt_results_fm.grid_Svt(inputs, interval_length=50, layer_thickness=5,exclude_below_line='xyz', exclude_above_line=5,new_gain='G_fave_fm')
        svf.write_grid_to_csv(svt_results_fm,freq,output_dir=output_dir+'svt/fave/')

        svt_results_fm.grid_Svt(inputs, interval_length=50, layer_thickness=5,exclude_below_line='xyz', exclude_above_line=5,new_gain='G_int_fm')
        svf.write_grid_to_csv(svt_results_fm,freq,output_dir=output_dir+'svt/int/')


        svt_results_cw.grid_Svt(inputs, interval_length=50, layer_thickness=5,exclude_below_line='xyz', exclude_above_line=5)
        svt_results_cw.get_noise(inputs, passive_lookup=passiveLookup)
        svf.write_grid_to_csv(svt_results_cw,freq,output_dir=output_dir+'cw/')