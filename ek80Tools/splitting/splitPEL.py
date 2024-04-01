# coding=utf-8

from echolab2.instruments import echosounder
from glob import glob
import os
import argparse
from tqdm import tqdm
import numpy as np

class splitFiles():
    '''
    The purpose of this code is to split multichannel EK80 .raw files into individual
    files, each containing the datagram and associated metadata for a  single channel, 
    regardless of pulse form.

    Just for organization I've set this up as a class, but the bare bones of it uses 
    pyEcholab (https://github.com/CI-CMG/pyEcholab) to read in a raw file, identify 
    the channel ids, channel pulse forms, start/end frequencies, and write out a new
    .raw file with an additional suffix to an output directory. Additional wrapping 
    provides the iterating through multiple channels/files.

    usage: splitPEL.py [-h] [--out_dir OUT_DIR] [--read_suffix READ_SUFFIX] in_dir

    positional arguments:
    in_dir             Input directory of raw files

    optional arguments:
    -h, --help           show this help message and exit
    --in_files IN_FILES  define a single or comma-separated list of filenames rather than iterating through all files. 
                            If a filename has whitespace you must use double quotes around it
    --out_dir OUT_DIR    Output directory for split files. Default is new directory 'splitFiles' in input directory
    --read_suffix READ_SUFFIX      File suffix used to identify multichannel/FM file type. E.g., if all target files end in 
                            ..._2.raw, use --read_suffix=_2
    --read_prefix READ_PREFIX     File prefix used to identify multichannel/FM file type. E.g., if all target files start in 
                            DY_FM..., use --read_prefix=DY_FM
    --group GROUP       If True, create new files containing all channels with same pulse form
    --within_channel WITHIN_CHANNEL If True, split files where multiple pulse forms exist in a single channel
    --overwrite OVERWRITE If True, process every file found, regardless of if basename already exists in output
    --ignore_last IGNORE_LAST If True, process every file except the one most recently modified. This is for real-time survey applications.
    
    A major thing to clean up is hitting the 'reset' button between each file. Some 
    of the init assignments need to be done for ever file in a file list so some stuff 
    should get moved around. Lots of cleaning to do.

    todo
    - add write suffix option

'''

    def __init__(self, args):
        self.in_dir = args.in_dir
        if args.out_dir is None: 
           self.out_dir = args.in_dir+'/splitFiles/'
        else:
            self.out_dir = args.out_dir
        if args.in_files:
            self.in_files = [self.in_dir+'/'+f for f in args.in_files.split(',')]
        else:
            self.in_files = glob(self.in_dir+'/'+args.read_prefix+'*'+args.read_suffix+'.raw')
        
        print('Input contains '+str(len(self.in_files))+' files')
        if not args.overwrite:
            for file in sorted(self.in_files):
                if glob(self.out_dir+(os.path.splitext(os.path.basename(file))[0])+'*'):
                    self.in_files.remove(file)
            print(str(len(self.in_files))+' files unprocessed')

        if args.ignore_last:
            latest_file = max(self.in_files, key=os.path.getctime)
            print('Ignoring latest file: '+latest_file+' (most recent modification time)')
            self.in_files.remove(latest_file)

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.group = args.group
        if self.group is True:
            self.channelGroups = {'CW':[],'FM':[]}
        self.within_channel=args.within_channel
        self.main()
        
    def openFile(self):
            self.nameBase = os.path.basename(self.file)
            self.ek_data = echosounder.read([self.file])[0]
    
    def splitFile(self):
        for self.channel in (pbar2 := tqdm(self.ek_data.channel_ids,leave=False)):
            pbar2.set_description('Processing '+str(len(self.ek_data.channel_ids))+' channels')
            self.getchannelFreq()
            if self.group is True:
                self.channelGroups[self.channelPulse].append(self.channel)
            else:
                self.writeChannelIDs = [self.channel]
                self.writeChannels() 
        if self.group:
            for k in self.channelGroups.keys():
                if self.channelGroups[k]:
                    self.writeChannelIDs = self.channelGroups[k]
                    self.groupName = k
                    self.writeChannels()
            self.channelGroups = {'CW':[],'FM':[]}
            

    def setChannelName(self):
        if self.channelPulse == 'CW':
            if self.within_channel:
                for k in self.ria_list['CW']:
                    if self.ria_list['CW'][k].size > 0:
                        f_s =k.frequency
                        break
                self.nameFreq = str(int(f_s[0]/1000))+'CW'
            elif self.group:
                self.nameFreq = self.groupName
            else:
                self.nameFreq = str(int(self.ek_data.get_channel_data()[self.channel][0].frequency[0]/1000))+'CW'
        elif self.channelPulse == 'FM':
            if self.within_channel:
                for k in self.ria_list['FM']:
                    if self.ria_list['FM'][k].size > 0:
                        f_s =k.frequency_start
                        f_e = k.frequency_end
                        break
                self.nameFreq = str(int(f_s[0]/1000))+'-'+str(int(f_e[0]/1000))+'FM'
            elif self.group:
                self.nameFreq = self.groupName
            else:
                self.nameFreq = str(int(self.ek_data.get_channel_data()[self.channel][0].frequency_start[0]/1000))+'-'+\
                        str(int(self.ek_data.get_channel_data()[self.channel][0].frequency_end[0]/1000))+'FM'

    def getchannelFreq(self):
        if self.ek_data.get_channel_data()[self.channel][0].is_cw():
            self.channelPulse = 'CW'
        elif self.ek_data.get_channel_data()[self.channel][0].is_fm():
            self.channelPulse = 'FM'
        else:
            print('No pulse form identified')
            pass
    
    def writeChannels(self):
        if self.within_channel:
            self.buildRIA()
            for key in self.ria_list:
                self.channelPulse = key
                self.setChannelName()
                out_files = {} # This has to be reset for every channel since the filename key remains the same but the value changes
                out_files[self.nameBase] = os.path.join(self.out_dir, '')+os.path.splitext(self.nameBase)[0]+'_'+self.nameFreq+'_split.raw'
                files_written = self.ek_data.write_raw(out_files, raw_index_array=self.ria_list[key],overwrite=True)
        else:
            self.setChannelName()
            out_files = {} # This has to be reset for every channel since the filename key remains the same but the value changes
            out_files[self.nameBase] = os.path.join(self.out_dir, '')+os.path.splitext(self.nameBase)[0]+'_'+self.nameFreq+'_split.raw'
            files_written = self.ek_data.write_raw(out_files, channel_ids=self.writeChannelIDs,overwrite=True)

    def buildRIA(self):
        self.ria_list = {}
        ria = {}
        for data in self.ek_data.get_channel_data()[self.channel]:
            if data.is_cw():
                ria[data] = np.arange(data.n_pings)
            else:
                ria[data] = np.array([])
        self.ria_list['CW'] = ria
        ria = {}
        for data in self.ek_data.get_channel_data()[self.channel]:
            if data.is_fm():
                ria[data] = np.arange(data.n_pings)
            else:
                ria[data] = np.array([])
        self.ria_list['FM'] = ria

    def main(self):
        print('Writing files to: '+self.out_dir)
        for self.file in (pbar1 := tqdm(self.in_files)):
            pbar1.set_description('Processing '+str(len(self.in_files))+' .raw files')
            self.openFile()
            self.splitFile()
            del self.ek_data # ekdata seems to be storing something that slows things down when it carries it over
        print('Complete')

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_dir", help="Input directory of raw files", type=str)
    parser.add_argument("--in_files", help="define a single or comma-separated list of filenames rather than \
            iterating through all files. If a filename has whitespace you must use double quotes around it", type=str, default=None)
    parser.add_argument("--out_dir", help="Output directory for split files. \
            Default is new directory 'splitFiles' in input directory", type=str, default=None)
    parser.add_argument("--read_prefix", help="File prefix used to identify multichannel/FM file type. \
            E.g., if all target files start in DY_FM..., use --read_prefix=DY_FM", type=str, default='')
    parser.add_argument("--read_suffix", help="File suffix used to identify multichannel/FM file type. \
            E.g., if all target files end in ..._2.raw, use --read_suffix=_2", type=str, default='')
    parser.add_argument("--group", help="If True, create new files containing all channels with same pulse form", type=bool, default=False)
    parser.add_argument("--within_channel", help="If True, split files where multiple pulse forms exist in a single channel", type=bool, default=False)
    parser.add_argument("--overwrite", help="If True, process every file found, regardless of if basename already exists in output", type=bool, default=False)
    parser.add_argument("--ignore_last", help="If True, process every file except the one most recently modified. This is for real-time survey applications", type=bool, default=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parseArguments()
    splitFiles(args)