# coding=utf-8

from echolab2.instruments import echosounder
from glob import glob
import os
import argparse
from tqdm import tqdm

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

    usage: splitPEL.py [-h] [--out_dir OUT_DIR] [--suffix SUFFIX] in_dir

    positional arguments:
    in_dir             Input directory of raw files

    optional arguments:
    -h, --help           show this help message and exit
    --in_files IN_FILES  define a single or comma-separated list of filenames rather than iterating through all files. 
                            If a filename has whitespace you must use double quotes around it
    --out_dir OUT_DIR    Output directory for split files. Default is new directory 'splitFiles' in input directory
    --suffix SUFFIX      File suffix used to identify multichannel/FM file type. E.g., if all target files end in 
                            ..._2.raw, use --suffix=_2
    --prefix PREFIX     File prefix used to identify multichannel/FM file type. E.g., if all target files start in 
                            DY_FM..., use --prefix=DY_FM
    --group GROUP       If True, create new files containing all channels with same pulse form
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
            self.in_files = glob(self.in_dir+'/'+args.prefix+'*'+args.suffix+'.raw')
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.group = args.group
        if self.group is True:
            self.channelGroups = {'CW':[],'FM':[]}
        self.main()
        
    def openFile(self):
            self.nameBase = os.path.basename(self.file)
            self.ek_data = echosounder.read([self.file])[0]
    
    def splitFile(self):
        for self.channel in self.ek_data.channel_ids:
            self.getchannelFreq()
            if self.group is True:
                self.channelGroups[self.channelPulse].append(self.channel)
            else:
                self.setChannelName()
                self.writeChannelIDs = [self.channel]
                self.writeChannels() 
        if self.group is True:
            for k in self.channelGroups.keys():
                if self.channelGroups[k]:
                    self.writeChannelIDs = self.channelGroups[k]
                    self.nameFreq = k
                    self.writeChannels()

    def setChannelName(self):
        if self.channelPulse == 'CW':
            self.nameFreq = str(int(self.ek_data.get_channel_data()[self.channel][0].frequency[0]/1000))+'CW'
        elif self.channelPulse == 'FM':
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
        out_files = {} # This has to be reset for every channel since the filename key remains the same but the value changes
        out_files[self.nameBase] = os.path.join(self.out_dir, '')+os.path.splitext(self.nameBase)[0]+'_'+self.nameFreq+'_split.raw'
        files_written = self.ek_data.write_raw(out_files, channel_ids=self.writeChannelIDs,overwrite=True)

    def main(self):
        print('Processing '+str(len(self.in_files))+' .raw files\nWriting files to: '+self.out_dir)
        for self.file in self.in_files:
            self.openFile()
            print('Splitting '+self.nameBase+': '+str(len(self.ek_data.channel_ids))+' channels')
            self.splitFile()
        print('Complete')

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_dir", help="Input directory of raw files", type=str)
    parser.add_argument("--in_files", help="define a single or comma-separated list of filenames rather than \
            iterating through all files. If a filename has whitespace you must use double quotes around it", type=str, default=None)
    parser.add_argument("--out_dir", help="Output directory for split files. \
            Default is new directory 'splitFiles' in input directory", type=str, default=None)
    parser.add_argument("--prefix", help="File prefix used to identify multichannel/FM file type. \
            E.g., if all target files start in DY_FM..., use --prefix=DY_FM", type=str, default='')
    parser.add_argument("--suffix", help="File suffix used to identify multichannel/FM file type. \
            E.g., if all target files end in ..._2.raw, use --suffix=_2", type=str, default='')
    parser.add_argument("--group", help="If True, create new files containing all channels with same pulse form", type=bool, default=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parseArguments()
    splitFiles(args)