from echolab2.instruments import echosounder
from glob import glob
import os
import argparse

class splitFiles():

    def __init__(self, args):
        if args.out_dir is None:
           self.out_dir = args.in_dir+'/splitFiles/'
        else:
            self.out_dir = args.out_dir
        self.in_dir = args.in_dir
        self.in_files = glob(self.in_dir+'/*'+args.suffix+'.raw')
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.main()
        
    def openFile(self):
            self.nameBase = os.path.basename(self.file)
            self.ek_data = echosounder.read([self.file])[0]
    
    def splitFile(self):
        for self.channel in self.ek_data.channel_ids:            
            self.getchannelFreq()
            self.writeChannel()        

    def getchannelFreq(self):
        if self.ek_data.get_channel_data()[self.channel][0].is_cw():
            self.nameFreq = str(int(self.ek_data.get_channel_data()[self.channel][0].frequency[0]/1000))+'CW'
        elif self.ek_data.get_channel_data()[self.channel][0].is_fm():
            self.nameFreq = str(int(self.ek_data.get_channel_data()[self.channel][0].frequency_start[0]/1000))+'-'+\
                        str(int(self.ek_data.get_channel_data()[self.channel][0].frequency_end[0]/1000))+'FM'
        else:
            print('No pulse form identified')
            pass
    
    def writeChannel(self):
        out_files = {}
        out_files[self.nameBase] = os.path.join(self.out_dir, '')+os.path.splitext(self.nameBase)[0]+'_'+self.nameFreq+'_split.raw'
        files_written = self.ek_data.write_raw(out_files, channel_ids=[self.channel],overwrite=True)

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
    parser.add_argument("--out_dir", help="Output directory for split files. \
            Default is new directory 'splitFiles' in input directory", type=str, default=None)
    parser.add_argument("--suffix", help="File suffix used to identify multichannel/FM file type. \
            E.g., if all target files end in ...-2.raw, use --suffix=-2", type=str, default='')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # Parse the arguments
    args = parseArguments()
    splitFiles(args)