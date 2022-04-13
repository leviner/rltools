%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% help sampledatavec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (This exmaple has 3 split beam transducers and 3 WBTs and 2 pings)
% (All fields have 6 entries in this example: ch1p1 ch2p1 ch3p1 ch1p2 ch2p2 ch3p2)
%
% ALL DATA ARE SAMPLED AT 1.5 MHz
% But there are two stages of decimation.
% Only the first stage of decimation has been applied to the data.
% The second stage allows you to get to the same as the software display.
%
% sampledatamat = 
% 3x2 struct array with fields:
%     channelno %% In this example = [1 2 3 1 2 3]
%     mode_low %%FORMAT OF RAW DATA AS IT IS STORED 
%     mode_high %%FORMAT OF RAW DATA AS IT IS STORED 
%     mode %%FORMAT OF RAW DATA AS IT IS STORED 
%     transducerdepth %% SET IN THE SOFTWARE TRANSDUCER INSTALLATION 
%     frequency %% ACTUAL START FREQUENCIES OF SWEEPS
%     transmitpower %% ACTUAL POWER TRANSMITTED in Watts
%     pulselength %% ACTUAL PulseLength in seconds
%     bandwidth %% ACTUAL BANDWIDTH in Hz
%     sampleinterval %%ACTUAL TIME BETWEEN SAMPLES IN A PING in s
%                   %%( = 1/sampling frequency after stage 1 decimation)
%     soundspeed
%     absorptioncoefficient %% At nominal frequency at 250 m at T and S set
%                               or default values
%     heave
%     roll
%     pitch
%     temperature %% Defaults to 10 degrees and default S is 30? psu  
%     heading
%     transmitmode %% 0=Active  1=Passive  2=Test Signal, Tone, *Input* to test receiver 
%     spare1
%     sweep %% Chirp sweep rate in Hz/s
%     offset 
%     count %% Number of counts in each ping determined from the range
%                   input in the software file save tab and sound speed in 
%                   the software based on temperature and salinity  
%     complexsamples %% The complex data samples
%
%  To access the complex data samples in ch1 ping 2 look at
%  sampledatamat(4).complexsamples
