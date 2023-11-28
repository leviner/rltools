% Ek80 initialization and configuration file
%
%     Sets up input and output data folders
%
%   Currently, user is prompted for all directories. 
%   This can be hardocded if desired for use or testing.
%
%
%  DO NOT CHANGE VARIABLE NAMES
%
%    Set Version
%    Set Num_parallel_workers
%    
% -------------- Set version and workers -------------------- %

clear; close all

% If needed, modify the version number for older EK80 raw formats
Version = 'V4';    % 2020+, 'V3' can be used for 2016-2019 EK80 filetypes

% If running in parallel, default is 4 workers
Num_parallel_workers = 0;% 4; % 32 ; % 0 can be used to specify no parallel

%----------------- Setup folders -------------------------%
WorkDir  = [pwd '\'];   % MUST end with slash 
addpath([WorkDir 'lib\'])                  % utility programs folder, end with slash

% Path to input .raw files for processing, MUST end with slash
%    if using GUI, will use as repository, not input files
DataFilePath = [uigetdir(pwd,'Select Input Data Directory') '\'];

% path to folder that will contain .mat files that hold processed data
MatOutDir = [uigetdir(DataFilePath,'Select Output Directory') '\']

% path to folder that will hold images and plots. End with a slash
PlotOutDir =  [MatOutDir 'proc_plots\'];

% check if data dir exists
if ~exist(DataFilePath,'dir')
       error(['Cannot find ' DataFilePath]);
end
 
% check if the others exist
if ~exist(MatOutDir,'dir')
       mkdir(MatOutDir);
end

if ~exist(PlotOutDir,'dir')
       mkdir(PlotOutDir);
end
% ----------------------------------------------------------------------- 