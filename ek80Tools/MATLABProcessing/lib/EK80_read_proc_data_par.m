%     EK80_read_proc_data_par
%
%    parallelized version
%
% Main routine that provides examples of how to read and process
% EK80 broadband data based on the codes provided by
% Andone Lavery, WHOI and
% Gavin Macaulay, IMR, Norway on 11-02-2014 (gavin.macaulay@imr.no)
%
%   Saves the data in .mat files for imaging and analysis
%
%   DATA FORMAT (subject to change):
%           echodata: [  struct]
%            filters: [  cFilter1Data]
%             config: [  struct]
%            environ: [  struct]
%               nmea: [  cNMEA]
%              param: [  struct]
%         channelMap: [  containers.Map]
%                mru: [  struct]
%         channelIDs: [  containers.Map]
%            version: 'V3'
%         parameters: [1x1 struct]
%
%
%   Run the init routine first to initialize data and output space
%
%
% Modified by Dezhang Chu, NOAA Fisheries, NWFSC (dezhang.chu@noaa.gov )
% March 7, 2016
% Modified by Matt perkins, URI, March 2017

% routine taken from EK80_read_files_Perkins.m
%   renamed, modified
%   Newhall, 2017
%
% Making more generic
%   Newhall 12/2018
%
% revised to work with empty data files
%   Newhall 2019
% -----------------------------------------------------------------

close all             % start fresh

% if these files are all ready added, then will give a warning
%

Pool = gcp('nocreate');
warning('off');
addAttachedFiles(Pool, 'EK80readRawFiles');
addAttachedFiles(Pool, 'EK80pulseCompress');
addAttachedFiles(Pool, 'EK80readRawV4');
addAttachedFiles(Pool, 'EK80calcSentSignal');
addAttachedFiles(Pool, 'parseconfxmlstruct');
addAttachedFiles(Pool, 'EK80_parse_nmea_GPRMC');
addAttachedFiles(Pool, 'EK80chirp');
addAttachedFiles(Pool, 'Save_vars');
addAttachedFiles(Pool, 'Save_data');


% -----------------------------------------------------------------
% Select files
% -----------------------------------------------------------------

[filename, filepath] = uigetfile(append(string(DataFilePath),'*.raw'),'Pick a raw data file','MultiSelect','on');

dlg.Cal = questdlg('Do you have xml calibration files?','Calibration','Yes','No', 'Yes');
if strfind(dlg.Cal,'Yes')
    [calfilename, calfilepath] = uigetfile(append(DataFilePath,'*.xml'),'select xml calibration files','MultiSelect','on');
    cals = parseCals(calfilename, calfilepath);
end

% GUI check...  comment out if hardwiring data files.
% if not a cell, turn it into one, if choose 1 file, not a cell
if ~iscell(filename)
    filename = {filename};
end

% returns 0 if cancel hit
if filename{1} == 0
    fprintf('\tError: File(s) not found...\n');
    return
end

% skip files that have already been processed (RK 2021)
match = ones(1,length(filename));
processedfiles = dir(MatOutDir');
filesdone = {processedfiles.name};
for f = 1:length(filename)
    s = filename(f);s = s{1};s = split(s,'.');s = s(1); s = s{1};
    matches = startsWith(string(filesdone),s);
    if isempty(find(matches == 1))
        match(f) = 0;
    end
end
filename = filename(find(match == 0));
if isempty(filename)
    fprintf('\t All files have already been processed. \n');
    return
end


%
%  Loop over all selected files
% ---------------------------------------------------------

tstart=tic;


%cc = 1;
%parfor i = 1:length(filename)
for i= 1:length(filename)
    rawfiles = [filepath filename{i}];
    StrName = char(filename{i});
    
    fprintf(['\t Processing file ' num2str(i) ' of ' num2str(length(filename)) '. \n']);
    %cc = cc + 1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % read data from single file
    data = EK80readRawFiles(rawfiles,'V4');
    
    % fixing partial files 2019
    if isempty(data)
        fprintf('\t Skipping file: %s\n',rawfiles);
        continue
    end
    
    
    %           DATA:  (see README for all definitions)
    %       echodata: [  struct]
    %        filters: [  cFilter1Data]
    %         config: [  struct]
    %        environ: [  struct]
    %           nmea: [  cNMEA]
    %          param: [  struct]
    %     channelMap: [  containers.Map]
    %            mru: [  struct]
    %     channelIDs: [  containers.Map]
    %        version: 'V3'
    %     parameters: [  struct]
    
    
    
    
    
    % ------------------------------------------------------------
    % make GPS vectors from NMEA strings
    %
    %   Lat, Lon, GPStime
    %
    %   Looks for nmea string GPRMC
    %
    %   if not complete, use NaNs
    % ------------------------------------------------------------
    
    GPStime=[]; Lat=[]; Lon=[];
    for k=1:length(data.nmea)
        
        
        % New for CTRiverDT  AN
        %    Just using this NEMA string only, has all the info...
        
        if strcmp(data.nmea(k).text(1:6),'$GPRMC')
            [lat,lon,YYYY,MM,DD,GPShh,GPSmm,GPSss] = EK80_parse_nmea_GPRMC(data.nmea(k).text);
            GPStime(end+1)=datenum([YYYY' MM' DD' GPShh' GPSmm' GPSss'])';
            Lat(end+1)=lat;
            Lon(end+1)=lon;
        end
        
        
    end  % end of parsing nmea strings
    
    
    Nch=size(data.echodata,1);                % number of transducers
    
    % ---------------------------------------------------------------
    % pulse compression
    %       chu's pulse compression
    % --------------------------------------------------------------
    %
    for channel=1:Nch
        
        data = EK80pulseCompress(data, channel);
        
        close all
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract parameters from data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Name of Channel
        ChannelID = data.channelIDs(channel);
        %number of pings per file
        Npings = size(data.echodata(channel,:),2);
        %Computer time
        ComputerTime=[data.echodata(channel,:).timestamp];
        %FrequencyStart
        FreqStart=data.param(channel,1).FrequencyStart;
        %FrequencyEnd
        FreqEnd=data.param(channel,1).FrequencyEnd;
        %Range
        %OldRange=(1:Npings)';
        %OldVolt=(1:Npings)';
        %OldComp=(1:Npings)';
        OldRange=(1:length(data.echodata(channel,1).range))';
        OldVolt=(1:length(data.echodata(channel,1).range))';
        OldComp=(1:length(data.echodata(channel,1).range))'; % changed by RK 7/21/2020
        
        for ping =1:Npings
            if ~isempty(data.echodata(channel,ping).channelID)
                
                datasz=size(data.echodata(channel, ping).range,1);
                NewRange=NaN*ones(max(datasz,size(OldRange,1)),1);
                NewRange(1:size(OldRange,1))=OldRange;
                NewRange(1:datasz) = data.echodata(channel, ping).range;
                OldRange=NewRange;
                
                
                if isfield(data.echodata,'compressed')
                    datasz=size(data.echodata(channel, ping).compressed,1);
                    NewComp=NaN*ones(max(datasz,size(OldComp,1)),Npings);
                    %NewComp = NaN*ones(datasz,Npings); % Changed by RK 7/21/2020
                    NewComp(1:size(OldComp,1),1:size(OldComp,2))=OldComp;
                    NewComp(1:datasz,ping) = mean(data.echodata(channel,ping).compressed,2);
                    OldComp=NewComp;
                    CompressedVoltage=NewComp;
                    
                end
                
                %Voltage  -  New code 12/26/2018 ES333 may have single sensor
                
                if isvector(data.echodata(channel, ping).complexsamples)
                    datasz=length(data.echodata(channel, ping).complexsamples);
                    NewVolt=NaN*ones(max(datasz,size(OldVolt,1)),Npings);
                    NewVolt(1:size(OldVolt,1),1:size(OldVolt,2))=OldVolt;
                    NewVolt(1:datasz,ping) = data.echodata(channel,ping).complexsamples;
                    OldVolt=NewVolt;
                else
                    datasz=size(data.echodata(channel, ping).complexsamples,1);
                    NewVolt=NaN*ones(max(datasz,size(OldVolt,1)),Npings);
                    NewVolt(1:size(OldVolt,1),1:size(OldVolt,2))=OldVolt;
                    NewVolt(1:datasz,ping) = mean(data.echodata(channel,ping).complexsamples,2);
                    OldVolt=NewVolt;
                end
                
            end
        end
        Range=NewRange;
        Voltage=NewVolt;
        
        
        
        
        % ----------------------------------------------------
        %     save data in mat files
        % ----------------------------------------------------
        
        ndx1 = strfind(ChannelID,'ES');
        %ndx2 = strfind(ChannelID,'-');
        %ndx2 = ndx2(min(find(ndx2>=ndx1)))-1;
        %ID = ChannelID(ndx1:ndx2);
        %ID = ChannelID(ndx1(1):ndx2(2)-1); % RK added and commented out previous 2 lines 7/21/2020
        ID = ChannelID(ndx1:end-2); % edited by Bassett on 2/14/22
        StrName = char(filename{i}); StrName = StrName(1:end-4);
        
        % COMMENT OUT THESE LINES IF RUNNING PARALLEL
        %          if ~exist('GPStime','var') % Added by RK 6/2021
        %              GPStime = [];
        %          elseif ~exist('Lat','var')
        %              Lat = [];
        %          elseif ~exist('Lon','var')
        %              Lon = [];
        %          end
        
        
        % New 2019  AEN
        %  Checking for Compressed Voltage and writing accordingly
        
        if isfield(data.echodata,'compressed')
            
            
            Npings=size(CompressedVoltage,2);
            
            S = [MatOutDir StrName '_',ID,'.mat' ];
            disp(['Saving ',S])
            
            Save_vars(S,ChannelID,Npings,...
                GPStime,Lat,Lon,ComputerTime,FreqStart,FreqEnd,...
                Range,Voltage,CompressedVoltage);
            
        else
            
            Npings=size(Voltage,2);
            CompressedVoltage = NaN*length(Voltage);
            
            S = [MatOutDir StrName '_',ID,'.mat' ];
            disp(['Saving ',S])
            
            Save_vars(S,ChannelID,Npings,...
                GPStime,Lat,Lon,ComputerTime,FreqStart,FreqEnd,...
                Range,Voltage,CompressedVoltage);
            
        end
        
    end
    if strfind(dlg.Cal,'Yes')
        data = addCals(data, cals);
    end
    
    S = [MatOutDir StrName '.mat' ];
    Save_data(S,data);
    
    fclose('all');
    disp(['     Done with ',StrName]);
    
    
    
end  % end loop over files


% wake up after long run with multiple files
yyy=load('chirp.mat');
sound(yyy.y*.5);

telapsed = toc(tstart);
fprintf('\n Done!   This run took %.2f minutes.... \n',telapsed/60);
% -----------------------------------------------------------------------
