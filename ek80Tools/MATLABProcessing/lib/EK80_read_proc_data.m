%     EK80_read_proc_data
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
%            version: 'V4'
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

% -----------------------------------------------------------------
% Select files
% ---------------------------------------------------------------

[filename, filepath] = uigetfile({'*.raw'},'Pick a raw data file','MultiSelect','on');

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

% ---------------------------------------------------------
%  Loop over all selected files
% ---------------------------------------------------------

tstart=tic;
for i= 1:length(filename)
    rawfiles = [filepath filename{i}];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Send info to screen for monitoring
    if i>1
        telapsed = toc(tstart);
        percentdone = (i-1)/length(filename)*100;
        remainingtime = round((100/percentdone*telapsed-telapsed)/60);  %remaining minutes
        disp([num2str(round(percentdone)) '% done'])
        disp(['Estimated ' num2str(remainingtime) ' minute(s) left until complete'])
    end
    disp(['(' num2str(i) '/' num2str(length(filename)) ')'])
    
    
    % read data from single file
    data = EK80readRawFiles(rawfiles, Version);     % V3-2019, V4-2020
    
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
    
    %%% COMMENTED OUT THIS SECTION 2/11/2021 by RK
    %     GPStime=[]; Lat=[]; Lon=[];
    %     for k=1:length(data.nmea)
    %
    %
    %         % New for CTRiverDT  AN
    %         %    Just using this NEMA string only, has all the info...
    %
    %         if strcmp(data.nmea(k).text(1:6),'$GPRMC')
    %             [lat,lon,YYYY,MM,DD,GPShh,GPSmm,GPSss] = EK80_parse_nmea_GPRMC(data.nmea(k).text);
    %             GPStime(end+1)=datenum([YYYY' MM' DD' GPShh' GPSmm' GPSss'])';
    %             Lat(end+1)=lat;
    %             Lon(end+1)=lon;
    %         end
    %
    %
    %     end  % end of parsing nmea strings
    
    
    
    
    
    
    clear DD MM YYYY
    Nch=size(data.echodata,1);                % number of transducers
    
    % ---------------------------------------------------------------
    % pulse compression
    %       chu's pulse compression
    % --------------------------------------------------------------
    
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
        %         OldRange=(1:Npings)';   % set up temp space
        %         OldVolt=(1:Npings)';
        %         OldComp=(1:Npings)';
        OldRange=(1:length(data.echodata(channel,1).range))';
        OldVolt=(1:length(data.echodata(channel,1).range))';
        OldComp=(1:length(data.echodata(channel,1).range))'; % changed by RK 7/21/2020
        
        for ping =1:Npings
            if ~isempty(data.echodata(channel,ping).channelID)
                
                datasz=size(data.echodata(channel, ping).range,1);
                NewRange=NaN*ones(max(datasz,size(OldRange,1)),1);
                NewRange(1:size(OldRange,1)) = OldRange;
                NewRange(1:datasz) = data.echodata(channel, ping).range;
                OldRange=NewRange;
                
                
                if isfield(data.echodata,'compressed')
                    datasz=size(data.echodata(channel, ping).compressed,1);
                    NewComp=NaN*ones(max(datasz,size(OldComp,1)),Npings);
                    NewComp(1:size(OldComp,1),1:size(OldComp,2))=OldComp;
                    % mean along each row
                    NewComp(1:datasz,ping) = mean(data.echodata(channel,ping).compressed,2);
                    OldComp=NewComp;
                    CompressedVoltage=NewComp;
                end
                
                %Voltage  -  New code 12/26/2018 ES333 may have single sensor, check
                
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
        
        clear OldRange NewRange OldVolt NewVolt OldComp NewComp
        
        
        % ----------------------------------------------------
        %     save data in mat files
        % ----------------------------------------------------
        ndx1 = strfind(ChannelID,'ES');
        ndx2 = strfind(ChannelID,'-');
        %ndx2 = ndx2(min(find(ndx2>=ndx1)))-1;
        %ID = ChannelID(ndx1:ndx2);
        %ID = ChannelID(ndx1(1):ndx2(2)-1); % RK added and commented out previous 2 lines 7/21/2020
        ID = ChannelID(ndx1(1):end); % RK added and commented out previous 2 lines 7/21/2020
        
        StrName = char(filename{i}); StrName = StrName(1:end-4);
        
        if ~exist('GPStime','var') % Added by RK 6/2021
            GPStime = [];
        elseif ~exist('Lat','var')
            Lat = [];
        elseif ~exist('Lon','var')
            Lon = [];
        end
        
        if exist('CompressedVoltage','var')
            Npings=size(CompressedVoltage,2);
            disp(['Saving ',StrName '_',ID,'.mat'])
            % filepath1='./CTRiverDT_matfiles/';
            
            
            filepath1 = MatOutDir;      % set in init file,
            %             save(fullfile(filepath1,[StrName,'_',ID,'.mat']),'ChannelID','Npings',...
            %                 'GPStime','Lat','Lon','ComputerTime','FreqStart','FreqEnd',...
            %                 'Range','Voltage','CompressedVoltage');
            save(fullfile(filepath1,[StrName,'_',ID,'.mat']),'ChannelID','Npings',...
                'ComputerTime','FreqStart','FreqEnd',...
                'Range','Voltage','CompressedVoltage');
        else
            Npings=size(Voltage,2);
            disp(['Saving ',StrName '_',ID,'.mat'])
            % filepath1='./CTRiverDT_matfiles/';
            
            
            
            filepath1 = MatOutDir;       % set in init file
            %            save(fullfile(filepath1,[StrName,'_',ID,'.mat']),'ChannelID','Npings',...
            %                'GPStime','Lat','Lon','ComputerTime','FreqStart','FreqEnd',...
            %                'Range','Voltage');
            save(fullfile(filepath1,[StrName,'_',ID,'.mat']),'ChannelID','Npings',...
                'ComputerTime','FreqStart','FreqEnd',...
                'Range','Voltage');
        end
        clear Range Voltage CompressedVoltage
    end
    
    if strfind(dlg.Cal,'Yes')
        data = addCals(data, cals);
    end
    
    disp('Saving Data Struct')
    save(fullfile(filepath1,[filename{i}(1:end-4),'.mat']),'data');
    clear data;
    fclose('all');
end  % end loop over files


% wake up after long run with multiple files
yyy=load('chirp.mat');
sound(yyy.y*.5);


telapsed = toc(tstart);
fprintf('\n Done!   This run took %.2f minutes.... \n',telapsed/60);


% -----------------------------------------------------------------------
