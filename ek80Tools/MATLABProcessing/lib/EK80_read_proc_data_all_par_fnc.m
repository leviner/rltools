function EK80_read_proc_CTRiver_r243_data_all_v2_par_fnc(DataFilePath, all_files, MatOutDir)
%  function data = EK80_read_proc_CTRiver_r243_data_all_v2_par_fnc(DataFilePath, all_files, MatOutDir)
%
%    all_files must be strcucture from dir command
%       ie   all_files = dir(*.raw);
%
%       if no argument, looks for all files
%
% This routine loops through ALL raw data files, 
%       DOES NOT USE GUI. Use EK80_read_proc_CTRiver_r243_data
%       to bring up GUI.
%
%   Initializes with configuration file.
%
% Main routine that provides examples of how to read and process
% EK80 broadband data based on the codes provided by
% Andone Lavery, WHOI and
% Gavin Macaulay, IMR, Norway on 11-02-2014 (gavin.macaulay@imr.no)
%
%   Saves the data in .mat files for imaging and analysis
%
%   DATA FORMAT :
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
% Modified by Dezhang Chu, NOAA Fisheries, NWFSC (dezhang.chu@noaa.gov )
% March 7, 2016
% Modified by Matt perkins, URI, March 2017
%
% routine taken from EK80_read_files_Perkins.m
%   renamed, modified, reconfigured  
%
%   This software takes into account multiple channels
%
%  1 transducer, multiple channels
%       200kHz  3 quadrant (3 channels)
%       333kHz  single beam (1 channel)
%
%
%
%   Newhall, 2017
% -----------------------------------------------------------------

    

% -----------------------------------------------------------------
% doing ALL files in folder
% ---------------------------------------------------------------

% just one for testing
% all_files = dir([DataFilePath 'R243-D20170626-T201225.raw']); 

if nargin == 0
    error('Please include filenames in an argument\n');   
end

% ---------------------------------------------------------
%  Loop over all selected files
% ---------------------------------------------------------



% if these files are all ready added, then will give a warning
%

Pool = gcp('nocreate'); 
warning('off');
addAttachedFiles(Pool, 'EK80readRawFiles_DT');
addAttachedFiles(Pool, 'EK80pulseCompressr_DT');
addAttachedFiles(Pool, 'EK80readRawV3_DT');
addAttachedFiles(Pool, 'EK80calcSentSignal_DT');
addAttachedFiles(Pool, 'parseconfxmlstruct');
addAttachedFiles(Pool, 'EK80_parse_nmea_GPRMC');
addAttachedFiles(Pool, 'EK80chirp');
addAttachedFiles(Pool, 'Save_vars');
addAttachedFiles(Pool, 'Save_data'); 



tstart=tic;

parfor i=1:length(all_files)    % skip first file, partial and no data
 
    filename = all_files(i).name;
    rawfiles = [DataFilePath filename];        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % not needed for parallel
% %     if i>1
% %         telapsed = toc(tstart);
% %         percentdone = (i-1)/length(all_files)*100;
% %         remainingtime = round((100/percentdone*telapsed-telapsed)/60);  %remaining minutes
% %         disp([num2str(round(percentdone)) '% done'])
% %         disp(['Estimated ' num2str(remainingtime) ' minute(s) left until complete'])
% %     end
% %     disp(['(' num2str(i) '/' num2str(length(all_files)) ')'])
    
    
    data = EK80readRawFiles_DT(rawfiles,'V3');
    
    
    %           DATA:
    %       echodata: [5x7 struct]
    %        filters: [5x2 cFilter1Data]
    %         config: [1x1 struct]
    %        environ: [1x1 struct]
    %           nmea: [1x58 cNMEA]
    %          param: [5x7 struct]
    %     channelMap: [5x1 containers.Map]
    %            mru: [1x7 struct]
    %     channelIDs: [5x1 containers.Map]
    %        version: 'V3'
    %     parameters: [1x1 struct]
    
    
    
    
   
    % ------------------------------------------------------------
    % make GPS vectors from NMEA strings
    %
    %   Lat, Lon, GPStime
    %   if not complete, will be empty
    %
    %     Using ONLY nmea string $GPRMC which contains
    %           all necessary information.    Newhall
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
    
    
    % %  clear DD MM YYYY             % no clears in parallel
    Ntr=size(data.echodata,1);                % number of transducers
    
    % ---------------------------------------------------------------
    % pulse compression
    %       chu's pulse compression
    %
    %    rewritten to handle multiple channels - Newhall
    % --------------------------------------------------------------
    
    for trans=1:Ntr
        
        data = EK80pulseCompressr_DT(data, trans);
        % close all
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract parameters from data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Name of Channel
        ChannelID = data.channelIDs(trans);
        % number of pings per file
        Npings = size(data.echodata(trans,:),2);
        % Computer time
        ComputerTime=[data.echodata(trans,:).timestamp];
        % FrequencyStart
        FreqStart=data.param(trans,1).FrequencyStart;
        % FrequencyEnd
        FreqEnd=data.param(trans,1).FrequencyEnd;
        % Range
        OldRange=(1:Npings)';
        OldVolt=(1:Npings)';
        OldComp=(1:Npings)';
        
        for ping =1:Npings
            if ~isempty(data.echodata(trans,ping).channelID)
                    datasz=size(data.echodata(trans, ping).range,1);
                    NewRange=NaN*ones(max(datasz,size(OldRange,1)),1);
                    NewRange(1:size(OldRange,1))=OldRange;
                    NewRange(1:datasz) = data.echodata(trans, ping).range;
                    OldRange=NewRange;


                if isfield(data.echodata,'compressed')
                    datasz=size(data.echodata(trans, ping).compressed,1);
                    NewComp=NaN*ones(max(datasz,size(OldComp,1)),Npings);
                    NewComp(1:size(OldComp,1),1:size(OldComp,2))=OldComp;
                    NewComp(1:datasz,ping) = mean(data.echodata(trans,ping).compressed,2);
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
                        
        % clear OldRange NewRange OldVolt NewVolt OldComp NewComp
        

        % ----------------------------------------------------
        %     save data in mat files
        % ----------------------------------------------------
            ndx1 = strfind(ChannelID,'ES');
            ndx2 = strfind(ChannelID,'-'); ndx2 = ndx2(min(find(ndx2>=ndx1)))-1;
            ID = ChannelID(ndx1:ndx2);
     
% % % %         if exist('CompressedVoltage','var')
            Npings=size(CompressedVoltage,2);
            % % disp(['Saving ',ChannelID(15:end)])
            disp(['Saving ',[filename(1:end-4) '_',ID,'.mat']])
 
            S = [MatOutDir [filename(1:end-4) '_',ID,'.mat'] ];
            Save_vars(S,ChannelID,Npings,...
                 GPStime,Lat,Lon,ComputerTime,FreqStart,FreqEnd,...
                 Range,Voltage,CompressedVoltage);
           
            
% % % %         else
% % % %             
% % % %             Npings=size(Voltage,2);
% % % %             disp(['Saving ',ChannelID(15:end)])
% % % %           
% % % %             save(fullfile(MatOutDir,[filename(1:end-4),'_',ChannelID(15:end),'.mat']),'ChannelID','Npings',...
% % % %                 'GPStime','Lat','Lon','ComputerTime','FreqStart','FreqEnd',...
% % % %                 'Range','Voltage');
% % % %         end
        
        % clear Range Voltage CompressedVoltage
    end

 
S = [MatOutDir filename(1:end-4) '_data.mat' ];
Save_data(S,data);

fclose('all');
fprintf('\t Done with %s \n',filename);

end  % end loop over files

yyy=load('gong.mat');
sound(yyy.y*.5);

telapsed = toc(tstart);
fprintf('\n Done!   This run took %.2f minutes.... \n',telapsed/60);


