%% Setup
clear all; close all; clc
addpath lib

dlgTitle    = 'Outputs';
dlg.Sv = questdlg('Output processed volume scattering?',dlgTitle,'Yes','No', 'Yes');
dlg.TS = questdlg('Output processed TS?',dlgTitle,'Yes','No', 'Yes');
if strfind(dlg.Sv,'Yes');
    dlg.Echo = questdlg('Output Sv Echograms?',dlgTitle,'Yes','No', 'Yes');
end

% prompt user for directory for results if not in workspace
if ~exist('MatOutDir')
    outdir = uigetdir(pwd,'Select Directory for Results');
else
    outdir = MatOutDir;
end

dlg.Cal = questdlg('Do you have calibration files?',dlgTitle,'Yes','No', 'Yes');
if strfind(dlg.Cal,'Yes');
    [filename, filepath] = uigetfile(append(outdir,'*.xml'),'select xml calibration files','MultiSelect','on');
else
    cal = 0; % This is a placeholde RML
end

files = dir([outdir '\*.mat']);
ct=1;
for jj = 1:length(files)
    if isempty(strfind(files(jj).name,'_ES'))
        nm = strsplit(files(jj).name,'.mat');
        fbase{ct} = nm{1};
        ct=ct+1;
    end
end
%%
for i = 1:length(fbase)
    load([outdir '\' fbase{i} '.mat'])
    
    [Nch,Npings] = size(data.echodata);
    fchannels = dir([outdir '\' fbase{1} '*ES*.mat']);
    
    % Set environment from header
    w.z = data.environ.Depth; % Hold depth
    w.S = data.environ.Salinity; % Hold salinity
    w.T = data.environ.Temperature;  % Hold Temp
    w.pH = data.environ.Acidity;  % Hold pH
    dens = gsw_rho(w.S,w.T,w.z); % sea water density at estimate of density
    % need to fix all of this for deep water application
    w.P = dens*9.81*w.z*1e-4; % pressure in decibars
    w.c = gsw_sound_speed(w.S,w.T,w.P);
    clear dens
        
    for channel=1:Nch
        %load([outdir '\'  fchannels(channel).name]);
        range{channel,1} = data.echodata(channel,1).range;
        
        %time
        time{channel,1}=[data.echodata(channel,:).timestamp];
        
        if ~cal
            para{channel,1}.taueff = data.param(channel, 1).tauEffective;
            para{channel,1}.ptx = data.param(channel, 1).TransmitPower;
            para{channel,1}.zet = 75;  % transducer impedance # impedance from the wbat and carry over from the datagram
            % There are some weird cases in the xml reading where things are
            % sometimes characters and sometimes doubles so I just make
            % everything strings RML
            para{channel,1}.zer = str2num(string(data.config.transceivers(channel).Impedance)); % transceiver imedpance
            para{channel,1}.psinom = str2num(string(data.config.transceivers(channel).channels.transducer.EquivalentBeamAngle));
            para{channel,1}.sensalong = str2num(string(data.config.transceivers(channel).channels.transducer.AngleSensitivityAlongship));
            para{channel,1}.sensathw = str2num(string(data.config.transceivers(channel).channels.transducer.AngleSensitivityAthwartship));
            para{channel,1}.G = getGain(data,channel); % This became complicated because of gain tables
            
            para{channel,1}.fsdec = 1/data.param(channel, 1).SampleInterval;
            para{channel,1}.fstart = data.param(channel,1).FrequencyStart;
            para{channel,1}.fend = data.param(channel, 1).FrequencyEnd;
            para{channel,1}.fnom = str2num(string(data.config.transceivers(channel).channels.transducer.Frequency));
            para{channel,1}.fc = (para{channel}.fstart + para{channel}.fend) ./2;
            para{channel,1}.alpha = alpha_sea(w.z,w.S,w.T,w.pH,para{channel}.fc/1000); % in dB m
        end
        
        for ping =1:Npings
            if ~isempty(data.echodata(channel,ping).channelID)
                
                % Now Process Sv
                if strmatch(dlg.Sv,'Yes')
                    [sv, cv] = SvCalc(channel, data, para, w, range, ping);
                    Sv{channel,1}(:,ping) = sv;
                    CV{channel,1}(:,ping) = cv;
                end
                
                if strmatch(dlg.TS,'Yes')
                    [sp, phialong, phiathw] = TSCalc(channel, data, para, w, range, ping);
                    Sp{channel}(:,ping) = sp;
                    PhiAlong{channel}(:,ping) = phialong;
                    PhiAthw{channel}(:,ping) = phiathw;
                end
                
            end
        end
        
        if strfind(dlg.Echo,'Yes')
            % [~,fnm,~] =fileparts(rawfiles);
            fout = ['Sv_' fbase{i}];
            h = figure('visible','off','units','inch','position',[2,2,6,3]);
            imagesc(1:Npings,range{channel}, Sv{channel})
            xlabel('Ping','fontweight','bold')
            ylabel('Range [m]','fontweight','bold')
            set(gca,'clim',[-70 -34])
            set(gca,'ylim',[0 25*floor(max(range{channel})/25)])
            hcb=colorbar;
            ylabel(hcb,'Sv [dB re 1/m]','fontweight','bold')
            set(gca,'linewidth',2)
            cmap = crameri('oslo',50);
            colormap(cmap);
            if ~exist([outdir '\Figures\'])
                mkdir([outdir '\Figures\'])
            end
            saveas(h,[outdir '\Figures\' fout],'png');
            clear fout
            
        end
    end
    
    
    % Save the files
    % If calculating Sv, save separately
    if strfind(dlg.Sv,'Yes')
        % [~,fnm,~] =fileparts(rawfiles);
        fout2 = ['Sv_' fbase{i} '.mat'];
        % save([outdir '\Data\' fout2],'Sv','ComplexVoltage','pings','time','range','para','w','gps');
        if ~exist([outdir '\Data\'])
            mkdir([outdir '\Data\'])
        end
        save([outdir '\Data\' fout2],'Sv','CV','time','range','para','w');
        % need to edit to bring in ship GPS data
        clear fout2
    end
    
    % If calculating TS, save separately
    if strfind(dlg.TS,'Yes')
        % [~,fnm,~] =fileparts(rawfiles);
        fout2 = ['Sp_' fbase{i} '.mat'];
        % save([outdir '\Data\' fout2],'Sv','ComplexVoltage','pings','time','range','para','w','gps');
        if ~exist([outdir '\Data\'])
            mkdir([outdir '\Data\'])
        end
        save([outdir '\Data\' fout2],'Sp','PhiAlong','PhiAthw','time','range','para','w');
        % need to edit to bring in ship GPS data
        clear fout2
    end
    
    fclose('all');
    % insert file save here
    clear Sv CV Sp PhiAlong PhiAthw
end