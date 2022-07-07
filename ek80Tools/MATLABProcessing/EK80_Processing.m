%% Setup
clear all; close all; clc
addpath lib
addpath('lib/bin')

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
    fprintf('Processing %s \n',fbase{i})
    load([outdir '\' fbase{i} '.mat'])
    
    [Nch,Npings] = size(data.echodata);
    fchannels = dir([outdir '\' fbase{1} '*ES*.mat']);
    
    for channel=1:Nch
        %load([outdir '\'  fchannels(channel).name]);
        range{channel,1} = data.echodata(channel,1).range;
        
        %time
        time{channel,1}=[data.echodata(channel,:).timestamp];
        
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
        para{channel,1}.G = getGain(data,channel); % Either grab gain out of gain tables or from the calibration
        
        para{channel,1}.fsdec = 1/data.param(channel, 1).SampleInterval;
        para{channel,1}.fstart = data.param(channel,1).FrequencyStart;
        para{channel,1}.fend = data.param(channel, 1).FrequencyEnd;
        para{channel,1}.fnom = str2num(string(data.config.transceivers(channel).channels.transducer.Frequency));
        para{channel,1}.fc = (para{channel}.fstart + para{channel}.fend) ./2;
        para{channel,1}.alpha = alpha_sea(data.environ.Depth,data.environ.Salinity,...
            data.environ.Temperature,data.environ.Acidity,para{channel}.fc/1000); % in dB m
        
        for ping =1:Npings
            if ~isempty(data.echodata(channel,ping).channelID)
                
                % Now Process Sv
                if strmatch(dlg.Sv,'Yes')
                    [sv, cv] = SvCalc(channel, data, para, range, ping);
                    Sv{channel,1}(:,ping) = sv;
                    CV{channel,1}(:,ping) = cv;
                end
                
                if strmatch(dlg.TS,'Yes')
                    [sp, phialong, phiathw] = TSCalc(channel, data, para, range, ping);
                    Sp{channel}(:,ping) = sp;
                    PhiAlong{channel}(:,ping) = phialong;
                    PhiAthw{channel}(:,ping) = phiathw;
                end
                
            end
        end
        
        if strfind(dlg.Echo,'Yes')
            % [~,fnm,~] =fileparts(rawfiles);
            fout = strcat('Sv_',fbase{i},'_',string(para{channel}.fnom/1000),'kHz');
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
            saveas(h,strcat(outdir,'\Figures\',fout),'png');
            clear fout
            
        end
    end
    
    
    % Save the files
    % If calculating Sv, save separately
    if strfind(dlg.Sv,'Yes')
        fout2 = ['Sv_' fbase{i} '.mat'];
        if ~exist([outdir '\Data\'])
            mkdir([outdir '\Data\'])
        end
        save([outdir '\Data\' fout2],'Sv','CV','time','range','para');
        clear fout2
    end
    
    % If calculating TS, save separately
    if strfind(dlg.TS,'Yes')
        fout2 = ['Sp_' fbase{i} '.mat'];
        if ~exist([outdir '\Data\'])
            mkdir([outdir '\Data\'])
        end
        save([outdir '\Data\' fout2],'Sp','PhiAlong','PhiAthw','time','range','para');
        clear fout2
    end
    
    fclose('all');
    % insert file save here
    clear Sv CV Sp PhiAlong PhiAthw
end