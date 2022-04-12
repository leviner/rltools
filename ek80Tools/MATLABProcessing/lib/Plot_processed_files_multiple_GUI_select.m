% function  plot_processed_files_multiple_GUI_select
%
%   plot processed files in matlab output directory
%       mean plus all the quandrants
%
%   Will append files together for plotting
%
% NOTE: MUST SELECT FILES WITH SAME FREQUENCIES
%
%
%   files hardwired at the moment. Run as script.
%    
%  Quick and dirty....  for first looks
%
%   Have to run init program first to set working space
%
%    User options/flags can be set first
%   at bridge at 10:45:53
%
%     axes for James River
%       start: 29-Nov-2018 10:21:30
%       end:   29-Nov-2018 10:48:48 
% 
%% -----------------------------------------------------------
 
close all;  

% User options ----------------------------------------------
 
save_plots = 1;     % save the plots,  or not...
     plotFilename = 'Nov29_1000L_es333.png';  % edit dir scan below
     % plotFilename = 'Nov29_1000L_es200.png';
     % plotFilename = 'Nov29_1000L_es70.png';
     
GUIselect = 0;      %  This flag (1) opens GUI 
Intrp_shading = 0;  %  0: shading flat, 1: shading interp
    
%% -----------------------------------------------------------------
% Select files 
% ---------------------------------------------------------------
clear Files;

if GUIselect
    [filename, filepath] = uigetfile({'*.mat'},'Pick matlab ES files to append and plot','MultiSelect','on');
    % % % % GUI check...  comment out if hardwiring data files.
    % if not a cell, turn it into one, if choose 1 file, not a cell
    if ~iscell(filename)
        filename = {filename};
    end
    
    % returns 0 if cancel hit
    if filename{1} == 0
        fprintf('\t No files chosen.  Exiting...\n');
        return
    end
    
    % make structs with filename and path
    fct = 0;
    for i=1:length(filename)
        
        % make sure it is a file that contains transducer data 'ES'
        if contains(filename{i},'ES')
            fct = fct + 1;
            Files(fct) = dir([filepath filename{i}]);
        end
        
    end
    
else
    
    Files = dir([MatOutDir '*T15*ES333*.mat']);   
    % Files = dir([MatOutDir '*T15*ES200*.mat']);
    % Files = dir([MatOutDir '*T15*ES70*.mat']);
    
end


tic;

plots_up = 1;    % for making long plot window, must keep if using if below

if(~exist('MatOutDir'))
    fprintf('\t ERROR: Initialize by running init file\n');
    return
end


%% --------------------------------------------------------------------------
Compressed = [];
Time = [];
Ranges = [];
     

first = 1;
proc_files = 1:length(Files);
% proc_files = 1:23;


for ifile = proc_files       
    
    filename = Files(ifile).name;
    load([MatOutDir '/' filename]);
    fprintf('Loading %d of %d - %s\n',ifile, length(proc_files), filename);

    
    % load corresponding data file
  
    ndx1 = strfind(filename,'ES');
    dfilename = [filename(1:ndx1-2) '.mat'];
    load([MatOutDir '/' dfilename]);
    
    % get ID freq from filename 
    ndx2 = strfind(filename,'.');
    ID = filename(ndx1:ndx2-1); 
    
    % find freq ID in data struct, not in same order
    for i = 1:length(data.channelIDs)
        if contains(data.channelIDs(i),ID)
            ID_ndx = i;
            break;
        end
    end
    
     % James River is 5 hours diff, not 4 like Conn River
     ComputerTime = ComputerTime - datenum(0,0,0,5,0,0);   % for local time
     
     
     Compressed = [Compressed'; CompressedVoltage']';
     Time = [Time ComputerTime];
     % Time(end+1,:) = ComputerTime;
    
     
     
end

%% -----------------------------------------------------------------    
    
    % 4 freqs
    switch ID
        case 'ES70'
   
            % CLimDB = [-120 0];
            CLimDB = ([-80 -40]);
            
        case 'ES120' 
        
             CLimDB = [-130 0];    
             
        case 'ES200' 
          
             % CLimDB = [-120 0];     
             CLimDB = ([-95 -55]);
             
        case 'ES333' 
             
             % CLimDB = [-140 0];
             CLimDB = ([-105 -65]) ;
             
        otherwise
            fprintf('\t ERROR: ID %s not found', ID)
    end
    
   %% ------------------------------------------------------------
    
   
    figure(202); clf
    if plots_up                 % only do this once
        pos = get(202,'position');
        set(202,'position',[pos(1) pos(2) pos(3)*2 pos(4) ])
        plots_up = 0;
    end
    
    L = length(Time);
   pcolor(Time, Range,20*log10((abs(Compressed(:,1:L))))); shading flat; 
   % pcolor(20*log10((abs(Compressed(:,1:end))))); shading flat;  
    
    %  sometimes the CompressedVoltage has 1 more entry, check so doesn't
    %  crash
%     if size(CompressedVoltage,2) == length(ComputerTime)
%         pcolor(ComputerTime, Range,20*log10((abs(CompressedVoltage(:,1:end))))); shading flat; 
%     else
%          N = [1:length(ComputerTime)];
%          N = min(N,size(CompressedVoltage,2));
%         pcolor(ComputerTime, Range,20*log10((abs(CompressedVoltage(:,N))))); shading flat;
%     end
%     
    if Intrp_shading
        shading interp;
    end
    
    % Hardwired to be on same scale as Rocky!
    tmin = datenum(2018,11,29,10,22,0);
    tmax = datenum(2018,11,29,10,48,48);  
    set(gca,'xlim',[tmin tmax]);
    
   
    
    % %  pcolor(ComputerTime, Range,20*log10((abs(CompressedVoltage(:,1:end-1))))); 
    datetick('x','keeplimits');
    axis('ij');
   
    Ylims = get(gca,'ylim');
    ylim([0 floor(Ylims(2))]);
    
    caxis(CLimDB)
    
    
    axis('ij')   % % pings downward
    hcb = colorbar;
    o_ax = get(gca,'position');
    hcb.Position = [.82 .11 .04  .5];
  
    hcb.Position = [o_ax(1) + o_ax(3)+.015 .11 .03  .5];
    
    % New color scheme
    % CMap = cmap_ek80;
    colormap(cmap_ek80);
    
    
    % caxis([-84 0])
    % hcb.Limits = [-80 0];
     
    set(gca, 'position',o_ax);
    set(get(hcb,'title'),'string','dB');
    
    
    title(['James River '  ID],'interpreter','none');
    S = sprintf('Time (local %s)', datestr(ComputerTime(1),'mm/dd/yy'));
    xlabel(S); ylabel('Range (m)')
    
    drawnow
    
    if save_plots
        figure(202);
        fprintf( '\t Printing: %s\n', plotFilename );
        print('-dpng', plotFilename);
    end
  

    
    fprintf('\tTime to here %.2f mins \n',toc/60 );
  %% --------------------------------------------------------------
      
