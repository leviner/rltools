% function  plot_processed_files_GUI_select
%
%   plot processed files in matlab output directory
%       mean plus all the quandrants
%
%
%   files hardwired at the moment. Run as script.
%    
%  Quick and dirty....  for first looks
%
%   Have to run init program first to set working space
%
%    User options/flags can be set first
% 
% -----------------------------------------------------------
 
close all;  

% User options ----------------------------------------------
 
save_plots = 1;     % save the plots,  or not...
GUIselect = 1;      %  This flag (1) opens GUI 
Intrp_shading = 1;  %  0: shading flat, 1: shading interp

% -----------------------------------------------------------------
% Select files 
% ---------------------------------------------------------------
clear Files;

if GUIselect
    [filename, filepath] = uigetfile({'*.mat'},'Pick matlab ES files to plot','MultiSelect','on');
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
    
    Files = dir([MatOutDir '*ES*.mat']);
    
end


tic;

plots_up = 1;    % for making long plot window, must keep if using if below

if(~exist('MatOutDir'))
    fprintf('\t ERROR: Initialize by running init file\n');
end


% --------------------------------------------------------------------------
 

first = 1;
proc_files = 1:length(Files);

for ifile = proc_files       
    
    filename = Files(ifile).name;
    load([MatOutDir '/' filename]);
    fprintf('Working on %d of %d - %s\n',ifile, length(Files), filename);
    
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
    
      % 4 freqs
    switch ID
        case 'ES70'
   
            ClimDB = [-120 0];
            CLimDB = ([-80 -40]);
            
        case 'ES120' 
        
             ClimDB = [-130 0];    % curious, looks best
             
        case 'ES200' 
          
             ClimDB = [-120 0];    % curious, looks best
             CLimDB = ([-95 -55]);
             
        case 'ES333' 
             
             ClimDB = [-140 0];
             ClimDB = ([-105 -65]) ;
             
        otherwise
            ClimDB = [-70 -30];
            fprintf('\t ERROR: ID %s not found', ID)
    end
   
    
    figure(100); clf
    plot(Range,20*log10((abs(CompressedVoltage(:,1)'))),'k')
    title(filename,'interpreter','none');
    xlabel('Range (m)');
    grid on;
    
    figure(101); clf 
    Range = data.echodata(ID_ndx).range;
    plot(Range,20*log10(abs(data.echodata(ID_ndx).complexsamples)))
    title(filename,'interpreter','none');
    xlabel('Range (m)');
    grid on;
    
    hold on 
    plot(Range,20*log10((abs(CompressedVoltage(:,1)'))),'k')
    hold off
    
    
    % james river is local = UTC -5
    ComputerTime = ComputerTime% - datenum(0,0,0,5,0,0);   % for local time
    
    
    figure(102); clf
    if plots_up                 % only do this once
        pos = get(102,'position');
        set(102,'position',[pos(1) pos(2) pos(3)*2 pos(4) ])
        plots_up = 0;
    end
    
    %  sometimes the CompressedVoltage has 1 more entry, check so doesn't
    %  crash
    if size(CompressedVoltage,2) == length(ComputerTime)
        pcolor(ComputerTime, Range,20*log10((abs(CompressedVoltage(:,1:end))))); shading flat; 
    else
         N = [1:length(ComputerTime)];
         N = min(N,size(CompressedVoltage,2));
        pcolor(ComputerTime, Range,20*log10((abs(CompressedVoltage(:,N))))); shading flat;
    end
    
    if Intrp_shading
        shading interp;
    end
   
    
    % %  pcolor(ComputerTime, Range,20*log10((abs(CompressedVoltage(:,1:end-1))))); 
    datetick('x');
    axis('ij');
   

    % for %ES200 file #2 ????
    % pcolor(ComputerTime, Range,20*log10((abs(CompressedVoltage(:,1:end))))); shading flat     
    % datetick('x','keeplimits');
    datetick('x','keepticks');
    
    Ylims = get(gca,'ylim');
    ylim([0 floor(Ylims(2))]);
    
    caxis(ClimDB)
    
    
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
    % set(get(hcb,'title'),'string','dB');
    
    
    title(filename,'interpreter','none');
    S = sprintf('Time (local %s)', datestr(ComputerTime(1),'mm/dd/yy'));
    xlabel(S); ylabel('Range (m)')
    
    drawnow
    
    if save_plots
        figure(102);
        fprintf( '\t Printing: %s.png\n', filename(1:end-4));
        print( '-dpng',[ PlotOutDir filename(1:end-4) '.png']);
    end
    
    % ------------------------------------------------------------------
    % checking for good GPS data in files.  Optional
    gps_check = 0;
    if (gps_check && ~isempty(GPStime) )
        fprintf('\tGPS nmea data found...\n');
    end
    
    % get time left and print
    
    if first
        Proc1 = toc/60;   % time to process 1 image in minutes
        first = 0;
    end
    
    T_left = (length(proc_files) * Proc1) - toc/60;
    if(T_left <= 0), T_left = 0; end
    
    fprintf('\tTime to here %.2f mins, %.2f mins left\n',toc/60,T_left);
    
end
