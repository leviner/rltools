clear all; close all; clc
addpath bin

prompt = {'Enter bin width in pings','Enter bin depth in meters'};
test = inputdlg(prompt);

win.l = str2num(test{1});       % window length in pings
win.step = str2num(test{2});   % window step size in meters
win.nfft = 2^9; % current FFT size = 512 points (about twice as long as number in longest window)
n = win.nfft;
bar = waitbar(0,'Getting ready...') ;
[fn, filepath] = uigetfile('*.mat','Pick a raw data file','MultiSelect','on');
outdir = uigetdir(pwd,'Select Directory for Results');
if isstr(fn), fn={fn}; end  % convert char string to cellstr
load([filepath '\' fn{1}])
%%

[nChannels,nPings] = size(data.echodata);
startPings = [1:win.l:nPings];

maxr = [];
for i=1:nChannels
    maxr = [maxr data.echodata(i,1).maxRange];
end
minr = 0;maxr = min(maxr);
win.r = [minr:win.step:maxr];
win.meanrange = win.r+(win.r(2)-win.r(1))/2;
%%
bar = waitbar(0,'Getting ready...') ;
for iii = 1:length(fn)
    load([filepath '\' fn{iii}])
    clear CVAll
    for jjj= 1:nChannels
        clear CV
        t = data.echodata(jjj,:);
        if isempty(t(1).compressed)
            continue
        else
        for s =1:nPings
            CV(:,s) = mean(t(s).compressed,2);
        end
        CVAll{jjj,1} = [CV];
        end
    end

    clear Spec f F
    for j = 1:length(startPings) % for each ping
        for jj = 1:length(win.r) % for each range bin
            waitbar(j/length([1:win.l:nPings]),bar,['Calculating spectra for each ' num2str(win.l) ' ping by ' num2str(win.step) ' m bin\n For file ' iii ' of ' length(fn)]) ;
            for jjj = 2:nChannels % for each channel
                if isempty(CVAll{jjj,1})
                    svtmp{jjj} = NaN;
                    f{jjj} = NaN;
                    continue
                end
                
                timestamp{jjj,j}=data.echodata(jjj,startPings(j)).timestamp;
                
                if startPings(j)+win.l-1 < nPings
                    CVwin = mean(CVAll{jjj,1}(:,startPings(j):startPings(j)+win.l-1),2);
                else
                    CVwin = mean(CVAll{jjj,1}(:,startPings(j):nPings),2);
                end
                
                CVwinR = CVwin((data.echodata(jjj,1).range < win.r(jj)+win.step) &(data.echodata(jjj,1).range >= win.r(jj)));
                specvec = CVwinR.*data.echodata(jjj,1).range(((data.echodata(jjj,1).range < win.r(jj)+win.step)&(data.echodata(jjj,1).range >= win.r(jj))));
                b = tukeywin(length(specvec),0.1)./(norm(tukeywin(length(specvec),0.1))./sqrt(length(specvec)));
                specvec = specvec.* b;
                specvec = fft(specvec,n);                
                
                fsdec = (data.param(jjj, 1).FrequencyStart+data.param(jjj, 1).FrequencyEnd)/2;
                
                
                if isstr(data.config.transceivers(jjj).channels.transducer.Frequency)
                    fnom = str2num(data.config.transceivers(jjj).channels.transducer.Frequency);
                else
                    fnom = data.config.transceivers(jjj).channels.transducer.Frequency;
                end
                
                [ftmp, FFTvec_tmp] = freqtransf(specvec,fsdec,fnom);
                alphaf =  alpha_sea(data.environ.Depth,data.environ.Salinity,data.environ.Temperature,data.environ.Acidity,fsdec/1000);
                calf = data.calibration(jjj).Frequency;
                calg = data.calibration(jjj).Gain;
                dens = gsw_rho(data.environ.Salinity,data.environ.Temperature,data.environ.Depth);
                c = gsw_sound_speed(data.environ.Salinity,data.environ.Temperature,dens*9.81*data.environ.Depth*1e-4);
                
                if isstr(data.config.transceivers(jjj).channels.transducer.EquivalentBeamAngle)
                    calpsi = str2num(data.config.transceivers(jjj).channels.transducer.EquivalentBeamAngle);
                else
                    calpsi = data.config.transceivers(jjj).channels.transducer.EquivalentBeamAngle;
                end
                
                G = interp1(calf,calg,ftmp);
                psi =calpsi + 20*log10(fnom./ftmp);
                dt =2*((win.r(jj)+win.step)-win.r(jj))/c;
                pr = abs(FFTvec_tmp).^2;
                zet = data.parameters.Ztrd;
                zer = str2num(string(data.config.transceivers(jjj).Impedance)); % transceiver imedpance
                P_tr = data.param(jjj,1).TransmitPower;
                sv = 10*log10(pr) + ...
                    2.*alphaf.*win.meanrange(jj) - 2.*G - psi - ...
                    10*log10(dt) + ...
                    10*log10(4./zet./P_tr./(2*sqrt(2)).^2) +...
                    10.*log10((zer+zet)/zer) - ...
                    10.*log10(c^3./(32.*pi^2.*ftmp.^2));
                svtmp{jjj} = sv';
                f{jjj} = ftmp';
                
                
            end
            Spec{jj,j} = [svtmp{:}];
            F{jj,j} =  [f{:}];
            
        end
    end
 
fout = [outdir '\Spectra_' char(fn(iii))];
save(fout, 'Spec','F','win','timestamp')
    
end
waitbar(1,bar,'Done') ;