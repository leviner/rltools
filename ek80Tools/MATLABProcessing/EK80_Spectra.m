clear all; close all; clc
addpath lib
addpath('lib/bin')
addpath('lib/bin/gsw')

prompt = {'Enter window width in pings','Enter window width overlap in pings','Enter window depth in meters',};
test = inputdlg(prompt);

win.l = str2num(test{1});       % window length in pings
win.overlap = str2num(test{2});  % window overlap in pings
win.step = str2num(test{3});   % window step size in meters
win.nfft = 2^10; % current FFT size = 1024 points (about twice as long as number in longest window)


[fn, filepath] = uigetfile('*.mat','Pick a raw data file','MultiSelect','on');
outdir = uigetdir(pwd,'Select Directory for Results');
if isstr(fn), fn={fn}; end  % convert char string to cellstr
load([filepath '\' fn{1}])
%%

[nChannels,nPings] = size(data.echodata);
startPings = 1:win.l-win.overlap:nPings;

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
        if isempty(t(1).complexsamples)
            continue
        else
            for s =1:nPings
                CV(:,s) = mean(t(s).complexsamples,2);
            end
            CVAll{jjj,1} = [CV];
        end
    end

    clear Spec f F SpecBins
    for j = 1:nPings % for each ping
        for jj = 1:length(win.r) % for each range bin
            waitbar(j/length([1:win.l:nPings]),bar,['Calculating spectra for each ' num2str(win.l) ' ping by ' num2str(win.step) ' m bin ' newline 'for file ' num2str(iii) ' of ' num2str(length(fn))]) ;
            for jjj = 1:nChannels % for each channel
                ranges = data.echodata(jjj,1).range+(0-min(data.echodata(jjj,1).range));
                if data.param(jjj,1).PulseForm == 0
                    svtmp{jjj} = NaN;
                    f{jjj} = NaN;
                    continue
                end

                timestamp{jjj,j}=data.echodata(jjj,j).timestamp;

                CVwin = CVAll{jjj,1}(:,j); % grab the ping

                CVwinR = CVwin((ranges < win.r(jj)+win.step) &(ranges >= win.r(jj)));
                specvec = CVwinR.*ranges(((ranges < win.r(jj)+win.step)&(ranges >= win.r(jj))));
                b = tukeywin(length(specvec),0.1)./(norm(tukeywin(length(specvec),0.1))./sqrt(length(specvec)));
                specvec = specvec.*b;
                specvec = fft(specvec,win.nfft);

                fsdec = 1/data.param(jjj, 1).SampleInterval;


                if isstr(data.config.transceivers(jjj).channels.transducer.Frequency)
                    fnom = str2num(data.config.transceivers(jjj).channels.transducer.Frequency);
                else
                    fnom = data.config.transceivers(jjj).channels.transducer.Frequency;
                end

                [ftmp, FFTvec_tmp] = freqtransf(specvec,fsdec,fnom);
                alphaf =  alpha_sea(data.environ.Depth,data.environ.Salinity,data.environ.Temperature,data.environ.Acidity,ftmp/1000);
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
                if j == 1
                    Gs{jjj}=G;
                end

            end
            Spec{jj,j} = [svtmp{:}];
            F{jj,j} =  [f{:}];


        end
    end

    % Now we make the average for the ping (horizontal) bins
    waitbar(1,bar,'Calcuating horizontal bin averages...');
    for p = 1:length(startPings)
        if startPings(p)+win.l-1 > nPings
            test = Spec(:,startPings(p):nPings);
        else
            test = Spec(:,startPings(p):startPings(p)+win.l-1);
        end

        for tt= 1:size(test,1)
            curBin = test(tt,:);
            h = 0;
            for v = 1:length(curBin)
                h = h+10.^(curBin{v}./10);
            end
            h2 = h./length(curBin);
            h2 = 10.*log10(h2);

            SpecBins{tt,p} = h2;
        end
    end

    waitbar(1,bar,['Saving file '  num2str(iii) ' of ' num2str(length(fn))]);
    fout = [outdir '\Spectra_Complex_' num2str(win.l)...
        'p' num2str(win.step) 'm' num2str(win.overlap) 'overlap_' char(fn(iii))];
    save(fout, 'Spec','SpecBins','F','win','timestamp')

end
waitbar(1,bar,'Done') ;