%load('C:\Users\cbassett\Desktop\GOABB\Processed\Data\Sv_DY2104-D20210620-T053246.mat')

clear all; close all; clc
% read in broadband GOA data and calculate Sv(f) in small spatial bins 
% with window lenghts, overlaps, etc. specified below
addpath bin
win.l = 1;       % window length in pings
win.step = 0.5; % window step size in meters
win.nfft = 2^9; % current FFT size = 512 points (about twice as long as number in longest window)
n = win.nfft;
load('D:\FMProcessing\EK80_matlab_software_2021\testin\GOA_BB_Cals.mat')
%load('C:\Users\cbassett\Desktop\GOABB\Processed\Data\Sv_DY2104-D20210620-T053246.mat')
fn = dir('D:\FMProcessing\EK80_matlab_software_2021\test\out\Data\*.mat')
load([fn(1).folder '\' fn(1).name])

dr = [diff(range{1}(1:2))];% diff(range{2}(1:2)) diff(range{3}(1:2)) diff(range{4}(1:2))];
maxr = [max(range{1})];% max(range{2}) max(range{3}) max(range{4})];
maxr = min(maxr); minr = 0; 
win.r = [minr:win.step:maxr];
win.rplot = win.r + win.step/2;
win.ninds = ceil(win.step./dr);

%%

% start loop for reload
%load('C:\Users\cbassett\Desktop\GOABB\Processed\Data\Sv_DY2104-D20210620-T053246.mat')
%load('C:\Users\cbassett\Desktop\GOABB\Processed\Data\Sv_DY2104-D20210620-T051523_FM_split.mat')

fn = dir('D:\FMProcessing\EK80_matlab_software_2021\test\out\Data\*.mat');

for iii = 1:length(fn)
load([fn(iii).folder '\' fn(iii).name]);    
clear Spec f   
nopings = length(CV{1}(1,:));
for j = 1%:nopings

    if j == 1 % on first loop create channel inds
        for ii = 1:length(win.r)
        [~,wind{1}(ii)] = min(abs(range{1}- win.r(ii)));
        %[~,wind{2}(ii)] = min(abs(range{2}- win.r(ii)));
        %[~,wind{3}(ii)] = min(abs(range{3}- win.r(ii)));
        %[~,wind{4}(ii)] = min(abs(range{4}- win.r(ii)));
        end  
    end
    
    for jj = 1%:length(win.r)
        
        
        
        for jjj = 1%:length(CV)
            st = wind{jjj}(jj); et = st+win.ninds(jjj);
            %y = zeros(win.nfft,1);
            %y(1:win.ninds(jjj)+1) = CV{1,jjj}(st:et,j);
            
            %specvec = zeros(1,n)  ;  
                 %   rtmp = ping.rind{jj,jcnt}(jjj);   % start window range - temp
                    p = j;        % ping number - temp
                    
                    specvec = CV{jjj}(1110:1200,1).* range{jjj}(1238:1328,1)%(wind{jjj}(jj):wind{jjj}(jj)+win.ninds(jjj)-1,j).* ...
                      %  range{jjj}(wind{jjj}(jdj):wind{jjj}(jj)+win.ninds(jjj)-1);
                    %b = tukeywin(n,0.1);
                    b = tukeywin(91,0.1)./(norm(tukeywin(91,0.1))./sqrt(91));
                    specvec = specvec.* b;
                    meanrange = win.rplot(jj); %%%%%%%%%%%%%%%%%%
                    specvec = fft(specvec,n);

                    [ftmp, FFTvec_tmp] = freqtransf(specvec,para{jjj}.fsdec,para{jjj}.fnom);
                    alphaf = alpha_sea(w.z,w.S,w.T,w.pH,ftmp./1000);
                    if jjj == 1
                       calf = Cal.f70;
                       calg = Cal.G70;
                       calpsi = para{1}.psinom;
                       calpsi = calpsi + 20*log10(70000./ftmp);
                    elseif jjj == 2
                       calf = Cal.f70;
                       calg = Cal.G70;   
                       calpsi = para{2}.psinom;
                       calpsi = calpsi + 20*log10(70000./ftmp);

                    elseif jjj == 3
                       calf = Cal.f120;
                       calg = Cal.G120;  
                       calpsi = para{3}.psinom;
                       calpsi = calpsi + 20*log10(120000./ftmp);
                    elseif jjj == 4
                       calf = Cal.f200;
                       calg = Cal.G200;  
                       calpsi = para{4}.psinom;
                       calpsi = calpsi + 20*log10(200000./ftmp);
                    end
                    G = interp1(calf,calg,ftmp);
                   % psi = interp1(calf,calpsi,ftmp,'linear','extrap');
                    psi = calpsi;
                    dt =2*( range{jjj}(wind{jjj}(jj)+win.ninds(jjj) ) - range{jjj}(wind{jjj}(jj)) )/w.c;
                    pr = abs(FFTvec_tmp).^2;
                    zet = para{jjj}.zet;
                    zer = para{jjj}.zer;
                    P_tr = para{jjj}.ptx;

                    svtmp{jjj} = 10*log10(pr) + ...
                                   2.*alphaf.*meanrange - 2.*G - psi - ...
                                   10*log10(dt) + ...
                                   10*log10(4./zet./P_tr./(2*sqrt(2)).^2) +...
                                   10.*log10((zer+zet)/zer) - ...
                                   10.*log10(w.c^3./(32.*pi^2.*ftmp.^2));
%                                svtmp(:,cnt2) = 10*log10(pr) + ...
%                                    2.*alphaf.*meanrange - 2.*G - psi - ...
%                                    10*log10(dt) + ...
%                                    10*log10(4./zet./P_tr./(2*sqrt(2)).^2) +...
%                                    10.*log10((zer+zet)/zer) - ...
%                                    10.*log10(w.c^3./(32.*pi^2.*ftmp.^2));
                    
         
            
            if j == 1
                % create a full frequency vector for storage
                if jjj == 1
                f = ftmp';
                else
                f = [f, NaN, ftmp'];    
                end  
            end
        
            if jjj == 1
              %fu70 = find(and(f(1:n) > 34999, f(1:n) < 43000.1)) ;
              fu70 =  f;%find(and(f(n+2:2*n+2) > 45999, f(n+2:2*n+2) < 86001));
              %fu120 = find(and(f(2*n+3:3*n+3) > 94999, f(2*n+3:3*n+3) < 157001));
              %fu200 = find(and(f(3*n+4:end) > 164999, f(3*n+4:end) < 255001));
              Spec{jj,j} = NaN(size(f));
              
              %svtmp38 = NaN(1,512); svtmp70 = NaN(1,512); svtmp120 = NaN(1,512); svtmp200 = NaN(1,512);
              %svtmp38(fu38) = svtmp{1}(fu38);
             % svtmp70 = NaN(1,512)
              svtmp70 =svtmp{1};
              %svtmp120(fu120) = svtmp{3}(fu120);
              %svtmp200(fu200) = svtmp{4}(fu200);
 
          %    Spec{jj,j} = [svtmp{1}(fu38); NaN; svtmp{2}(fu70); NaN; svtmp{3}(fu120); NaN; svtmp{4}(fu200)];    
                         Spec{jj,j} = [svtmp70'];    
     
            end  
        %Spectmp = [spectmp{1}', NaN, spectmp{2}', NaN, spectmp{3}', NaN, spectmp{4}']; 
        %Spec{j,jj} = Spectmp;   
    
    end
    
    end
j
end
 
fout = ['D:\FMProcessing\EK80_matlab_software_2021\test\out\Data\' fn(iii).name(1:end-4) '_Spectra.mat'];
r = win.rplot; 
notes = 'r - the middle of the window used for the processing';
save(fout, 'Spec','f','r','win','notes')

end
%%

%for j = 20:10:70
figure
subplot(211)
% use for initial data
%     plot(f./1000,Spec{56,39},'b'), hold on
%     plot(f./1000,Spec{56,40},'r')
%     plot(f./1000,Spec{56,41},'b')
%     plot(f./1000,Spec{56,42},'c')
% 
%    spec = 10.*log10(10.^(Spec{56,39}./10)+10.^(Spec{56,40}./10)+10.^(Spec{56,41}./10)+10.^(Spec{56,42}./10));
    plot(f./1000,Spec{70,11},'b'), hold on
    plot(f./1000,Spec{70,12},'r')
    plot(f./1000,Spec{70,13},'b')
    plot(f./1000,Spec{70,14},'c')

   spec = 10.*log10(10.^(Spec{70,30}./10)+10.^(Spec{70,31}./10)+10.^(Spec{70,32}./10)+10.^(Spec{70,33}./10)); 
hold on
   plot(f./1000,spec,'k','linewidth',3)
   ylabel('Sv [dB re 1/m]')
   xlabel('Freq [kHz]')
   axis([0 300 -70 -20])
%end
title('Fish')

subplot(212)
    plot(f./1000,Spec{47,15},'b'), hold on
    plot(f./1000,Spec{47,12},'r')
    plot(f./1000,Spec{47,13},'b')
    plot(f./1000,Spec{47,14},'c')

   spec = 10.*log10(10.^(Spec{47,11}./10)+10.^(Spec{47,12}./10)+10.^(Spec{47,13}./10)+10.^(Spec{47,14}./10));
   hold on
   plot(f./1000,spec,'k','linewidth',3)
   ylabel('Sv [dB re 1/m]')
   xlabel('Freq [kHz]')
      axis([0 300 -120 -60])
title('Krill')
%end
%%
figure(2)
pcolor(1:length(Sv{1}(1,:)),-range{1},Sv{1}), shading flat
set(gca,'clim',[-70 -40])

%%
   plot(f./1000,spec,'k','linewidth',3)
