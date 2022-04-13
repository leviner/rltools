function [Psi,BP,freq] = Beam_Pattern_ftheta(sys_indx, theta_c, plot_it);
% compute beam pattern bp as a function of frequency at a particular theta_c
% and compute equivalent beam angles for Airmar and Simrad transducers

%clear 
%close all
%hold on

%Channel
%sys_indx=3;    
% 1 = Airmar 25-45 kHz, "38 kHz"
% 2 = Simrad 70kHz, 55-90 kHz, ES70-7CD, two-way beam angle = -21dB
% 3 = Simrad 120kHz, 90-160 kHz, ES120-7CD
% 4 = Simrad 200kHz, 160-260 kHz, ES200-7CD
% 5 = Simrad 333kHz, 260-410 kHz, ES333-7CD

%Sound Speed
cw=1500;

switch sys_indx
    case 1
        f1=25000;
        f2=55000;
        fc=38000;
        a=87.3e-3./2;       % radius
        sys_str='A38';
    case 2
        f1=45000;
        fc=70000;
        f2=95000;
        a=185e-3./2;      % radius
        sys_str='S70';
    case 3
        f1=90000;
        fc=120000;
        f2=160000;
        a=112e-3./2;      % radius
        sys_str='S120';
   case 4
        f1=160000;
        fc=200000;
        f2=260000;
        a=60e-3./2;      % radius
        sys_str='S200';
    case 5
        f1=260000;
        fc=333000;
        f2=420000;
        a=90e-3./2;        % radius of A1 AirMar M159 33.5 kHz 12 degree half beamwidth (nominal center frequency 41 kHz?)
        sys_str='S333';
end

%Number of frequencies
nf=181;
%Frequency
freq=linspace(f1,f2,nf);
%Number of angles for equivalent beam angle calculation
nth=400;
%Angles for equivalent beam angle calculation
th=linspace(0,pi/3,nth);
dth=th(2)-th(1);

Threshold=-30;   % in dB
threshold=10.^(Threshold/10);

for ll=1:nf                     % frequency loop
  f=freq(ll);
  k=2*pi*f/cw;
  arg=k*a*sin(th)+eps;
  arg1=k*a*sin(theta_c*pi./180)+eps;

  %% 
  bp(ll)=abs(2*besselj(1,arg1)./arg1).^4;

  %% composite beam pattern (TX & RCV)
  bp_c=abs(2*besselj(1,arg)./arg).^4;
  
  %% equivalent beam angle
  indx_ij=find(bp_c > threshold);
  psi(ll)=2*pi*nansum(bp_c(indx_ij).*sin(th(indx_ij)))*dth;
end

BP=10.*log10(bp);
Psi=10*log10(psi);
if plot_it==1
    figure(1)
    plot(freq*1e-3,BP,'linewidth',1.5)
    xlabel('Frequency (kHz)','fontsize',12,'fontweight','bold')
    figure(2)
    hold on
    plot(freq*1e-3,Psi,'linewidth',1.5)
    xlabel('Frequency (kHz)','fontsize',12,'fontweight','bold')
end

[min_tmp,ij]=min(abs(freq-fc));
freq(ij);
Psi(ij)
clear ij
clear min_tmp