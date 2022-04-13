clear
close all
clc

%% Set-up signal
dec=1;
T=1e-3;
fs=1.5e6;
fi=95e3;
fband=100e3;
ff=fi+fband;

%% Transmit Signal
t=0:1./fs:(T-1./fs);
tx=chirp(t,fi,T,ff);
tx1=0.5*chirp(t,fi+10e3,T,ff-10e3)+tx;

t=0:1./fs:(2*T-1./fs);
tx=[tx zeros(size(tx))];
tx1=[tx1 zeros(size(tx1))];

tx=hilbert(tx);
tx1=hilbert(tx1);

%t=decimate(t,dec);
%tx=decimate(tx,dec);
%tx1=decimate(tx1,dec);

%% Receive Signal
%rx=1*(tx+0.6*randn(size(tx)));
rx=1*((tx1)+0.3*randn(size(tx1)));

figure(1)
plot(t,real(tx))
ylim([-1.1 1.1])
hold on
plot(t,real(rx),'r')
hold off

%% FFT of transmit signal
n=4096*2;%2.^13
n=2.^11;
FFT_tx=fft(tx,n);
f=(0:length(FFT_tx)-1)*fs/n/dec;

figure(2)
%plot(f/1000,20.*log10(abs(FFT_tx)./max(abs(FFT_tx))))
plot(f/1000,20.*log10(abs(FFT_tx)))

%% FFT of receive signals
FFT_rx=fft(rx,n);

figure(2)
hold on
%plot(f/1000,20.*log10(abs(FFT_rx)./max(abs(FFT_rx))),'r')
plot(f/1000,20.*log10(abs(FFT_rx)),'r')
hold off

%% Cross-correlataion of transmit and receive signals
MF=xcorr(rx,tx);
figure(3)
plot(abs(MF)./max(abs(MF)),'b')

%% Pick out peak of cross-correlation
[tmp,ij]=max(MF);clear tmp
ii=(ij-150):(ij+150);
MF_peak=MF(ii);clear ij

figure(3)
hold on
plot(ii,(abs(MF_peak)./max(abs(MF_peak))),'k');clear ij
hold off;clear ii

%% FFT of cross-correlated signal
FFT_MF=fft(MF_peak,n);
f_MF=(0:length(FFT_MF)-1)*fs/n/dec;

figure(2)
hold on
%plot(f_MF/1000,20.*log10(abs(FFT_MF)./max(abs(FFT_MF))),'g')
plot(f_MF/1000,20.*log10(abs(FFT_MF))-37.5,'g')
hold off

%% Convolution of transmit and receive signals
CP=conv(rx,conj(fliplr(tx)));

figure(3)
hold on
plot((abs(CP)./max(abs(CP))),'r')
hold off

%% Pick out peak of convolution
[tmp,ij]=max(CP);clear tmp
ii=(ij-150):(ij+150);
CP_peak=CP(ii);clear ij

figure(3)
hold on
plot(ii,(abs(CP_peak)./max(abs(CP_peak))),'m');clear ij
hold off;clear ii


%% FFT of convolved signal
FFT_CP=fft(CP_peak,n);
f_CP=(0:length(FFT_CP)-1)*fs/n/dec;
figure(2)
hold on
%plot(f_CP/1000,20.*log10(abs(FFT_CP)./max(abs(FFT_CP))),'g')
plot(f_CP/1000,20.*log10(abs(FFT_CP))-37.5,'k')
hold off

whos
