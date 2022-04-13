clear
%close all
clc


%% Set-up signal
dec=2;
T=1e-3;
fs=1.5e6/1;
fi=95e3;
fband=100e3;
ff=fi+fband;


%% Transmit Signal
t=0:1./fs:(T-1./fs);
tx=chirp(t,fi,T,ff);
tx1=0.5*chirp(t,fi+10e3,T,ff-10e3)+tx;
tx2=tx+0.1*randn(size(tx));

t=0:1./fs:(2*T-1./fs);
tx=[tx zeros(size(tx))];
tx1=[tx1 zeros(size(tx1))];
tx2=[tx2 zeros(size(tx2))];

%tx=hilbert(tx);
%tx1=hilbert(tx1);

t=decimate(t,dec);
tx=decimate(tx,dec);
tx1=decimate(tx1,dec);
tx2=decimate(tx2,dec);

%% Receive Signal
%rx=1*(tx+0.6*randn(size(tx)));%
%rx=1*((tx1)+0.3*randn(size(tx1)));
rx=tx2;

figure
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

figure
%plot(f/1000,20.*log10(abs(FFT_tx)./max(abs(FFT_tx))))
plot(f/1000,20.*log10(abs(FFT_tx)))


%% FFT of receive signals
FFT_rx=fft(rx,n);
hold on
%plot(f/1000,20.*log10(abs(FFT_rx)./max(abs(FFT_rx))),'r')
plot(f/1000,20.*log10(abs(FFT_rx)),'r')
hold off