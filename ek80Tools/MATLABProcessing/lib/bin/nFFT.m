function [Fft,f] = nFFT(data,Fs)

n=4096*2;
FFTX=fft(data,n);

NumUniquePts = ceil((n+1)/2);
f=(0:NumUniquePts-1)*Fs/n;
FFTX=FFTX(1:NumUniquePts,:);
Fft=abs(FFTX)/length(data);