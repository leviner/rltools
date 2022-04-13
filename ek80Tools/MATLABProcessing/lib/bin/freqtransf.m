function [fvec, FFTvec ] = freqtransf( FFTvecin,fsdec,fc)

%   Detailed explanation goes here

nfft        = length(FFTvecin);
FFTvec      = [FFTvecin; FFTvecin; FFTvecin];
fvec        = fsdec*linspace(0,1-1/nfft,nfft)';
fvec        = [fvec;fsdec+fvec;2*fsdec+fvec];

if (fc>fsdec/2)
    idxmin      = round((fc-fsdec/2)/fsdec*nfft);
else
    idxmin      = 1;
end
%fprintf('The frequency shift is %d.\n',round(fvec(idxmin)))

FFTvec      = FFTvec(idxmin:idxmin+nfft-1);
fvec        = fvec(idxmin:idxmin+nfft-1);

end