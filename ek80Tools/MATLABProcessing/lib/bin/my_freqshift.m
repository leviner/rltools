function [fvecout, FFTvec]=my_freqshift(FFTin,fsdec,op);

nfft        = length(FFTin)
FFTvec      = [FFTin; FFTin; FFTin];
fvec        = fsdec*linspace(0,1-1/nfft,nfft)';
fvec        = [fvec;fsdec+fvec;2*fsdec+fvec];

f1=op.fc;
f2=op.fstart;
f3=(op.fstart+op.fstop)./2;

round((op.fc-fsdec/2)/fsdec*nfft);
round((op.fstart-fsdec/2)/fsdec*nfft);
round(((op.fstart+op.fstop)./2-fsdec/2)/fsdec*nfft);

f=f3;

if (f>fsdec/2)
    idxmin      = round((f-fsdec/2)/fsdec*nfft);
else
    idxmin      = 1;
end

FFTvec      = FFTvec(idxmin:idxmin+nfft-1);
fvecout        = fvec(idxmin:idxmin+nfft-1);

end