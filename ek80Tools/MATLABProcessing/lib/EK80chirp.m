function [t, y] = EK80chirp(txpower, fstart, fstop, slope, tau, z, fs)
        
    % z is the tranducer quadrant nominal impedance [Ohms]
    % fs is the sampling frequency [Hz]
    
    % WBT signal generator
    dt = 1/fs;  % sample interval
    a  = sqrt((txpower/4) * (2*z));
    
    % Create transmit signal
    t  = (0:dt:tau-dt)';
    nt = length(t);
    
    nwtx    = 2*floor(slope*nt);
    wtxtmp  = hann(nwtx);
    nwtxh   = round(nwtx/2);
    wtx     = [wtxtmp(1:nwtxh); ones(nt-nwtx,1); wtxtmp(nwtxh+1:end)];

    y = a * chirp(t, fstart, tau, fstop) .* wtx;
    
    % The transmit signal needs to have a max amplitude of 1, as per
    % the documentation from Lars.
    y = y / max(abs(y));

end
