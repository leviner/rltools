function ytx = EK80calcSentSignal_DT(data, channel, settingsPing)
    %
    % Generates a simulated sent pulse waveform using the parameters from
    % the given channel and ping number.
    %
    %   Revised for CTRiver - Perkins
        
    % generate transmitted pulse
    if isfield(data.param,'PulseLength')
        tx.tau = data.param(channel, settingsPing).PulseLength;
    elseif isfield(data.param,'PulseDuration')
        tx.tau = data.param(channel, settingsPing).PulseDuration;
    end
    
    tx.sf  = 1/data.param(channel, settingsPing).SampleInterval;
    tx.f0  = data.param(channel, settingsPing).FrequencyStart;
    tx.f1  = data.param(channel, settingsPing).FrequencyEnd;
    
    % Transmit pulse output parameters
    [tx.t, tx.ey] = EK80chirp(data.param(channel, settingsPing).TransmitPower, ...
        tx.f0, tx.f1, data.param(channel, settingsPing).Slope, tx.tau, ...
        data.parameters.Ztrd, data.parameters.fs);

    % Filter and decimate to create matched signal
    figure
    dt=tx.t(2)-tx.t(1);
    T0=length(tx.ey)*dt;
    df0=1e-3/T0;
    f0=[0:length(tx.ey)-1]*df0;
    % % plot(f0,20*log10(abs(fft(tx.ey))),'b') % original
    plot(f0,20*log10(fftshift(abs(fft(tx.ey)))),'b')
    
    
    % WBT filter and decimation
    ytx  = conv(tx.ey,     data.filters(channel, 1).FilterData);
    dt=tx.t(2)-tx.t(1);
    T=length(ytx)*dt;
    df=1e-3/T;
    f=[0:length(ytx)-1]*df;
    hold on
    plot(f,20*log10(fftshift(abs(fft(ytx)))),':b')
    
    ytx  = downsample(ytx, data.filters(channel, 1).Decimation);
    dt1=dt*data.filters(channel, 1).Decimation;
    T1=length(ytx)*dt1;
    df1=1e-3/T1;
    f1=[0:length(ytx)-1]*df1;
    % plot(f1,20*log10(abs(fft(ytx))),'r')  % original
    plot(f1,20*log10(fftshift(abs(fft(ytx)))),'r')  % original  
    
    
    % PC filter and decimation (optional)
    if strcmp(data.version, 'V4')
        ytx = conv(ytx,       data.filters(channel, 2).FilterData);
        T1a=length(ytx)*dt1;
        df1a=1e-3/T1a;
        f1a=[0:length(ytx)-1]*df1a;
        % % plot(f1a,20*log10(abs(fft(ytx))),':r')
        plot(f1a,20*log10(fftshift(abs(fft(ytx)))),':r')        
        ytx = downsample(ytx, data.filters(channel, 2).Decimation);
    end
    
    
    dt2=dt1*data.filters(channel, 2).Decimation;
    T2=length(ytx)*dt2;
    df2=1e-3/T2;
    f2=[0:length(ytx)-1]*df2;
    hold on
    % % plot(f2,(20*log10(abs(fft(ytx)))),'g')
    plot(f2,(20*log10(fftshift(abs(fft(ytx))))),'g')
%     axis([10 50 -20 60])
    title(data.echodata(channel,1).channelID)
end






