function [gain] = getGain(data,channel)
fnom = str2num(string(data.config.transceivers(channel).channels.transducer.Frequency));
fc = (data.param(channel,1).FrequencyStart+data.param(channel, 1).FrequencyEnd)/2;

if ~isempty(data.calibration(channel).Gain)
    
    if data.param(channel,1).PulseForm == 1
        gain = interp1(data.calibration(channel).Frequency,data.calibration(channel).Gain,fnom);
        gain = gain + 20*log10(fc/fnom);
    else
        gain = data.calibration(channel).Gain;
    end
    
else
    if ~isa(data.config.transceivers(channel).channels.PulseDurationFM,'double')
        pd = strsplit(data.config.transceivers(channel).channels.PulseDurationFM,';');
        for i=1:length(pd)
            pd_table(i) = round(str2num(cell2mat(pd(i))),4);
        end
    else
        pd_table = data.config.transceivers(channel).channels.PulseDurationFM;
    end
    if ~isa(data.config.transceivers(channel).channels.transducer.Gain,'double')
        g = strsplit(data.config.transceivers(channel).channels.transducer.Gain,';');
        for i=1:length(g)
            g_table(i) = round(str2num(cell2mat(g(i))),4);
        end
    else
        g_table = data.config.transceivers(channel).channels.transducer.Gain;
    end
    gain = g_table(find(pd_table==round(data.param(channel,1).PulseDuration,4)));
    gain = gain + 20*log10(fc/fnom);
    
end
end

