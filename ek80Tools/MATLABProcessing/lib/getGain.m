
function [gain] = getGain(data,channel)


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
    for i=1:length(pd)
        g_table(i) = round(str2num(cell2mat(g(i))),4);
    end
else
    g_table = data.config.transceivers(channel).channels.transducer.Gain;
end
gain = g_table(find(pd_table==round(data.param(channel,1).PulseDuration,4)));

end
