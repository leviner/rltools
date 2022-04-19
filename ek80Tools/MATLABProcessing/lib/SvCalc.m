function [sv, yc] = SvCalc(channel, data, para, range, ping)

TVG = real(abs(20*log10(range{channel}) + 2*para{channel}.alpha*range{channel})); % b
dens = gsw_rho(data.environ.Salinity,data.environ.Temperature,data.environ.Depth);
c = gsw_sound_speed(data.environ.Salinity,data.environ.Temperature,dens*9.81*data.environ.Depth*1e-4);

if data.param(channel,ping).PulseForm == 1
    yc = mean(data.echodata(channel,ping).compressed,2);
elseif data.param(channel,ping).PulseForm == 0
    yc = mean(data.echodata(channel,ping).complexsamples,2);
end

prx = (abs(yc)/(sqrt(2))).^2 * ((para{channel}.zer+para{channel}.zet)/para{channel}.zer).^2*(1/para{channel}.zet);

sv = 10*log10(prx)  + TVG - 10*log10(para{channel}.ptx*(c/para{channel}.fc)^2 * ...
    c/(32*pi^2)) - 2*para{channel}.G - ...
    10*log10(para{channel}.taueff) - para{channel}.psinom + 20*log10(para{channel}.fnom/para{channel}.fc);
sv(find(sv < -999)) = 0;
