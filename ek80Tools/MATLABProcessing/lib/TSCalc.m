function [sp, phialong, phiathw] = TSCalc(channel, data, para, range, ping)
%if strmatch(dlg.TS,'Yes')
TVG = real(40*log10(range{channel}) + 2*para{channel}.alpha*range{channel});
dens = gsw_rho(data.environ.Salinity,data.environ.Temperature,data.environ.Depth);
c = gsw_sound_speed(data.environ.Salinity,data.environ.Temperature,dens*9.81*data.environ.Depth*1e-4);

if data.param(channel,ping).PulseForm == 1
    cv = data.echodata(channel,ping).compressed;
elseif data.param(channel,ping).PulseForm == 0
    cv = data.echodata(channel,ping).complexsamples;
end

% single beam processing
if size(cv,2) == 1
    yc = nanmean(cv,2);
    prx = (abs(yc)/(sqrt(2))).^2 * ((para{channel}.zer+para{channel}.zet)/para{channel}.zer).^2*(1/para{channel}.zet);
    phialong = NaN;
    phiathw = Nan;
    
% stanard 4-quadrant split beam processing
elseif size(cv,2) ==4
    yc_i = cv;
    yc      = sum(yc_i,2)/4;
    yfore   = sum(yc_i(:,3:4),2)/2;
    yaft    = sum(yc_i(:,1:2),2)/2;
    ystar   = (yc_i(:,1) + yc_i(:,4))/2;
    yport   = sum(yc_i(:,2:3),2)/2;
    
    phialong  = angle( yfore.*conj(yaft)) *180/pi / (para{channel}.sensalong * para{channel}.fc/para{channel}.fnom);
    phiathw  = angle( ystar.*conj(yport)) *180/pi / (para{channel}.sensathw  * para{channel}.fc/para{channel}.fnom);
    phi = sqrt(phialong.^2 + phiathw.^2);
    prx = (abs(yc)/(sqrt(2))).^2 * ((para{channel}.zer+para{channel}.zet)/para{channel}.zer).^2*(1/para{channel}.zet);

% 3-quadrant, not done yet
elseif size(data.echodata(channel,ping).compressed,2) ==3
    yc = NaN;
    phialong = NaN;
    phiathw = Nan;
end

sp = 10*log10(prx)  + TVG - ...
    10*log10(para{channel}.ptx *(c/para{channel}.fc)^2/(16*pi^2)) - ...
    2*(para{channel}.G) + 20*log10(para{channel}.fc/para{channel}.fnom);
sp(find(sp < -999)) = 0; % replace -Inf values
end