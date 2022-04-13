
function [sp, phialong, phiathw] = TSCalc(channel, data, para, w, range, ping)
%if strmatch(dlg.TS,'Yes')
TVG = real(40*log10(range{channel}) + 2*para{channel}.alpha*range{channel});

% single beam processing
if size(data.echodata(channel,ping).compressed,2) == 1
    disp('Single beam TS')
    yc = nanmean(data.echodata(channel,ping).compressed,2);
    prx = (abs(yc)/(sqrt(2))).^2 * ((para{channel}.zer+para{channel}.zet)/para{channel}.zer).^2*(1/para{channel}.zet);
    sp = 10*log10(prx)  + TVG - ...
        10*log10(para{channel}.ptx *(w.c/para{channel}.fc)^2/(16*pi^2)) - ...
        2*(para{channel}.G) + 20*log10(para{channel}.fc/para{channel}.fnom);
    sp(find(sp < -999)) = 0; % replace -Inf values
    phialong = NaN;
    phiathw = Nan;
    
elseif size(data.echodata(channel,ping).compressed,2) ==4
    disp('4-sector TS')
    % stanard 4-quadrant split beam processing
    yc_i = data.echodata(channel,ping).compressed;
    yc      = sum(yc_i,2)/4;
    yfore   = sum(yc_i(:,3:4),2)/2;
    yaft    = sum(yc_i(:,1:2),2)/2;
    ystar   = (yc_i(:,1) + yc_i(:,4))/2;
    yport   = sum(yc_i(:,2:3),2)/2;
    
    phialong  = angle( yfore.*conj(yaft)) *180/pi / (para{channel}.sensalong * para{channel}.fc/para{channel}.fnom);
    phiathw  = angle( ystar.*conj(yport)) *180/pi / (para{channel}.sensathw  * para{channel}.fc/para{channel}.fnom);
    phi = sqrt(phialong.^2 + phiathw.^2);
    prx = (abs(yc)/(sqrt(2))).^2 * ((para{channel}.zer+para{channel}.zet)/para{channel}.zer).^2*(1/para{channel}.zet);
    
    sp = 10*log10(prx)  + TVG - ...
        10*log10(para{channel}.ptx *(w.c/para{channel}.fc)^2/(16*pi^2)) - ...
        2*(para{channel}.G) + 20*log10(para{channel}.fc/para{channel}.fnom);
    sp(find(sp < -999)) = 0; % replace -Inf values
    
elseif size(data.echodata(channel,ping).compressed,2) ==3
    %3-sector. not functional yet
end
end