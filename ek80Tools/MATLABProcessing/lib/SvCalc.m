function [sv] = SvCalc(channel, ComplexVoltage, para, w, range, ping)
TVG = real(abs(20*log10(range{channel}) + 2*para{channel}.alpha*range{channel})); % b
yc = ComplexVoltage(:,ping);
prx = (abs(yc)/(sqrt(2))).^2 * ((para{channel}.zer+para{channel}.zet)/para{channel}.zer).^2*(1/para{channel}.zet);
sv = 10*log10(prx)  + TVG - 10*log10(para{channel}.ptx*(w.c/para{channel}.fc)^2 * ...
    w.c/(32*pi^2)) - 2*para{channel}.G - ...
    10*log10(para{channel}.taueff) - para{channel}.psinom + 20*log10(para{channel}.fnom/para{channel}.fc);
sv(find(sv < -999)) = 0; 
