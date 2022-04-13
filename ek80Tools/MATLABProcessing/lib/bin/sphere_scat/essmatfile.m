clear
load c:\BIOMAPER2\hti2mat\data\EN331\ESSfiles\essEN331
II=49000:52560;
YD=ess(II,1);dd = ess(II,2);
lat = ess(II,21);  lon = ess(II,22);
T = ess(II,4); S= ess(II,6);

depth = sw_dpth(dd,lat);
dens = sw_dens(S,T,dd);
svel = sw_svel(S,T,dd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ess parameters.
%lat = ess(:,21);  lon = ess(:,22);  yrday = ess(:,1);
%temp = ess(:,4); salin = ess(:,6); flour = ess(:,14);
%lite = ess(:,19); theta = ess(:,5); sigma = ess(:,7);
%tempco = ess(:,20); press = ess(:,2); ptran = ess(:,15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ess;
subplot(3,1,1)
plot(YD,-dd)
subplot(3,1,2)
plot(350./svel,-dd)
axis([.23,0.24,-250,0])
subplot(3,1,3)
%plot(dens,-dd)
hold on;
plot(0.0012*dens,-dd)
hold off

