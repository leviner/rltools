clear
% function     [outx, outy]=fluid_fs(proc_flag,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of form function
%		 2: complex form function
%		 3: modular of scattering amplitude normalized by (a*sqrt(pi))
%		 4: complex scattering amplitude normalized by (a*sqrt(pi))
% para ()   :   parameter array
%  para(1) = ns   :  number of data points
%  para(2) = x0   :  starting value for varible 1
%  para(3) = xe   :  ending value for varible 1
%  para(4) = g    :  rho2/rho1
%  para(5) = h    :  c2/c1
%  para(6) = theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 


proc_flag = 1;S=1;out_flag=3;
g0=.0012;h0=.22; %material properties for air bubble at surface
z_moc=[0,25,50,75,100,125,150,175];  %depth
g=.0012*(1+0.1*z_moc);
x0=[100,1e-3,3,g0,h0,180];
[ka0,fbs0]=fluid_fs(proc_flag,S,out_flag,x0);
RTS0=10*log10((abs(fbs0)).^2);
plot(ka0,RTS0)
load c:\WORK\SBworkshop\look_up_tables\siphon_tbl.dat
hold on
plot(siphon_tbl(:,1),siphon_tbl(:,3),'r','linewidth',[2])
hold off

hold on
for ii=1:8
   g_temp=g(ii);
   x=[150,-4.5,-0.5,g_temp,h0,180];
   [ka_temp1,fbs_temp1]=fluid_fs(proc_flag,2,out_flag,x);
   x=[18,-0.49,0.65,g_temp,h0,180];
   [ka_temp2,fbs_temp2]=fluid_fs(proc_flag,2,out_flag,x);
   ka_temp=[ka_temp1, ka_temp2];
   fbs_temp=[fbs_temp1, fbs_temp2];
   
   RTS_temp=10*log10((abs(fbs_temp)).^2);
   ka(:,ii)=ka_temp';
   fbs(:,ii)=fbs_temp';
   RTS(:,ii)=RTS_temp';
   plot(ka_temp,RTS_temp)
end
hold off
h=h0;
save siphon_EN331M5_tbl ka RTS fbs g h z_moc

%December 1999 EN331 MOC tow 5 
%net    depth     nominal depth
%5.1 175-211      175
%5.2 150-174      150 
%5.3 125-149      125
%5.4 100-123      100
%5.5  75-100      75 
%5.6   50-74      50
%5.7   25-50      25
%5.8    0-24      0
