% scattering by an ellipsoidal weakly scatterer 
% based on note 7/20/02 (Note Book: object scattering, p25)

clear

addpath c:/bin/matlab

freq=120e3;
cw=1450;
a=14e-3;                % mm to m -- semi-major of the cross section
b=8.5e-3;               % mm to m -- semi-minor of the cross section
c=69e-3;                % mm to m -- 1/2 of length
e_ba=b/a;
e_ca=c/a;
g=1.026;
h=1.017;

k=2*pi*freq/cw;
ka=k*a;
Cb=(1-g*h*h)/(g*h*h)-(g-1)/g;
coef=ka.*ka*a*e_ba*e_ca;

th=pi/2;phi=0;        % dorsal aspect
th=pi/2;phi=pi/4;     % lateral aspect

n=50;
th=pi/2;phi=linspace(0,pi/2,n);

mux=2*sin(th)*cos(phi);
muy=2*sin(th)*sin(phi);
muz=2*cos(th);

arg=(ka/h).*sqrt(mux.*mux+e_ba*e_ba*muy.*muy+e_ca*e_ca*muz.*muz)+eps;
fs=Cb*coef.*sphbeslj(1,arg)'./arg;

TS=20*log10(abs(fs));

if 1
    
if n < 10
    disp([phi(:)*180/pi TS(:)]);
else
    plot(phi*180/pi,TS,'o-')
end
grid
xlabel('SCATTERING ANGLE (deg)','fontsize',16,'fontweight','bold');
ylabel('TARGET STRENGTH (dB)','fontsize',16,'fontweight','bold');
%title('BACKSCATTERING BY SWIMBLADDERLESS FISH (38 kHz)','fontsize',16,'fontweight','bold');
ht1=text(0.0,0.05,'<-- dorsal aspect','sc');
ht2=text(0.84,0.05,'lateral aspect-->','sc');
ht3=text(0.05,0.95,sprintf('%d kHz',freq/1000),'sc');
set([ht1 ht2],'fontsize',8,'fontweight','bold');
set(ht3,'fontweight','bold');


end


