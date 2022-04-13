% check formula 10.1.45-47 of Abramowitz & Stegun (pp 440)
% 8/1/96
% theta in 10.1.47 is different from those in 10.1.45 and in 10.1.46
%  in 10.1.45-46
%           thete = 0    backscattering
%           theta = 180  forward scattering
% in 10.1.47
%           thete = 180   backscattering
%           theta =   0   forward scattering

clear
theta=pi;
x=cos(theta);
kr=1;
krho=1.5;
kR=sqrt(krho^2+kr*kr-2*krho*kr*x);

n=0:65;
pn=Pn(n,x);
if krho > kr
  jn=sphbeslj(n,kr);
  hn=sphhn(n,krho);
else
  jn=sphbeslj(n,krho);
  hn=sphhn(n,kr);
end

im=i.^n;

% plane wave
S=sum( (2*n+1).*im.*pn.*jn);
S0=exp(i*min([kr krho])*x);
disp('Plane Wave:')
disp([S0;S])

% spherical wave
Ss0=exp(i*kR)/kR;
Ss=i*sum((2*n+1).*pn.*hn.*jn);
disp('Spherical Wave:')
disp([Ss0;Ss])
