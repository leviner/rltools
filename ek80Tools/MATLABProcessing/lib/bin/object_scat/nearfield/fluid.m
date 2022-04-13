% scattering by a fluid sphere due to a plane wave

clear

DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
scale=1;		% 1: linear spacing; 2: log spacing
scat_ang=180;		% scattering angle in degree
%g=0.0012;h=0.22;
g=1.04;h=1.04;


if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(0.1,30,300);
%     ka1=logspace(-1,0.5,800);
   else
     ka1=logspace(-3,1,400);
   end
   ka2=ka1/h;
   m=length(ka1);
   Nmax=max(ka1)+10;
   theta=scat_ang*DEG2RAD;
   x=cos(theta);
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2=sphbeslj(n,ka2);
   yn2=sphbesly(n,ka2);
   djn2=sphbesldj(n,ka2,1,jn2);
   term1=djn2.*yn1./(jn2.*djn1)-g*h*dyn1./djn1;
   term2=djn2.*jn1./(jn2.*djn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   for j=1:m
      s=2*nl.*pn.*bn(j,:)/ka1(j);
      f(j)=sum(s);
   end
   if (scale == 1)
%     plot(ka1,20*log10(abs(f)*0.0043/2),'r')
      a=0.0043;
      k=ka1/a;
      L=10;
      rhod=1000;
      freq=ka1/(2*pi*a)*1500*1e-3;
      plot(freq,180*a*real(f)*L*rhod./k)
   else
     loglog(freq,abs(f),'y')
   end
   xlabel('freq. (kHz)')
   ylabel('phase diff (deg)')
else
   ka1=36;
   ka2=ka1/h;
   theta=[0:2:360]*DEG2RAD;
   x=cos(theta);
   m=length(x);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2=sphbeslj(n,ka2);
   yn2=sphbesly(n,ka2);
   djn2=sphbesldj(n,ka2,1,jn2);
   term1=djn2.*yn1./(jn2.*djn1)-g*h*dyn1./djn1;
   term2=djn2.*jn1./(jn2.*djn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   for j=1:m
      s=(nl.*pn(j,:)).*bn/ka1;
      f(j)=abs(sum(s));
      fc(j)=sum(s);
   end
%   polar(theta,f/max(abs(f)))
   plot(theta*180/pi,20*log10(abs(fc)*0.0043/2))
   xlabel('Scat. Angle (deg)')
   ylabel('TS (dB)')
end
subplot(211)
plot(ka1,real(f))
ylabel('Real Part of Form function')
subplot(212)
plot(ka1,imag(f))
ylabel('Imaginary Part of Form function')
xlabel('ka')
