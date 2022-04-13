% scattering by a fluid cylinder

clear

DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
scale=2;		% 1: linear spacing; 2: log spacing
g=1.0357;
h=1.0279;
g=0.0012;h=0.22;
ka0=5.0;
np=300;
if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(0.1,5,np);
   else
     ka1=logspace(-4,1.3,np);
   end
   ka2=ka1/h;
   m=length(ka1);
   Nmax=ceil(max(ka1))+10;
   theta=0*180*DEG2RAD;
   n=0:Nmax;
   pn=cos(theta*n);
   em=2*[0.5 ones(1,Nmax)];
   Jn1=cylbeslj(n,ka1);
   Yn1=cylbesly(n,ka1);
   dJn1=cylbesldj(n,ka1,1,Jn1);
   dYn1=cylbesldy(n,ka1,1,Yn1);
   Jn2=cylbeslj(n,ka2);
   Yn2=cylbesly(n,ka2);
   dJn2=cylbesldj(n,ka2,1,Jn2);
   term1=dJn2.*Yn1./(Jn2.*dJn1)-g*h*dYn1./dJn1;
   term2=dJn2.*Jn1./(Jn2.*dJn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   for j=1:m
      s=2*em.*pn.*bn(j,:)/sqrt(pi*ka1(j));
      f(j)=abs(sum(s));
   end
   if (scale == 1)
     plot(ka1,20*log10(f),'.-')
   else
     loglog(ka1,f./sqrt(ka1),'.-')
   end
else
   ka1=ka0;
   ka2=ka1/h;
   theta=linspace(0,360,np)*DEG2RAD;
   m=length(theta);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   pn=cos(theta'*n);
   em=2*[0.5 ones(1,Nmax)];
   Jn1=cylbeslj(n,ka1);
   Yn1=cylbesly(n,ka1);
   dJn1=cylbesldj(n,ka1,1,Jn1);
   dYn1=cylbesldy(n,ka1,1,Yn1);
   Jn2=cylbeslj(n,ka2);
   Yn2=cylbesly(n,ka2);
   dJn2=cylbesldj(n,ka2,1,Jn2);
   term1=dJn2.*Yn1./(Jn2.*dJn1)-g*h*dYn1./dJn1;
   term2=dJn2.*Jn1./(Jn2.*dJn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   for j=1:m
      s=2*(em.*pn(j,:)).*bn/sqrt(pi*ka1);  
      f(j)=abs(sum(s));
   end
   polar(theta,f)
   hold off
end
