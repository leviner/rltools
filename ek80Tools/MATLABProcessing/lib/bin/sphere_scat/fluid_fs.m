% scattering by a fluid sphere
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

function     [outx, outy]=fluid_fs(proc_flag,scale,out_flag,para)

DEG2RAD=pi/180;
ns=para(1);
x0=para(2);
xe=para(3);
g=para(4);		
h=para(5);		
if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(x0,xe,ns);
     Nmax=max(ka1)+10;
  else
     ka1=logspace(x0,xe,ns);
     if ns == 1
        Nmax=10^max(ka1)+10;
     else
        Nmax=max(ka1)+10;
	  end
   end
   ka2=ka1/h;
   m=length(ka1);
   theta=para(6)*DEG2RAD;
   x=cos(theta);
   n=0:round(Nmax);
   pn1=Pn(n,x);
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
   bn=-1./(1+i*cn);
%   term1a=djn2.*yn1-g*h*dyn1.*jn2;
%   term2a=djn2.*jn1-g*h*jn2.*djn1;
%   cn_a=term1a./term2a;
%   bn_a=-term2a./(term2a+i*term1a);
   for j=1:m
      s=nl.*pn1.*bn(j,:);
      f(j)=sum(s);
%      disp(ka1(j));plot(abs(cumsum(s)),'o-');pause
   end
%   loglog(ka1,abs(term1(:,1)),ka1,abs(term2(:,1)));pause
   outx=ka1;
else
   ka1=para(6);
   ka2=ka1/h;
   theta=linspace(x0,xe,ns)*DEG2RAD;
   x=cos(theta);
   m=length(x);
   Nmax=max(ka1)+10;
   n=0:round(Nmax);
   pn1=Pn(n,x);
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
   bn=-1./(1+i*cn);
   for j=1:m
      s=(nl.*pn1(j,:)).*bn;
      f(j)=sum(s);
   end
   F=(1-g*h*h)/(g*h*h)+(1-g)/g*x;
   outx=theta;
end
ka=ka1;
if ( out_flag == 1 )			% modular of form function
   outy=abs(2*f./ka);		
elseif ( out_flag == 2 )		% complex form function
   outy=-i*2*f./ka;
elseif ( out_flag == 3 )		% modular of scattering amplitude
   outy=abs(f)./(sqrt(pi)*ka);
elseif ( out_flag == 4 )		% complex scattering amplitude
   outy=-i*f./(sqrt(pi)*ka);
else
   disp(' out_flag must be within 1-4, try again !');
end
%plot(outx,abs(outy))
