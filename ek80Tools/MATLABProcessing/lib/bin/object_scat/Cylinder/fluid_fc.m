% scattering by a fluid cylinder
% function     [outx, outy]=fluid_fc(proc_flag,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of form function
%		 2: complex form function
%		 3: modular of scattering amplitude normalized by L
%		 4: complex scattering amplitude normalized by L
% para ()   :   parameter array
%  para(1) = ns   :  number of data points
%  para(2) = x0   :  starting value for varible 1
%  para(3) = xe   :  ending value for varible 1
%  para(4) = g    :  rho2/rho1
%  para(5) = h    :  c2/c1
%  para(6) = theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 

function     [outx, outy]=fluid_fc(proc_flag,scale,out_flag,para)

DEG2RAD=pi/180;
ns=para(1);
x0=para(2);
xe=para(3);
g=para(4);		
h=para(5);	

if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(x0,xe,ns);
   else
     ka1=logspace(x0,xe,ns);
   end
   ka2=ka1/h;
   m=length(ka1);
   Nmax=round(max(ka1))+10;
   theta=para(6)*DEG2RAD;
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
      s=-em.*pn.*bn(j,:);
      f(j)=sum(s);
   end
   outx=ka1;
else
   ka1=para(6);
   ka2=ka1/h;
   theta=linspace(x0,xe,ns)*DEG2RAD;
   m=length(theta);
   Nmax=round(max(ka1))+10;
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
      s=(em.*pn(j,:)).*bn;
      f(j)=sum(s);
   end
   F=(1-g*h*h)/(2*g*h*h)-(1-g)/(1+g)*cos(theta);
   outx=theta;
end
ka=ka1;
if ( out_flag == 1 )			% modular of form function
   outy=abs(2*f./sqrt(pi*ka));		
elseif ( out_flag == 2 )		% complex form function
   outy=2*exp(-i*pi/4)*f./sqrt(pi*ka);
elseif ( out_flag == 3 )		% modular of scattering amplitude
   outy=abs(f)/pi;
elseif ( out_flag == 4 )		% complex scattering amplitude
   outy=-i*f/pi;
else
   disp(' out_flag must be within 1-4, try again !');
end
