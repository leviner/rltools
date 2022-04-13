% scattering pressure field by a fluid sphere due to a plane wave 
% function	[ka,kr,ang,ps]=fluid_sw(proc_flag,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of the normalized scattering pressure field by 1/R
%		 			  2: complex form of the normalized scattering pressure field by 1/R
%		 			  3: modular of the normalized total pressure field by 1/R
%		 			  4: complex form of the normalized total pressure field by 1/R
% para ()   :   parameter array
%  para(1) = ns   :  number of data points
%  para(2) = x0   :  starting value for varible 1
%  para(3) = xe   :  ending value for varible 1
%  para(4) = g    :  rho2/rho1
%  para(5) = h    :  c2/c1
%  para(6) = theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 
%  para(7) = ra  	:  r/a  - r: distance from the receiver to the scatterer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Dezhang Chu,			7-14-99
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function	[var,ps]=fluid_pw(proc_flag,scale,out_flag,para)

DEG2RAD=pi/180;
ns=para(1);
x0=para(2);
xe=para(3);
g=para(4);		
h=para(5);
ra=para(7);

if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(x0,xe,ns);
   else
     ka1=logspace(log10(x0),log10(xe),ns);
   end
   ka2=ka1/h;
   m=length(ka1);
   Nmax=max(ka1)+10;
   theta=para(6)*DEG2RAD;
   x=cos(theta);
   kr=ra*ka1;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   hn_r=sphhn(n,kr);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2=sphbeslj(n,ka2);
   yn2=sphbesly(n,ka2);
   djn2=sphbesldj(n,ka2,1,jn2);
   term1=djn2.*yn1./(jn2.*djn1)-g*h*dyn1./djn1;
   term2=djn2.*jn1./(jn2.*djn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   sign_term=(i).^n;
   for j=1:m
      s=-nl.*pn.*bn(j,:).*hn_r(j,:).*sign_term;
      pc(j)=sum(s);
   end
   var=ka1;
else
   ka1=para(6);
   ka2=ka1/h;
   ang=linspace(x0,xe,ns);
   theta=ang*DEG2RAD;
   x=cos(theta);
   kr=ra*ka1;
   m=length(x);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   hn_r=sphhn(n,kr);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2=sphbeslj(n,ka2);
   yn2=sphbesly(n,ka2);
   djn2=sphbesldj(n,ka2,1,jn2);
   term1=djn2.*yn1./(jn2.*djn1)-g*h*dyn1./djn1;
   term2=djn2.*jn1./(jn2.*djn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   sign_term=(i).^n;
   rest=-nl.*bn.*hn_r.*sign_term;
   for j=1:m
      s=(pn(j,:)).*rest;
      fc(j)=sum(s);
   end
   pc=fc;
   var=ang;
end

pinc=exp(i*kr*cos(theta));
switch out_flag
  case 1
    ps=abs(pc);
  case 2
    ps=pc;
  case 3
    ps=abs(pinc+pc);
  case 4
    ps=pinc+pc;
end
    
