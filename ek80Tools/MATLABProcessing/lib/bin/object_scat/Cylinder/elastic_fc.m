% scattering by an solid elastic cylinder
% function     [outx, outy]=elastic_fc(proc_flag,scale,out_flag,para)
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
%  para(5) = hc   :  c2L/c1
%  para(6) = hs   :  c2T/c1
%  para(7) = theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 


function     [outx, outy]=elastic_fc(proc_flag,scale,out_flag,para)
DEG2RAD=pi/180;
ns=para(1);
x0=para(2);
xe=para(3);
g=para(4);		
hc=para(5);
hs=para(6);


if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(x0,xe,ns);
   else
     ka1=logspace(x0,xe,ns);
   end
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   m=length(ka1);
   Nmax=round(max(ka1))+10;
   theta=para(7)*DEG2RAD;
   n=0:Nmax;
   em=2*[0.5 ones(1,Nmax)];
   pn=cos(theta'*n);
   Jn1=cylbeslj(n,ka1);
   Yn1=cylbesly(n,ka1);
   dJn1=cylbesldj(n,ka1,1,Jn1);
   dYn1=cylbesldy(n,ka1,1,Yn1);
   Jn2L=cylbeslj(n,ka2L);
   dJn2L=cylbesldj(n,ka2L,1,Jn2L);
   Jn2s=cylbeslj(n,ka2s);
   dJn2s=cylbesldj(n,ka2s,1,Jn2s);
   for j=1:m
      nn=n.*n;
      tan1=-ka2L(j)*dJn2L(j,:)./Jn2L(j,:);
      tan2=-ka2s(j)*dJn2s(j,:)./Jn2s(j,:);
      tan3=-ka1(j)*dJn1(j,:)./Jn1(j,:);
      tan_beta=-ka1(j)*dYn1(j,:)./Yn1(j,:);
      tan_del=-Jn1(j,:)./Yn1(j,:);
      d1=tan1+1;
      d2=nn-ka2s(j)^2/2+tan2;
      term1a=tan1./d1;
      term1b=nn./d2;
      term2a=(nn-ka2s(j)^2/2+tan1)./d1;
      term2b=nn.*(tan2+1)./d2;
      td=-0.5*ka2s(j)^2*(term1a-term1b)./(term2a-term2b);
      tan_phi=-td/g;
      tan_eta=tan_del.*(tan_phi+tan3)./(tan_phi+tan_beta);
      cos_eta=1./sqrt(1+tan_eta.*tan_eta);
      sin_eta=tan_eta.*cos_eta;
      bn=-sin_eta.*(cos_eta+i*sin_eta);
      s=em.*pn.*bn;
      f(j)=sum(s);
   end
   outx=ka1;
else
   ka1=para(7);
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   theta=linspace(x0,xe,ns)*DEG2RAD;
   m=length(theta);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   pn=cos(theta'*n);
   em=2*[0.5 ones(1,Nmax)];
   Jn1=cylbeslj(n,ka1);
   Yn1=cylbesly(n,ka1);
   dJn1=cylbesldj(n,ka1,1,Jn1);
   dYn1=cylbesldy(n,ka1,1,Yn1);
   Jn2L=cylbeslj(n,ka2L);
   dJn2L=cylbesldj(n,ka2L,1,Jn2L);
   Jn2s=cylbeslj(n,ka2s);
   dJn2s=cylbesldj(n,ka2s,1,Jn2s);
   nn=n.*n;
   tan1=-ka2L*dJn2L./Jn2L;
   tan2=-ka2s*dJn2s./Jn2s;
   tan3=-ka1*dJn1./Jn1;
   tan_beta=-ka1*dYn1./Yn1;
   tan_del=-Jn1./Yn1;
   d1=tan1+1;
   d2=nn-ka2s^2/2+tan2;
   term1a=tan1./d1;
   term1b=nn./d2;
   term2a=(nn-ka2s^2/2+tan1)./d1;
   term2b=nn.*(tan2+1)./d2;
   td=-0.5*ka2s^2*(term1a-term1b)./(term2a-term2b);
   tan_phi=-td/g;
   tan_eta=tan_del.*(tan_phi+tan3)./(tan_phi+tan_beta);
   cos_eta=1./sqrt(1+tan_eta.*tan_eta);
   sin_eta=tan_eta.*cos_eta;
   bn=-sin_eta.*(cos_eta+i*sin_eta);
   for j=1:m
      s=em.*pn(j,:).*bn;
      f(j)=sum(s);
   end
   F=(1-g*hc*hc)/(2*g*hc*hc)-(1-g)/(g+1)*cos(theta);
   outx=theta;
end
ka=ka1;
if ( out_flag == 1 )			% modular of form function
   outy=abs(2*f./sqrt(pi*ka));		
elseif ( out_flag == 2 )		% complex form function
   outy=2*exp(i*pi/4)*f./sqrt(pi*ka);
elseif ( out_flag == 3 )		% modular of scattering amplitude
   outy=abs(f)/pi;
elseif ( out_flag == 4 )		% complex scattering amplitude
   outy=-i*f;
else
   disp(' out_flag must be within 1-4, try again !');
end
