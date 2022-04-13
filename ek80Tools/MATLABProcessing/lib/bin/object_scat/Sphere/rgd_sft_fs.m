% scattering by rigid/fixed or soft sphere
% function     [outx, outy]=rgd_sft_fs(proc_flag,obj_type,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% obj_type  =	 1: rigid/fixed;    2: soft
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of form function
%		 2: complex form function
%		 3: modular of scattering amplitude normalized by (a*sqrt(pi))
%		 4: complex scattering amplitude normalized by (a*sqrt(pi))
% para ()   :   parameter array
%  para(1) = ns   :  number of data points
%  para(2) = x0   :  starting value for varible 1
%  para(3) = xe   :  ending value for varible 1
%  para(4) = theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 

function     [outx, outy]=rgd_sft_fs(proc_flag,obj_type,scale,out_flag,para)

DEG2RAD=pi/180;
ns=para(1);
x0=para(2);
xe=para(3);

if ( proc_flag == 1)
   if (scale == 1)
     ka=linspace(x0,xe,ns);
   else
     ka=logspace(x0,xe,ns);
   end
   m=length(ka);
   Nmax=max(ka)+10;
   theta=para(4)*DEG2RAD;
   x=cos(theta);
   n=0:Nmax;
   pn=Pn(n,x);
   nl=-(2*n+1);
   if (obj_type == 1)
      djn=sphbesldj(n,ka,0);
      dyn=sphbesldy(n,ka,0);
      for j=1:m
        s=(nl.*pn).*(djn(j,:)./(djn(j,:)+i*dyn(j,:)));
        f(j)=sum(s);
      end
   else
      jn=sphbeslj(n,ka);
      yn=sphbesly(n,ka);
      for j=1:m
         s=(nl.*pn).*(jn(j,:)./(jn(j,:)+i*yn(j,:)));
         f(j)=sum(s);
      end      
   end
   outx=ka;
else
   ka=para(4);
   theta=linspace(x0,xe,ns)*DEG2RAD;
   x=cos(theta);
   m=length(x);
   Nmax=max(ka)+10;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=-(2*n+1);
   if (obj_type == 1)
      djn=sphbesldj(n,ka,0);
      dyn=sphbesldy(n,ka,0);
      bn=djn./(djn+i*dyn);
      for j=1:m
        s=(nl.*pn(j,:)).*bn;
        f(j)=sum(s);
      end
   else
      jn=sphbeslj(n,ka);
      yn=sphbesly(n,ka);
      bn=jn./(jn+i*yn);
      for j=1:m
        s=(nl.*pn(j,:)).*bn;
        f(j)=sum(s);
      end
   end    
   outx=theta;  
end
if ( out_flag == 1 )			% modular of form function
   outy=abs(2*f./ka);		
elseif ( out_flag == 2 )		% complex form function
   outy=-i*2*f./ka;
elseif ( out_flag == 3 )		% modular of scattering amplitude
   outy=abs(f)./ka;
elseif ( out_flag == 4 )		% complex scattering amplitude
   outy=-i*f./ka;
else
   disp(' out_flag must be within 1-4, try again !');
end
