% scattering by a rigid/fixed or soft  cylinder
% function     [outx, outy]=rgd_sft_fc(proc_flag,obj_type,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% obj_type  =	 1: rigid/fixed;    2: soft
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of form function
%		 2: complex form function
%		 3: modular of scattering amplitude normalized by L
%		 4: complex scattering amplitude normalized by L
% para ()   :   parameter array
%  para(1) = ns   :  number of data points
%  para(2) = x0   :  starting value for varible 1
%  para(3) = xe   :  ending value for varible 1
%  para(4) = theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 

function     [outx, outy]=rgd_sft_fc(proc_flag,obj_type,scale,out_flag,para)

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
   n=0:Nmax;
   pn=cos(theta*n);
   em=2*[0.5 ones(1,Nmax)];
   if (obj_type == 1)
      dJn=cylbesldj(n,ka,0);
      dYn=cylbesldy(n,ka,0);
      for j=1:m
        s=(em.*pn).*(dJn(j,:)./(dJn(j,:)+i*dYn(j,:)));
        f(j)=sum(s);
      end
   else
      Jn=cylbeslj(n,ka);
      Yn=cylbesly(n,ka);
      for j=1:m
         s=(em.*pn).*(Jn(j,:)./(Jn(j,:)+i*Yn(j,:)));
%disp(conj(s'));pause
         f(j)=sum(s);
      end
   end
   outx=ka;
else
   ka=para(4);
   theta=linspace(x0,xe,ns)*DEG2RAD;
   m=length(theta);
   Nmax=max(ka)+10;
   n=0:Nmax;
   pn=cos(theta'*n);
   em=2*[0.5 ones(1,Nmax)];
   if (obj_type == 1)
      dJn=cylbesldj(n,ka,0);
      dYn=cylbesldy(n,ka,0);
      bn=dJn./(dJn+i*dYn);
      for j=1:m
        s=(em.*pn(j,:)).*bn;
        f(j)=sum(s);
      end
   else
      Jn=cylbeslj(n,ka);
      Yn=cylbesly(n,ka);
      bn=Jn./(Jn+i*Yn);
      for j=1:m
        s=(em.*pn(j,:)).*bn;
        f(j)=sum(s);
      end
   end      
   F=-0.5+cos(theta);
   outx=theta;
end
if ( out_flag == 1 )			% modular of form function
   outy=abs(2*f./sqrt(pi*ka));		
elseif ( out_flag == 2 )		% complex form function
   outy=2*exp(i*pi/4)*f./sqrt(pi*ka);
   outy=2*f./sqrt(pi*ka);
elseif ( out_flag == 3 )		% modular of scattering amplitude
   outy=abs(f)/pi;
elseif ( out_flag == 4 )		% complex scattering amplitude
   outy=-i*f;
else
   disp(' out_flag must be within 1-4, try again !');
end
