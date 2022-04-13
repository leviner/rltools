function     [outx, outy]=fluid_fluid_shell_cyl(proc_flag,scale,out_flag,para)

% scattering by a fluid-filled fluid shell
% function     [outx, outy]=shell_fs(proc_flag,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of form function
%		 2: complex form function
%		 3: modular of scattering amplitude
%		 4: complex scattering amplitude
% para ()   :   parameter array
%  para(1) = ns   :  number of data points
%  para(2) = x0   :  starting value for varible 1
%  para(3) = xe   :  ending value for varible 1
%  para(4) = rho1 :  rho1
%  para(5) = rho2 :  rho2
%  para(6) = rho3 :  rho3
%  para(7) = c1   :  c1
%  para(8) = c2   :  c2
%  para(9) = c3   :  c3
%  para(10) = a    :  inner radius
%  para(11) = b    :  outer radius (b>a)
%  para(12)= theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 

DEG2RAD=pi/180;
rho12=para(4)./para(5);
rho23=para(5)./para(6);		
c12=para(7)/para(8);
c23=para(8)/para(9);
r=para(11)/para(10);
a=para(10);
b=para(11);
ns=para(1);
x0=para(2);
xe=para(3);

if ( proc_flag == 1)
   if (scale == 1)
     k1b=linspace(x0,xe,ns);
   else
     k1b=logspace(x0,xe,ns);
   end

   k2a=k1b.*c12.*a/b;
   k2b=k1b.*c12;
   k3a=k1b.*c12.*c23.*a/b;
   k3b=k1b.*c12.*c23;
  
   m=length(k1b);
   Nmax=max(k1b)+10;
   theta=para(12)*DEG2RAD;
   n=round(0:Nmax);
   pn=cos(n.*theta);
   nl=neumann_N(n);   

   Jn1b=besselj(n,k1b');
   dJn1b=besselJd(n,k1b');
   Jn2b=besselj(n,k2b');
   dJn2b=besselJd(n,k2b');
   Jn3a=besselj(n,k3a');
   dJn3a=besselJd(n,k3a');
   Jn2a=besselj(n,k2a');
   dJn2a=besselJd(n,k2a');
  
   Yn1b=bessely(n,k1b');
   dYn1b=besselYd(n,k1b');
   Yn2b=bessely(n,k2b');
   dYn2b=besselYd(n,k2b');
   Yn2a=bessely(n,k2a');
   dYn2a=besselYd(n,k2a');
   
   for jj=1:m
       
       d11=Jn1b(jj,:)+i*Yn1b(jj,:);
       d12=0*ones(size(d11));
       d13=-Jn2b(jj,:);
       d14=-Yn2b(jj,:);
       d21=dJn1b(jj,:)+i*dYn1b(jj,:);
       d22=0*ones(size(d11));
       d23=-rho12*c12*dJn2b(jj,:);
       d24=-rho12*c12*dYn2b(jj,:);
       d31=0*ones(size(d11));
       d32=-Jn3a(jj,:);
       d33=Jn2a(jj,:);
       d34=Yn2a(jj,:);
       d41=0*ones(size(d11));
       d42=-rho23*c23*dJn3a(jj,:);
       d43=dJn2a(jj,:);
       d44=dYn2a(jj,:);
       
       a11=-Jn1b(jj,:);
       a21=-dJn1b(jj,:);
       
       for pp=1:length(n)
         
          D=[d12(pp) d13(pp) d14(pp)
             d22(pp) d23(pp) d24(pp)
             d32(pp) d33(pp) d34(pp)
             d42(pp) d43(pp) d44(pp)];
         
          A=[a11(pp) a21(pp) 0 0];
          G=[d11(pp) d21(pp) d31(pp) d41(pp)];
          An=[A(:) D];
          Gn=[G(:) D];

          bn(pp)=det(An)/det(Gn);
      end
      s=nl.*pn.*bn;
      f(jj)=sum(s);
   end
   
   outx=k1b;
else
    fprintf('NOT YET PROGRAMMED\n')
    return
end

if ( out_flag == 1 )			% modular of form function
   outy=abs(2*exp(-i*pi/4).*f./sqrt(pi*k1b));		
elseif ( out_flag == 2 )		% complex form function
   outy=2*exp(-i*pi/4).*f./sqrt(pi*k1b);
elseif ( out_flag == 3 )		% modular of form function
   outy=abs(f);
elseif ( out_flag == 4 )		% complex form function
   outy=f;
else
   disp(' out_flag must be within 1-4, try again !');
end