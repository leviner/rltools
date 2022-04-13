function     [outx, outy]=fluid_fluid_shell(proc_flag,scale,out_flag,para)

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
   x=cos(theta);
   n=round(0:Nmax);
   pn=Pn(n,x);
   nl=(2*n+1);
   
   jn1b=sphbeslj(n,k1b);
   djn1b=sphbesldj(n,k1b,1,jn1b);
   jn2b=sphbeslj(n,k2b);
   djn2b=sphbesldj(n,k2b,1,jn2b);
   jn3a=sphbeslj(n,k3a);
   djn3a=sphbesldj(n,k3a,1,jn3a);
   jn2a=sphbeslj(n,k2a);
   djn2a=sphbesldj(n,k2a,1,jn2a);
  
   yn1b=sphbesly(n,k1b);
   dyn1b=sphbesldy(n,k1b,1,yn1b);
   yn2b=sphbesly(n,k2b);
   dyn2b=sphbesldj(n,k2b,1,yn2b);
   yn2a=sphbesly(n,k2a);
   dyn2a=sphbesldy(n,k2a,1,yn2a);
   for jj=1:m
       
       d11=jn1b(jj,:)+i*yn1b(jj,:);
       d12=0*ones(size(d11));
       d13=-jn2b(jj,:);
       d14=-yn2b(jj,:);
       d21=djn1b(jj,:)+i*dyn1b(jj,:);
       d22=0*ones(size(d11));
       d23=-rho12*c12*djn2b(jj,:);
       d24=-rho12*c12*dyn2b(jj,:);
       d31=0*ones(size(d11));
       d32=-jn3a(jj,:);
       d33=jn2a(jj,:);
       d34=yn2a(jj,:);
       d41=0*ones(size(d11));
       d42=-rho23*c23*djn3a(jj,:);
       d43=djn2a(jj,:);
       d44=dyn2a(jj,:);
       
       a11=jn1b(jj,:);
       a21=djn1b(jj,:);
     
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
   outy=abs(2*f./k1b);		
elseif ( out_flag == 2 )		% complex form function
   outy=-i*2*f./k1b;
elseif ( out_flag == 3 )		% modular of scattering amplitude
   outy=abs(f)./k1b;
elseif ( out_flag == 4 )		% complex scattering amplitude
   outy=-i*f./k1b;
else
   disp(' out_flag must be within 1-4, try again !');
end
