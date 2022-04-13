% scattering by rigid/fixed or soft sphere

clear

DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
obj_type=1;		% 1: rigid/fixed;    2: soft
scale=1;		% 1: linear spacing; 2: log spacing

if ( proc_flag == 1)
   if (scale == 1)
     ka=linspace(0.01,10,100);
   else
     ka=logspace(-2,1,100);
   end
   m=length(ka);
   Nmax=max(ka)+10;
   theta=180*DEG2RAD;
   x=cos(theta);
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   if (obj_type == 1)
      djn=sphbesldj(n,ka,0);
      dyn=sphbesldy(n,ka,0);
      for j=1:m
        s=(nl.*pn).*(djn(j,:)./(djn(j,:)+i*dyn(j,:)))./ka(j);
        f(j)=2*abs(sum(s));
      end
      plot(ka,f,'r')
   else
      jn=sphbeslj(n,ka);
      yn=sphbesly(n,ka);
      for j=1:m
         s=(nl.*pn).*(jn(j,:)./(jn(j,:)+i*yn(j,:)))./ka(j);
         f(j)=2*abs(sum(s));
      end      
      semilogx(ka,f)
   end
else
   ka=10;
   theta=[0:359]*DEG2RAD;
   x=cos(theta);
   m=length(x);
   Nmax=max(ka)+10;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=2*(2*n+1);
   if (obj_type == 1)
      djn=sphbesldj(n,ka,0);
      dyn=sphbesldy(n,ka,0);
      bn=djn./(djn+i*dyn)/ka;
      for j=1:m
        s=(nl.*pn(j,:)).*bn;
        f(j)=abs(sum(s));
      end
      polar(theta,f)
   else
      jn=sphbeslj(n,ka);
      yn=sphbesly(n,ka);
      bn=jn./(jn+i*yn)/ka;
      for j=1:m
        s=(nl.*pn(j,:)).*bn;
        f(j)=abs(sum(s));
      end
      polar(theta,f)
   end      
end
   
