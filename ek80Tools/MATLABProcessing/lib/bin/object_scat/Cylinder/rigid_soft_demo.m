% scattering by a rigid/fixed or soft  cylinder

clear

DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
obj_type=1;		% 1: rigid/fixed;    2: soft
scale=1;		% 1: linear spacing; 2: log spacing

if ( proc_flag == 1)
   if (scale == 1)
     ka=linspace(0.2,4.5,100);
   else
     ka=logspace(-0.2,1,100);
   end
   m=length(ka);
   Nmax=round(max(ka))+10;
   theta=180*DEG2RAD;
   n=0:Nmax;
   pn=cos(theta*n);
   em=2*[0.5 ones(1,Nmax)];
   if (obj_type == 1)
      dJn=cylbesldj(n,ka,0);
      dYn=cylbesldy(n,ka,0);
      for j=1:m
        s=(em.*pn).*(dJn(j,:)./(dJn(j,:)+i*dYn(j,:)))./sqrt(pi*ka(j));
        f(j)=2*abs(sum(s));
      end
      plot(ka,f)
   else
      Jn=cylbeslj(n,ka);
      Yn=cylbesly(n,ka);
      for j=1:m
         s=(em.*pn).*(Jn(j,:)./(Jn(j,:)+i*Yn(j,:)))./sqrt(pi*ka(j));
         f(j)=2*abs(sum(s));
      end
      semilogx(ka,f)
   end
else
   ka=0.1;
   theta=[0:359]*DEG2RAD;
   m=length(theta);
   Nmax=max(ka)+10;
   n=0:Nmax;
   pn=cos(theta'*n);
   em=2*[0.5 ones(1,Nmax)];
   if (obj_type == 1)
      dJn=cylbesldj(n,ka,0);
      dYn=cylbesldy(n,ka,0);
      bn=dJn./(dJn+i*dYn)/sqrt(pi*ka);
      for j=1:m
        s=(em.*pn(j,:)).*bn;
        f(j)=2*abs(sum(s));
      end
   else
      Jn=cylbeslj(n,ka);
      Yn=cylbesly(n,ka);
      bn=Jn./(Jn+i*Yn)/sqrt(pi*ka);
      for j=1:m
        s=(em.*pn(j,:)).*bn;
        f(j)=2*abs(sum(s));
      end
   end      
   F=-0.5+cos(theta);
   polar(theta,abs(F)/max(abs(F)),'r+')
   hold on
   polar(theta,f/max(abs(f)))
   hold off
end
   
