% scattering by an solid elastic sphere

clear

DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
scale=1;		% 1: linear spacing; 2: log spacing
g=7.8;
hs=2.08;hc=3.74;

if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(0.1,1,20);
   else
     ka1=logspace(-3,0.5,100);
   end
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   m=length(ka1);
   Nmax=max(ka1)+10;
   theta=180*DEG2RAD;
   x=cos(theta);
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2L=sphbeslj(n,ka2L);
   djn2L=sphbesldj(n,ka2L,1,jn2L);
   jn2s=sphbeslj(n,ka2s);
   djn2s=sphbesldj(n,ka2s,1,jn2s);
   for j=1:m
      nn=n.*n+n;
      tan1=-ka2L(j)*djn2L(j,:)./jn2L(j,:);
      tan2=-ka2s(j)*djn2s(j,:)./jn2s(j,:);
      tan3=-ka1(j)*djn1(j,:)./jn1(j,:);
      tan_beta=-ka1(j)*dyn1(j,:)./yn1(j,:);
      tan_del=-jn1(j,:)./yn1(j,:);
      d1=tan1+1;
      d2=nn-1-ka2s(j)^2/2+tan2;
      term1a=tan1./d1;
      term1b=nn./d2;
      term2a=(nn-ka2s(j)^2/2+2*tan1)./d1;
      term2b=nn.*(tan2+1)./d2;
      td=-0.5*ka2s(j)^2*(term1a-term1b)./(term2a-term2b);
      tan_phi=-td/g;
      tan_eta=tan_del.*(tan_phi+tan3)./(tan_phi+tan_beta);
      cos_eta=1./sqrt(1+tan_eta.*tan_eta);
      sin_eta=tan_eta.*cos_eta;
      bn=-sin_eta.*(cos_eta+i*sin_eta);
      s=nl.*pn.*bn;
      f(j)=2*sum(s)/ka1(j);
   end
   if (scale == 1)
     plot(ka1,imag(f),'g')
   else
     semilogx(ka1,f)
   end
else
   ka1=0.2;
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   theta=[0:2:359]*DEG2RAD;
   x=cos(theta);
   m=length(x);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2L=sphbeslj(n,ka2L);
   djn2L=sphbesldj(n,ka2L,1,jn2L);
   jn2s=sphbeslj(n,ka2s);
   djn2s=sphbesldj(n,ka2s,1,jn2s);
   nn=n.*n+n;
      tan1=-ka2L*djn2L./jn2L;
      tan2=-ka2s*djn2s./jn2s;
      tan3=-ka1*djn1./jn1;
      tan_beta=-ka1*dyn1./yn1;
      tan_del=-jn1./yn1;
      d1=tan1+1;
      d2=nn-1-ka2s^2/2+tan2;
      term1a=tan1./d1;
      term1b=nn./d2;
      term2a=(nn-ka2s^2/2+2*tan1)./d1;
      term2b=nn.*(tan2+1)./d2;
      td=-0.5*ka2s^2*(term1a-term1b)./(term2a-term2b);
      tan_phi=-td/g;
      tan_eta=tan_del.*(tan_phi+tan3)./(tan_phi+tan_beta);
      cos_eta=1./sqrt(1+tan_eta.*tan_eta);
      sin_eta=tan_eta.*cos_eta;
      bn=-sin_eta.*(cos_eta+i*sin_eta);
   for j=1:m
      s=nl.*pn(j,:).*bn;
      f(j)=2*abs(sum(s))/ka1;
   end
   polar(theta,f/abs(max(f)))
   hold on
   coef=ka1^2/3;
   F=coef*((1-g*hc*hc)/(g*hc*hc)+(1-g)/g*x);
   polar(theta,F/max(abs(F)),'r+')
   hold off
end
   
