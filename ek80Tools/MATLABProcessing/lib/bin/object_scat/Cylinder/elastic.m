% scattering by an solid elastic cylinder

clear

DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
scale=1;		% 1: linear spacing; 2: log spacing
g=7.9;
hs=2.08;hc=3.74;

if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(0.1,15,600);
   else
     ka1=logspace(-3,0.5,100);
   end
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   m=length(ka1);
   Nmax=max(ka1)+10;
   theta=180*DEG2RAD;
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
      f(j)=2*abs(sum(s))/sqrt(pi*ka1(j));
   end
   if (scale == 1)
     plot(ka1,f,'g')
   else
     semilogx(ka1,f)
   end
else
   ka1=0.1;
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   theta=[0:2:359]*DEG2RAD;
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
      f(j)=2*abs(sum(s))/sqrt(pi*ka1);
   end
   F=(1-g*hc*hc)/(2*g*hc*hc)-(1-g)/(g+1)*cos(theta);
   polar(theta,abs(F)/max(abs(F)),'r+')
   hold on
   polar(theta,f/max(abs(f)))
   hold off
end
   
