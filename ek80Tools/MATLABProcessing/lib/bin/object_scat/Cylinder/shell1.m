% scattering by a sphere with an elastic shell and fluid inside

clear

DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
scale=1;		% 1: linear spacing; 2: log spacing
rho12=0.3;		% rho12 = rho1 / rho2
rho32=0.28;             % rho23 = rho3 / rho2
r=0.8;			% b/a
ft=1-r;			% fractal ratio (a-b)/a
hs=2.7;hc=3.5;h=1.1;

if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(0.1,30,100);
   else
     ka1=logspace(-3,0.5,100);
   end
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   kb2L=ka2L*r;
   kb2s=ka2s*r;
   kb3=ka1*r/h;
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
   jn2La=sphbeslj(n,ka2L);
   djn2La=sphbesldj(n,ka2L,1,jn2La);
   yn2La=sphbesly(n,ka2L);
   dyn2La=sphbesldy(n,ka2L,1,yn2La);
   jn2sa=sphbeslj(n,ka2s);
   djn2sa=sphbesldj(n,ka2s,1,jn2sa);
   yn2sa=sphbesly(n,ka2s);
   dyn2sa=sphbesldy(n,ka2s,1,yn2sa);
   jn2Lb=sphbeslj(n,kb2L);
   djn2Lb=sphbesldj(n,kb2L,1,jn2Lb);
   yn2Lb=sphbesly(n,kb2L);
   dyn2Lb=sphbesldy(n,kb2L,1,yn2Lb);
   jn2sb=sphbeslj(n,kb2s);
   djn2sb=sphbesldj(n,kb2s,1,jn2sb);
   yn2sb=sphbesly(n,kb2s);
   dyn2sb=sphbesldy(n,kb2s,1,yn2sb);
   jn3=sphbeslj(n,kb3);
   djn3=sphbesldj(n,kb3,1,jn3);
   for j=1:m
      disp(j)
      nn=2*n.*(n+1);
      a1=-rho12*ka2s(j)*ka2s(j)*jn1(j,:);
      a2=ka1(j)*djn1(j,:);
      d11=rho12*ka2s(j)*ka2s(j)*(jn1(j,:)+i*yn1(j,:));
      d12=(nn-ka2s(j)*ka2s(j)).*jn2La(j,:)-4*ka2L(j)*djn2La(j,:);
      d13=(nn-ka2s(j)*ka2s(j)).*yn2La(j,:)-4*ka2L(j)*dyn2La(j,:);
      d14=nn.*(ka2s(j)*djn2sa(j,:)-jn2sa(j,:));
      d15=nn.*(ka2s(j)*dyn2sa(j,:)-yn2sa(j,:));
      d21=-ka1(j)*(djn1(j,:)+i*dyn1(j,:));
      d22=ka2L(j)*djn2La(j,:);
      d23=ka2L(j)*dyn2La(j,:);
      d24=n.*(n+1).*jn2sa(j,:);
      d25=n.*(n+1).*yn2sa(j,:);
      d32=2*(jn2La(j,:)-ka2L(j)*djn2La(j,:));
      d33=2*(yn2La(j,:)-ka2L(j)*dyn2La(j,:));
      d34=2*ka2s(j)*djn2sa(j,:)+(ka2s(j)*ka2s(j)-nn+2).*jn2sa(j,:);
      d35=2*ka2s(j)*dyn2sa(j,:)+(ka2s(j)*ka2s(j)-nn+2).*yn2sa(j,:);
      d42=-4*kb2L(j)*djn2Lb(j,:)+(nn-kb2s(j)*kb2s(j)).*jn2Lb(j,:);
      d43=-4*kb2L(j)*dyn2Lb(j,:)+(nn-kb2s(j)*kb2s(j)).*yn2Lb(j,:);
      d44=nn.*(kb2s(j)*djn2sb(j,:)-jn2sb(j,:));
      d45=nn.*(kb2s(j)*dyn2sb(j,:)-yn2sb(j,:));
      d46=rho32*kb2s(j)*kb2s(j)*jn3(j,:);
      d52=kb2L(j)*djn2Lb(j,:);
      d53=kb2L(j)*dyn2Lb(j,:);
      d54=n.*(n+1).*jn2sb(j,:);
      d55=n.*(n+1).*yn2sb(j,:);
      d56=-kb3(j)*djn3(j,:);
      d62=2*(jn2Lb(j,:)-kb2L(j)*djn2Lb(j,:));
      d63=2*(yn2Lb(j,:)-kb2L(j)*dyn2Lb(j,:));
      d64=2*kb2s(j)*djn2sb(j,:)+(kb2s(j)*kb2s(j)-nn+2).*jn2sb(j,:);
      d65=2*kb2s(j)*dyn2sb(j,:)+(kb2s(j)*kb2s(j)-nn+2).*yn2sb(j,:);
      for l=1:length(n)
          D=[d12(l) d13(l) -d14(l) -d15(l) 0
             -d22(l) -d23(l) d24(l) d25(l) 0
             -d32(l) -d33(l) d34(l) d35(l) 0
             d42(l) d43(l) -d44(l) -d45(l) -d46(l)
             -d52(l) -d53(l) d54(l) d55(l) d56(l)
             -d62(l) -d63(l) d64(l) d65(l) 0];
          A=[-a1(l) a2(l) 0 0 0 0];
          E=[-d11(l) d21(l) 0 0 0 0];
          Bn=[A(:) D];
          Dn=[E(:) D];
          bn(l)=det(Bn)/det(Dn);
      end
      s=nl.*pn.*bn;
      f(j)=2*abs(sum(s))/ka1(j);
   end
   if (scale == 1)
     plot(ka1,f,'r')
   else
     semilogx(ka1,f)
   end
else
   ka1=10;
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   kb2L=ka2L*r;
   kb2s=ka2s*r;
   kb3=ka1*r/h;
   theta=[0:2:359]*DEG2RAD;
   x=cos(theta);
   m=length(x);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   bn=nl;
   for j=1:m
      s=nl.*pn(j,:).*bn;
      f(j)=2*abs(sum(s))/ka1;
   end
   F=(1-g*hc*hc)/(g*hc*hc)+(1-g)/g*x;
   polar(theta,F/max(abs(F)),'r+')
   hold on
   polar(theta,f/max(abs(f)))
   hold off
end
   
