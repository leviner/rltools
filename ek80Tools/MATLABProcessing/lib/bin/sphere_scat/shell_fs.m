% scattering by a sphere with an elastic shell with fluid inside
% function     [outx, outy]=shell_fs(proc_flag,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of form function
%		 2: complex form function
%		 3: modular of scattering amplitude normalized by (a*sqrt(pi))
%		 4: complex scattering amplitude normalized by (a*sqrt(pi))
% para ()   :   parameter array
%  para(1) = ns   :  number of data points
%  para(2) = x0   :  starting value for varible 1
%  para(3) = xe   :  ending value for varible 1
%  para(4) = rho12:  rho1/rho2
%  para(5) = rho32:  rho3/rho2
%  para(6) = hc   :  c2L/c1
%  para(7) = hs   :  c2T/c1
%  para(8) = h    :  c3/c1
%  para(9) = r    :  b/a
%  para(10)= theta0  for proc_flag = 1  (degree)
%          = ka0     for proc_flag = 2 


function     [outx, outy]=shell_fs(proc_flag,scale,out_flag,para)

DEG2RAD=pi/180;
ns=para(1);
x0=para(2);
xe=para(3);
rho12=para(4);
rho32=para(5);		
hc=para(6);
hs=para(7);
h=para(8);
r=para(9);

ft=1-r;			% fractal ratio (a-b)/a
tic
if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(x0,xe,ns);
   else
     ka1=logspace(x0,xe,ns);
   end
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   kb2L=ka2L*r;
   kb2s=ka2s*r;
   kb3=ka1*r/h;
   m=length(ka1);
   Nmax=max(ka1)+10;
   theta=para(10)*DEG2RAD;
   x=cos(theta);
   n=round(0:Nmax);
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
          D=[d12(l) d13(l) d14(l) d15(l) 0
             d22(l) d23(l) d24(l) d25(l) 0
             d32(l) d33(l) d34(l) d35(l) 0
             d42(l) d43(l) d44(l) d45(l) d46(l)
             d52(l) d53(l) d54(l) d55(l) d56(l)
             d62(l) d63(l) d64(l) d65(l) 0];
          A=[a1(l) a2(l) 0 0 0 0];
          E=[d11(l) d21(l) 0 0 0 0];
          Bn=[A(:) D];
          Dn=[E(:) D];
          detBn=det(Bn);
          detDn=det(Dn);
          bn(l)=detBn/detDn;
      end
      s=nl.*pn.*bn;
      f(j)=sum(s);
   end
   outx=ka1;
else
   ka1=para(10);
   ka2L=ka1/hc;
   ka2s=ka1/hs;
   kb2L=ka2L*r;
   kb2s=ka2s*r;
   kb3=ka1*r/h;
   theta=linspace(x0,xe,ns)*DEG2RAD;
   x=cos(theta);
   m=length(x);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   nl=(2*n+1);
   pn=Pn(n,x);
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
   nn=2*n.*(n+1);
   a1=-rho12*ka2s*ka2s*jn1;
   a2=ka1*djn1;
   d11=rho12*ka2s*ka2s*(jn1+i*yn1);
   d12=(nn-ka2s*ka2s).*jn2La-4*ka2L*djn2La;
   d13=(nn-ka2s*ka2s).*yn2La-4*ka2L*dyn2La;
   d14=nn.*(ka2s*djn2sa-jn2sa);
   d15=nn.*(ka2s*dyn2sa-yn2sa);
   d21=-ka1*(djn1+i*dyn1);
   d22=ka2L*djn2La;
   d23=ka2L*dyn2La;
   d24=n.*(n+1).*jn2sa;
   d25=n.*(n+1).*yn2sa;
   d32=2*(jn2La-ka2L*djn2La);
   d33=2*(yn2La-ka2L*dyn2La);
   d34=2*ka2s*djn2sa+(ka2s*ka2s-nn+2).*jn2sa;
   d35=2*ka2s*dyn2sa+(ka2s*ka2s-nn+2).*yn2sa;
   d42=-4*kb2L*djn2Lb+(nn-kb2s*kb2s).*jn2Lb;
   d43=-4*kb2L*dyn2Lb+(nn-kb2s*kb2s).*yn2Lb;
   d44=nn.*(kb2s*djn2sb-jn2sb);
   d45=nn.*(kb2s*dyn2sb-yn2sb);
   d46=rho32*kb2s*kb2s*jn3;
   d52=kb2L*djn2Lb;
   d53=kb2L*dyn2Lb;
   d54=n.*(n+1).*jn2sb;
   d55=n.*(n+1).*yn2sb;
   d56=-kb3*djn3;
   d62=2*(jn2Lb-kb2L*djn2Lb);
   d63=2*(yn2Lb-kb2L*dyn2Lb);
   d64=2*kb2s*djn2sb+(kb2s*kb2s-nn+2).*jn2sb;
   d65=2*kb2s*dyn2sb+(kb2s*kb2s-nn+2).*yn2sb;
   for l=1:length(n)
       D=[d12(l) d13(l) d14(l) d15(l) 0
          d22(l) d23(l) d24(l) d25(l) 0
          d32(l) d33(l) d34(l) d35(l) 0
          d42(l) d43(l) d44(l) d45(l) d46(l)
          d52(l) d53(l) d54(l) d55(l) d56(l)
          d62(l) d63(l) d64(l) d65(l) 0];
       A=[a1(l) a2(l) 0 0 0 0];
       E=[d11(l) d21(l) 0 0 0 0];
       Bn=[A(:) D];
       Dn=[E(:) D];
       bn(l)=det(Bn)/det(Dn);
   end
   for j=1:m
      s=nl.*pn(j,:).*bn;
      f(j)=sum(s);
   end
   F=(rho12-hc*hc)/(hc*hc)+(rho12-1)*x;
   outx=theta;
end
ka=ka1;
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
