% scattering by a cylinder with an elastic shell and fluid inside
% function     [outx, outy]=shell_fc(proc_flag,scale,out_flag,para)
% proc_flag =	 1: vs. ka;         2: vs. angle
% scale	    =	 1: linear spacing; 2: log spacing
% out_flag  =    1: modular of form function
%		 2: complex form function
%		 3: modular of scattering amplitude normalized by L
%		 4: complex scattering amplitude normalized by L
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


function     [outx, outy]=shell_fc(proc_flag,scale,out_flag,para)

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
   n=0:Nmax;
   pn=cos(theta'*n);
   em=2*[0.5 ones(1,Nmax)];
   Jn1=cylbeslj(n,ka1);
   Yn1=cylbesly(n,ka1);
   dJn1=cylbesldj(n,ka1,1,Jn1);
   dYn1=cylbesldy(n,ka1,1,Yn1);
   Jn2La=cylbeslj(n,ka2L);
   dJn2La=cylbesldj(n,ka2L,1,Jn2La);
   Yn2La=cylbesly(n,ka2L);
   dYn2La=cylbesldy(n,ka2L,1,Yn2La);
   Jn2sa=cylbeslj(n,ka2s);
   dJn2sa=cylbesldj(n,ka2s,1,Jn2sa);
   Yn2sa=cylbesly(n,ka2s);
   dYn2sa=cylbesldy(n,ka2s,1,Yn2sa);
   Jn2Lb=cylbeslj(n,kb2L);
   dJn2Lb=cylbesldj(n,kb2L,1,Jn2Lb);
   Yn2Lb=cylbesly(n,kb2L);
   dYn2Lb=cylbesldy(n,kb2L,1,Yn2Lb);
   Jn2sb=cylbeslj(n,kb2s);
   dJn2sb=cylbesldj(n,kb2s,1,Jn2sb);
   Yn2sb=cylbesly(n,kb2s);
   dYn2sb=cylbesldy(n,kb2s,1,Yn2sb);
   Jn3=cylbeslj(n,kb3);
   dJn3=cylbesldj(n,kb3,1,Jn3);
   for j=1:m
      nn=2*n.*n;
      a1=rho12*ka2s(j)*ka2s(j)*Jn1(j,:);
      a2=ka1(j)*dJn1(j,:);
      d11=-rho12*ka2s(j)*ka2s(j)*(Jn1(j,:)+i*Yn1(j,:));
      d12=(nn-ka2s(j)*ka2s(j)).*Jn2La(j,:)-2*ka2L(j)*dJn2La(j,:);
      d13=(nn-ka2s(j)*ka2s(j)).*Yn2La(j,:)-2*ka2L(j)*dYn2La(j,:);
      d14=-2*n.*(ka2s(j)*dJn2sa(j,:)-Jn2sa(j,:));
      d15=-2*n.*(ka2s(j)*dYn2sa(j,:)-Yn2sa(j,:));
      d21=-ka1(j)*(dJn1(j,:)+i*dYn1(j,:));
      d22=-ka2L(j)*dJn2La(j,:);
      d23=-ka2L(j)*dYn2La(j,:);
      d24=n.*Jn2sa(j,:);
      d25=n.*Yn2sa(j,:);
      d32=-2*n.*(Jn2La(j,:)-ka2L(j)*dJn2La(j,:));
      d33=-2*n.*(Yn2La(j,:)-ka2L(j)*dYn2La(j,:));
      d34=2*ka2s(j)*dJn2sa(j,:)+(ka2s(j)*ka2s(j)-nn).*Jn2sa(j,:);
      d35=2*ka2s(j)*dYn2sa(j,:)+(ka2s(j)*ka2s(j)-nn).*Yn2sa(j,:);
      d42=-2*kb2L(j)*dJn2Lb(j,:)+(nn-kb2s(j)*kb2s(j)).*Jn2Lb(j,:);
      d43=-2*kb2L(j)*dYn2Lb(j,:)+(nn-kb2s(j)*kb2s(j)).*Yn2Lb(j,:);
      d44=-2*n.*(kb2s(j)*dJn2sb(j,:)-Jn2sb(j,:));
      d45=-2*n.*(kb2s(j)*dYn2sb(j,:)-Yn2sb(j,:));
      d46=-rho32*kb2s(j)*kb2s(j)*Jn3(j,:);
      d52=-kb2L(j)*dJn2Lb(j,:);
      d53=-kb2L(j)*dYn2Lb(j,:);
      d54=n.*Jn2sb(j,:);
      d55=n.*Yn2sb(j,:);
      d56=-kb3(j)*dJn3(j,:);
      d62=-2*n.*(Jn2Lb(j,:)-kb2L(j)*dJn2Lb(j,:));
      d63=-2*n.*(Yn2Lb(j,:)-kb2L(j)*dYn2Lb(j,:));
      d64=2*kb2s(j)*dJn2sb(j,:)+(kb2s(j)*kb2s(j)-nn).*Jn2sb(j,:);
      d65=2*kb2s(j)*dYn2sb(j,:)+(kb2s(j)*kb2s(j)-nn).*Yn2sb(j,:);
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
      s=em.*pn.*bn;
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
   m=length(theta);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   em=2*[0.5 ones(1,Nmax)];
   pn=cos(theta'*n);
   Jn1=cylbeslj(n,ka1);
   Yn1=cylbesly(n,ka1);
   dJn1=cylbesldj(n,ka1,1,Jn1);
   dYn1=cylbesldy(n,ka1,1,Yn1);
   Jn2La=cylbeslj(n,ka2L);
   dJn2La=cylbesldj(n,ka2L,1,Jn2La);
   Yn2La=cylbesly(n,ka2L);
   dYn2La=cylbesldy(n,ka2L,1,Yn2La);
   Jn2sa=cylbeslj(n,ka2s);
   dJn2sa=cylbesldj(n,ka2s,1,Jn2sa);
   Yn2sa=cylbesly(n,ka2s);
   dYn2sa=cylbesldy(n,ka2s,1,Yn2sa);
   Jn2Lb=cylbeslj(n,kb2L);
   dJn2Lb=cylbesldj(n,kb2L,1,Jn2Lb);
   Yn2Lb=cylbesly(n,kb2L);
   dYn2Lb=cylbesldy(n,kb2L,1,Yn2Lb);
   Jn2sb=cylbeslj(n,kb2s);
   dJn2sb=cylbesldj(n,kb2s,1,Jn2sb);
   Yn2sb=cylbesly(n,kb2s);
   dYn2sb=cylbesldy(n,kb2s,1,Yn2sb);
   Jn3=cylbeslj(n,kb3);
   dJn3=cylbesldj(n,kb3,1,Jn3);
      nn=2*n.*n;
      a1=rho12*ka2s*ka2s*Jn1;
      a2=ka1*dJn1;
      d11=rho12*ka2s*ka2s*(Jn1+i*Yn1);
      d12=(nn-ka2s*ka2s).*Jn2La-2*ka2L*dJn2La;
      d13=(nn-ka2s*ka2s).*Yn2La-2*ka2L*dYn2La;
      d14=-2*n.*(ka2s*dJn2sa-Jn2sa);
      d15=-2*n.*(ka2s*dYn2sa-Yn2sa);
      d21=-ka1*(dJn1+i*dYn1);
      d22=-ka2L*dJn2La;
      d23=-ka2L*dYn2La;
      d24=n.*Jn2sa;
      d25=n.*Yn2sa;
      d32=-2*n.*(Jn2La-ka2L*dJn2La);
      d33=-2*n.*(Yn2La-ka2L*dYn2La);
      d34=2*ka2s*dJn2sa+(ka2s*ka2s-nn).*Jn2sa;
      d35=2*ka2s*dYn2sa+(ka2s*ka2s-nn).*Yn2sa;
      d42=-2*kb2L*dJn2Lb+(nn-kb2s*kb2s).*Jn2Lb;
      d43=-2*kb2L*dYn2Lb+(nn-kb2s*kb2s).*Yn2Lb;
      d44=-2*n.*(kb2s*dJn2sb-Jn2sb);
      d45=-2*n.*(kb2s*dYn2sb-Yn2sb);
      d46=-rho32*kb2s*kb2s*Jn3;
      d52=-kb2L*dJn2Lb;
      d53=-kb2L*dYn2Lb;
      d54=n.*Jn2sb;
   d55=n.*Yn2sb;
   d56=-kb3*dJn3;
   d62=-2*(Jn2Lb-kb2L*dJn2Lb);
   d63=-2*(Yn2Lb-kb2L*dYn2Lb);
   d64=2*kb2s*dJn2sb+(kb2s*kb2s-nn).*Jn2sb;
   d65=2*kb2s*dYn2sb+(kb2s*kb2s-nn).*Yn2sb;
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
      s=em.*pn(j,:).*bn;
      f(j)=sum(s);
   end
   F=(rho12-hc*hc)/(2*hc*hc)-(rho12-1)/(rho12+1)*cos(theta);
   outx=theta;
end
ka=ka1;
if ( out_flag == 1 )			% modular of form function
   outy=abs(2*f./sqrt(pi*ka));		
elseif ( out_flag == 2 )		% complex form function
   outy=2*exp(i*pi/4)*f./sqrt(pi*ka);
elseif ( out_flag == 3 )		% modular of scattering amplitude
   outy=abs(f)/pi;
elseif ( out_flag == 4 )		% complex scattering amplitude
   outy=-i*f;
else
   disp(' out_flag must be within 1-4, try again !');
end
