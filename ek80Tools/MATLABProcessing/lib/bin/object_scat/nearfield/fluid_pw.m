% scattering by a fluid sphere due to a plane wave
% pressure field

clear

n=300;
DEG2RAD=pi/180;
proc_flag=1;		% 1: vs. ka;         2: vs. angle
scale=1;		% 1: linear spacing; 2: log spacing
scat_ang=180;		% scattering angle in degree
ka0=0.1;		% fixed ka where the scat. directivity is computed
kae=10;
g=2.0;h=2.0;		% fluid 1
g=1000.;h=1000;		% rigid fixed
g=1.05;h=1.05;		% fluid 1
%g=0.0012;h=0.22;	% gas bubble
r_to_a=100;
title_str=sprintf('Fluid Sphere (Backscattering): g = %g, h = %g',g,h);

if ( proc_flag == 1)
   if (scale == 1)
     ka1=linspace(ka0,kae,n);
   else
     ka1=logspace(-4,1.0,400);
   end
   ka2=ka1/h;
   m=length(ka1);
   Nmax=max(ka1)+10;
   theta=scat_ang*DEG2RAD;
   x=cos(theta);
   kr=r_to_a*ka1;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   hn_r=sphhn(n,kr);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2=sphbeslj(n,ka2);
   yn2=sphbesly(n,ka2);
   djn2=sphbesldj(n,ka2,1,jn2);
   term1=djn2.*yn1./(jn2.*djn1)-g*h*dyn1./djn1;
   term2=djn2.*jn1./(jn2.*djn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   sign_term=i.^n;
   for j=1:m
      s=-nl.*pn.*bn(j,:).*hn_r(j,:).*sign_term;
      s1=i*2*nl.*pn.*bn(j,:)/ka1(j);
      p(j)=sum(s);
      p1(j)=exp(i*kr(j))*sum(s1)/(2*r_to_a);
%      [amp, engy_r(j)]=energy(j,s,ka1(j)); engy(:,j)=amp(:);
   end
   if scale == 2
   subplot(311)
   loglog(ka1,abs(real(p1)),ka1,abs(real(p)),'.')
   text(1,10000,'solid: Plane wave')
   text(1,200,'dotted: Spherical')
   ylabel('Re(p)')
   title(title_str)
   subplot(312)
   loglog(ka1,abs(imag(p1)),ka1,abs(imag(p)),'.')
   ylabel('Im(p)')
   subplot(313)
   loglog(ka1,abs(p1),ka1,abs(p),'.')
   ylabel('|p|')
  xlabel('ka')
   end
else
   ka1=ka0;
   ka2=ka1/h;
   theta=[0:2:360]*DEG2RAD;
   x=cos(theta);
   kr=r_to_a*ka1;
   m=length(x);
   Nmax=max(ka1)+10;
   n=0:Nmax;
   pn=Pn(n,x);
   nl=(2*n+1);
   jn1=sphbeslj(n,ka1);
   yn1=sphbesly(n,ka1);
   hn_r=sphhn(n,kr);
   djn1=sphbesldj(n,ka1,1,jn1);
   dyn1=sphbesldy(n,ka1,1,yn1);
   jn2=sphbeslj(n,ka2);
   yn2=sphbesly(n,ka2);
   djn2=sphbesldj(n,ka2,1,jn2);
   term1=djn2.*yn1./(jn2.*djn1)-g*h*dyn1./djn1;
   term2=djn2.*jn1./(jn2.*djn1)-g*h;
   cn=term1./term2;
   bn=1./(1+i*cn);
   sign_term=(-1).^n;
   rest=nl.*bn.*hn_r.*sign_term/ka1;
   for j=1:m
      s=(pn(j,:)).*rest;
      f(j)=abs(sum(s));
      fc(j)=sum(s);
   end
   h=polar(theta,f/max(abs(f)));
   set(h,'linewidth',1.5)
   xlabel(['Scat. Angle (deg), ' sprintf(' ka = %g',ka1)])
   ylabel('Normalized Scattering Pressure')
   title(title_str)
end

subplot(311)
plot(ka1,real(p1),ka1,real(p),'.r')
ylabel('Re(p_{scat})')
legend('Approx.','Exact')
subplot(312)
plot(ka1,imag(p1),ka1,imag(p),'.r')
ylabel('Im(p_{scat})')
subplot(313)
plot(ka1,abs(p1),ka1,abs(p),'.r')
ylabel('|p_{scat}|')
xlabel('ka')

break

figure(1)
plot(kr,engy_r)
xlabel('kr')
ylabel('ENERGY RATIO AT CUTOFF MODE n=ka')
axis([0 max(ka1) 0.5 1.1])
title(title_str)
grid

figure(2)
mesh(ka1,n,engy)
%axis([0 Nmax 0 Nmax])
ylabel('Order N')
xlabel('ka')
title(title_str)

break
figure(2)
colormap(jet)
imagesc(n,ka1,engy)
axis('square')
%axis([0 Nmax 0 Nmax])
ylabel('Order N')
xlabel('ka')
title(title_str)
colorbar


