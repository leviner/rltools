% test cross section

clear

proc_flag=1;
scale=1;
out_flag=2;
para=[300 0.1 10.0 2.5 2.5 1.2 0];
%para=[300 0.1 10.0 2.5 2.5 0];
para=[300 0.1 5.0 0];

[outx, outy]=rgd_sft_fc(proc_flag,1,scale,out_flag,para);
%[outx, outy]=fluid_fc(proc_flag,scale,out_flag,para);
%[outx, outy]=elastic_fc(proc_flag,scale,out_flag,para);
plot(outx,0.5*sqrt(pi)*real(outy)./sqrt(outx));hold on

DEG2RAD=pi/180;
ns=50;
x0=0.1;
xe=5.0;
theta0=0;

   ka=linspace(x0,xe,ns);
   m=length(ka);
   Nmax=max(ka)+5;
   theta=theta0*DEG2RAD;
   n=0:Nmax;
   pn=cos(theta*n);
   em=2*[0.5 ones(1,Nmax)];
   dJn=cylbesldj(n,ka,0);
   dYn=cylbesldy(n,ka,0);
   for j=1:m
      gama=atan(dJn(j,:)./dYn(j,:));
      s=(em).*sin(gama).^2;
      sin_gama=dJn(j,:)./abs(dJn(j,:)+i*dYn(j,:));
      s1=-em.*sin(gama).*exp(-i*gama);
      s2=em.*sin_gama.^2;
      s3=em.*dJn(j,:)./(dJn(j,:)+i*dYn(j,:));
%disp(conj(s3'));pause
      I(j)=sum(s);
      I1(j)=sum(imag(s1));
      I2(j)=sum(s2);
      I3(j)=sum(real(s3));
   end
plot(ka,I./ka,ka,I1./ka,'+',ka,I2./ka,'o',ka,I3./ka,'r-')
hold off










break
proc_flag=1;
scale=1;
out_flag=2;
para=[300 0.1 10.0 2.5 2.5 1.2 0];
para=[300 0.1 10.0 0];

%[outx, outy]=elastic_fc(proc_flag,scale,out_flag,para);
%[outx, outy]=rgd_sft_fc(proc_flag,1,scale,out_flag,para);
[outx, outy]=rgd_sft_fc(proc_flag,1,scale,out_flag,para);
plot(outx,abs(4*pi*imag(outy./sqrt(outx))))
