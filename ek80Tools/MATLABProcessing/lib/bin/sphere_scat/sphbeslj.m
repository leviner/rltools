% spherical bessel function of the first kind
% function 	jn=sphbeslj(n,x)

function 	jn=sphbeslj(n,x)

indx=find(abs(x) < 1e-14);
x(indx)=1e-10*ones(size(indx));
x=x(:);
[n1,x1] = meshgrid(n,x);
size(n1);
size(x1);
Jm_1_2=besselj(n1+1/2,x1);
if (length(n) == 1)
  jn=sqrt(pi./(2*x1)).*Jm_1_2;
  return
else
  for ij=1:max(n+1)-min(n+1)+1
    jn(:,ij)=sqrt(pi./(2*x1(:,ij))).*Jm_1_2(:,ij);
  end
end
