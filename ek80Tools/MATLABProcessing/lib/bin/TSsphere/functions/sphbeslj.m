% spherical bessel function of the first kind
% function 	jn=sphbeslj(n,x)

function 	jn=sphbeslj(n,x)

indx=find(abs(x) < 1e-14);
x(indx)=1e-10*ones(size(indx));
x=x(:);
%Jm_1_2=bessel(n+1/2,x);
Jm_1_2=besselJ(n+1/2,x);

if (length(n) == 1)
  jn=sqrt(pi./(2*x)).*Jm_1_2;
  return
else
  for i=1:max(n+1)-min(n+1)+1
    jn(:,i)=sqrt(pi./(2*x)).*Jm_1_2(:,i);
  end
end
