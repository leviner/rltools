% derivative of the bessel function of the first kind
% function 	dJn=beslJd(n,x,flag,jn0)

%  flag = 0: compute jn0

%       = 1: do not compute jn0, with jn0 as input

%  djn(i)=i*jn(i)/x-jn(i+1)

function 	dJn=besselJd(n,x)

if nargin < 3

  flag=0;

end

indx=find(abs(x) < 1e-14);

x(indx)=1e-10*ones(size(indx));

x=x(:);

m=length(x);

Jn=besselj(min(n):max(n)+1,x);

if (length(n) == 1)
  dJn=n*Jn(:,1)./x-Jn(:,2);
  return
else
  for ii=1:max(n)-min(n)+1
    dJn(:,ii)=n(ii)*Jn(:,ii)./x-Jn(:,ii+1);
  end
end

