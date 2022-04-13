% energy distribution

function   [engy, ratio]=energy(j,f,x)

f2=f.*conj(f);
indx=round(x)+1;
engy_tot=sum(f2);
n=length(f);
p(1)=f2(1);
for i=2:n
   p(i)=p(i-1)+f2(i);
end
engy=f2(:)/engy_tot;
if indx == 1
  ratio=p(1)/engy_tot;
else
 ratio=p(indx)/engy_tot;
end
