function epsilon_m=neumann_N(m);

epsilon_m = ones(1,length(m)) ; 

if length(m)>1
    epsilon_m(2:length(m))=2 ;
end