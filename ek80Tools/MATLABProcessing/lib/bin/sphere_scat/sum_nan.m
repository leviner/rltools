%% function    psum=sum_nan(A,k)
%% summation which ignores all nan's 
%% if A is an 1D array, k is not necessary, if A is a matrix
%% k is optional. Without k or k = 1, A is summed over column,
%% and k = 2, summation is over rows

function    psum=sum_nan(A,k)

D=size(A);
if D(1) == 1 | D(2) == 1						% 1-D array
   [indx]=find(~isnan(A));
   psum=sum(A(indx));
else
   if nargin == 1
     k=1;				% default direction: sum over each colume
   end
   if k == 1
    for i=1:D(2)
       [indx]=find(~isnan(A(:,i)));
       if length(indx) > 0          
          psum(i)=sum(A(indx,i));
       else
          psum(i)=nan;
       end
    end
   else
    for i=1:D(1)
       [indx]=find(~isnan(A(i,:)));
       if length(indx) > 0
          psum(i)=sum(A(i,indx));
       else
          psum(i)=nan;
       end
   end
   psum=psum(:);
  end
end
