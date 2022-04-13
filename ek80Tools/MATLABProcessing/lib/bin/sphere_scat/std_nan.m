%% function    pstd=std_nan(A,k)
%% compute the standard deviation of A which ignores all nan's 
%% if A is an 1D array, k is not necessary, if A is a matrix
%% k is optional. Without k or k = 1, A is averaged over column,
%% and k = 2, average is over rows

function    pstd=std_nan(A,k)

D=size(A);
if D(1) == 1 | D(2) == 1						% 1-D array
   [indx]=find(~isnan(A));
   if ~isempty(indx)
      pstd=std(A(indx));
   else
      pstd=nan;
   end
else
   if nargin == 1
     k=1;				% default direction: average over each colume
  end
  if k == 1
    for i=1:D(2)
       [indx]=find(~isnan(A(:,i)));
       if ~isempty(indx)
          pstd(i)=std(A(indx,i));
       else
          pstd(i)=nan;
       end
    end
  else
    for i=1:D(1)
       [indx]=find(~isnan(A(i,:)));
       if ~isempty(indx)
          pstd(i)=std(A(i,indx));
       else
          pstd(i)=nan;
       end
   end
   pstd=pstd(:);
  end
end
  
