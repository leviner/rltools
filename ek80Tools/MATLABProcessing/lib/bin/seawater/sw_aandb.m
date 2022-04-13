function [a,b]=sw_aandb(S,T,P);

[m_S,n_S]=size(S);S=S(:);S=S';
if (m_S>1 & n_S>1)
    fprintf('No matrix S\n')
    return
end

[m_T,n_T]=size(T);T=T(:);T=T';
if (m_T>1 & n_T>1)
    fprintf('No matrix T\n')
    return
end

for ii=1:length(S)
    cw(ii,:) = sw_svel(S(ii)*ones(1,length(T)),T,P);
end


for ii=1:length(S)
    a(ii,:) = gradient(cw(ii,:),mean(diff(T)))./mean(mean(cw));
end
for ii=1:length(T)
    b(:,ii) = gradient(cw(:,ii),mean(diff(S)))./mean(mean(cw));
end
    

