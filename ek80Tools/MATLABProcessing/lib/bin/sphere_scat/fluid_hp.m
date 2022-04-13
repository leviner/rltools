function [sigma_bs]=fluid_hp(ka,g,h,G,F)
%Calculate sigma_bs/a.^2 for Stanton's 1989 high pass fluid sphere model

if isempty(G)
    G=1;
end
if isempty(F)
    F=1;
end

alpha_pis=(1-g*h.^2)./(3*g*h.^2) + (1-g)./(1+2*g);
alpha_pis2=alpha_pis.^2;
R=(g.*h-1)./(g.*h+1);

sigma_num=(ka).^4.*alpha_pis2;
sigma_den=1+4*(ka).^4.*alpha_pis2./R.^2;

sigma_bs=sigma_num./sigma_den;