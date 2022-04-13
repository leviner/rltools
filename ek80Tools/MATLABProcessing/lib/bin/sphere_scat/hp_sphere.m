function [sigmabs]=hp_sphere(a,ka,g,h,G,F);
%Stanton 1989

alpha_pis=(1-g.*h.^2)./3./g./h.^2  + (1-g)./(1+2.*g);
R=(g*h-1)./(g*h+1);

sigma_bs_num=a.^2.*(ka).^4.*alpha_pis.^2.*G;
sigma_bs_den=1+4.*(ka).^4.*alpha_pis.^2./R.^2./F;

sigmabs=sigma_bs_num./sigma_bs_den;

