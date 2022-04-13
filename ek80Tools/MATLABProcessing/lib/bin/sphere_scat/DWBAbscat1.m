% reduced scattering amplitude by a wealy scattering uniformly bent finite 
% inhomogenous (segmented) g and h along the body
% fluid cylinder  4/27/92 notes
% study  ka (frequency) dependence using solution based on integral solution
% no restriction on incident angle: theta function  [ka0, ang, f]=DWBAscat2(fname,simu_para, shape_para,phy_para,prof_fname)
function  [ka0, ang, f]=DWBAscat1(para)

% shape parameters
L_a=para.shape.L_a;
a=para.shape.L/para.shape.L_a*1e-3;		% in mm

% simulation parameter
n=length(para.simu.ka);				% # of ka point
n_int=para.simu.ni;		% # of integration points

ang=para.simu.ang; 		% different incident angle
m=length(ang);
th=ang(:)*pi/180;

ka0=para.simu.ka;

if para.phy.gFlag == 1   % regression
    g=para.phy.g_L;
    h=para.phy.h_L;
else
    g=para.phy.g1;
    h=para.phy.hL;
end
Cb=(1-g.*h.*h)./(g.*h.*h)-(g-1)./g;

% construct postion vectors
[r_pos,th_tilt,dr,gamma_t,taper,xp,zp]=buildpos1(para);

% construct other Matrices
X1=ka0*taper;
Tmp=h*ones(n,para.simu.ni);
X2=X1./Tmp;
Dtheta=th_tilt(:,ones(1,m))'- th(:,ones(1,n_int));
Cos_dtheta=abs(cos(Dtheta));   % choose different local coordinates to
                               % avoid negtive argument of Bessel function
Gamma_t=gamma_t(:,ones(1,m))';
Theta=th(:,ones(1,n_int));
Dgamma=Gamma_t-Theta;
Cos_dgamma=cos(Dgamma);
term0=L_a*ka0*(r_pos./h)';
term1=h.*h.*Cb.*dr/4;
Jarr=1:m;					  % angle index for checking purpose
for J=1:length(Jarr)        			  % angle loop
% disp(sprintf(' j = %g, stop =%g',J,status.stop))
    j=Jarr(J);
    cos_dtheta=Cos_dtheta(j,:);
    Cos_th=cos_dtheta(ones(n,1),:);
    cos_dgamma=Cos_dgamma(j,:);
    Cos_gamma=cos_dgamma(ones(1,n),:);
    Arg=2*X2.*Cos_th+eps;
    J1x=besselj(1,Arg)./Arg;
    EXP=exp(i*term0.*Cos_gamma);
    term2=(X2.*X2).*J1x.*EXP;
    f(:,j)=term2*term1+eps;
end

