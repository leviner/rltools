% scattering by an ellipsoidal weakly scatterer 
% f = complex scattering amplitude
% based on note 7/20/02 (Note Book: object scattering, p25)
% written by Dezhang Chu, Jan. 27, 2003

function   [ka0, ang, f]=DWBA_ellipsoid_fun(para)

% shape parameters
e_ca=para.shape.L_a/2;
e_cb=para.shape.L_b/2;
e_ba=e_cb/e_ca;

para.shape.L_a_ave=2*sqrt(e_ca*e_cb);

a=para.shape.L/para.shape.L_a*1e-3;
a_ave=para.shape.L/para.shape.L_a_ave*1e-3;	
para.simu.ka=para.simu.ka.*a_ave.*para.shape.L_a./(para.shape.L*1e-3);

% simulation parameter
n=length(para.simu.ka);				% # of ka point
ang=para.simu.ang; 		            % different incident angle
m=length(ang);
th=ang*pi/180;
phi=para.shape.phi;                 % azimuthal angle in radians

ka0=para.simu.ka;
coef=ka0.*ka0*e_ba*e_ca/para.shape.L_a;

g=para.phy.g1;
h=para.phy.hL;
Cb=(1-g.*h.*h)./(g.*h.*h)-(g-1)./g;


mux=2*sin(th)*cos(phi);
muy=2*sin(th)*sin(phi);
muz=2*cos(th);

%% prepare matrices to speed up computations
Ka=ka0(:,ones(1,m));
Coef=coef(:,ones(1,m));
Mux=mux(ones(n,1),:);
Muy=muy(ones(n,1),:);
Muz=muz(ones(n,1),:);

Mu=sqrt(Mux.*Mux+e_ba*e_ba*Muy.*Muy+e_ca*e_ca*Muz.*Muz);
Arg=(Ka/h).*Mu+eps;
for J=1:m
    Tmp(:,J)=sphbeslj(1,Arg(:,J))./Arg(:,J);
end
f=Cb*Coef.*Tmp;
ang=ang(:);

