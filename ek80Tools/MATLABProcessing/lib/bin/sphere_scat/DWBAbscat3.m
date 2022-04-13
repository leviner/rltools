% scattering amplitude by a wealy scattering fluid umbrella-shaped
% jellyfish 
% fluid shell 3/10/01 notes
% study  ka (frequency) dependence using solution based on integral solution
% no restriction on incident angle: theta function  [ka0, ang, f]=DWBAscat2(fname,simu_para, shape_para,phy_para,prof_fname)
function	[ka0, ang, f]=umbrella_int(para)
% shape parameters
L_a=para.shape.L_a;
%a=para.shape.L/para.shape.L_a*1e-3;		% in mm

% simulation parameter
n=length(para.simu.ka);				% # of ka point
n_int=para.simu.ni;		% # of integration points

ang=para.simu.ang; 		% different incident angle
m=length(ang);
th=ang(:)*pi/180;

ka0=para.simu.ka;

g=para.phy.g1;
h=para.phy.hL;
Cb=(1-g.*h.*h)./(g.*h.*h)-(g-1)./g;
% construct other Matrices
Ka1=ka0/h;
Kb1=ka0*para.shape.b1_a1/h;
Th=th(:,ones(1,n_int));
Mu_z=2*cos(Th);   
Mu_r=2*sin(Th);   
term01=-i*0.25*Cb*Ka1(:,ones(1,m)).*(para.shape.L(:,ones(1,m))*1e-3 ...
        /(2*para.shape.a_a1))./(cos(th(:,ones(1,n))')+eps);
% variables related to the outer spheroid
r1=linspace(0,para.shape.a_a1,n_int);
dr1=r1(2)-r1(1);
R1=r1(ones(n,1),:);
Sr1=sqrt(1-R1.*R1);
% variables related to the inner spheroid
Ka2=ka0*para.shape.a2_a1/h;
Kb2=ka0*para.shape.b2_a2*para.shape.a2_a1/h;
r2=linspace(0,para.shape.a_a2,n_int);
dr2=r2(2)-r2(1);
R2=r2(ones(n,1),:);
Sr2=sqrt(1-R2.*R2);
term02=-i*0.25*Cb*Ka2(:,ones(1,m)).*(para.shape.a2_a1* ...
      para.shape.L(:,ones(1,m))*1e-3/(2*para.shape.a_a1))./ ...
      (cos(th(:,ones(1,n))')+eps);

Jarr=1:m;					  % angle index for checking purpose
for J=1:length(Jarr)        			  % angle loop
% disp(sprintf(' j = %g, stop =%g',J,status.stop))
    j=Jarr(J);
    mu_z=Mu_z(j,:);
    mu_r=Mu_r(j,:);
    Muz=mu_z(ones(n,1),:);
    Mur=mu_r(ones(n,1),:);
    Arg_r1=Mur.*Ka1(:,ones(1,n_int));
    Arg_r2=Mur.*Ka2(:,ones(1,n_int));
    Arg_b1=Muz.*Kb1(:,ones(1,n_int));
    Arg_b2=Muz.*Kb2(:,ones(1,n_int));
    J0r1=besselj(0,Arg_r1);
    J0r2=besselj(0,Arg_r2);
    EXP_b1=exp(i*Arg_b1.*Sr1);
    EXP_b2=exp(i*Arg_b2.*Sr2);
    term1=R1.*J0r1.*EXP_b1*dr1;
    term2=R2.*J0r2.*EXP_b2*dr2;
    f1(:,j)=sum(term1,2);
    f2(:,j)=sum(term2,2);
end
f=term01.*f1-term02.*f2;
