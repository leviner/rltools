% build vectors related to postoins of the animal
% body axis is alone z-axis
% construct shape coordinates
% 
%function    [r_pos,th_tilt,dr,gamma_t,taper,x,z]=buildpos(n_int,prof_no, order,rho_L,fname,rough_axis)
function    [r_pos,th_tilt,dr,gamma_t,taper,x,z]=buildpos1(para)

n_int=para.simu.ni;
order=para.shape.order;
rho_L=para.shape.rho_L;
L_a=para.shape.L_a;

% uniformly bent cylinder and regularly tapered
  gamma=0.5/rho_L;
  ratio=2*rho_L;
  z=linspace(-1,1,n_int);
  taper=sqrt(1-z(:).^order);
% normalized by rho -  radius of curvature
  z=sin(gamma)*z(:);
  x=(1-sqrt(1-z.*z));
% normalized by L/2
  x=ratio*x;
  z=ratio*z;

x=x(:);z=z(:);taper=reshape(taper,1,n_int);
th_tilt=zeros(n_int,1);
r_pos=sqrt(x.*x+z.*z);
gamma_t=atan2(z,x);
dx=diff(x)+eps;
dz=diff(z);
alpha_t=[atan(dz./dx); atan(dz(n_int-1)/dx(n_int-1))];
indx1=find(alpha_t < 0);
if length(indx1) > 0
  th_tilt(indx1)=alpha_t(indx1)+pi/2;
end
indx2=find(alpha_t >= 0);
if length(indx2) > 0
  th_tilt(indx2)=alpha_t(indx2)-pi/2;
end
dr1=sqrt(dx.*dx+dz.*dz);
dr=[dr1(1); dr1];

disp_prof=0;
if disp_prof == 1
  itis=666;
  figure(2)
  z1=z+taper(:)/L_a;
  z2=z-taper(:)/L_a;
  z10=z+taper0(:)/L_a;
  z20=z-taper0(:)/L_a;
  plot(x,z,x,z1,'r',x,z2,'r',x,z10,'g',x,z20,'g');
  axis('equal');grid;
  disp('press any key to cont. ..');pause
 % plot(z,taper0,z,taper,'r');
  close(2)
end
%taper