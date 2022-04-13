function		[f2,para]=zoo_bscat(para)
%% compute reduced backscattering differential cross section 
%% based on the parameters specified in structure "para"
%para.simu.model_index = 1;
para.simu.freq0=para.simu.f0:para.simu.df:para.simu.fe;
if para.simu.aveL_flag == 1
   ns=para.simu.min_ni;								% number of sample points per wave length
%   if length(para.shape.L) == 1
   para.simu.ka0=para.simu.ka;
   if length(para.simu.ka) == 1
   	para.simu.freq1=para.simu.f0*(1-3.1*para.shape.Lstd):para.simu.df ...
      	:para.simu.fe*(1+3.1*para.shape.Lstd)+para.simu.df;   
   	para.simu.k=2*pi*para.simu.freq1(:)*1e3/para.phy.cw;   
   	kLmax=max(para.simu.k)*para.shape.L*1e-3*1.3;		% 1+ 3 std(dL/L)
   	para.simu.ni=max(para.simu.min_ni,ceil(kLmax*ns/(2*pi)));		% integration points along z-axis
 %  para.simu.ka0=2*pi*para.simu.freq0*para.shape.L/para.shape.L_a/para.phy.cw;
   	para.simu.ka1=para.simu.freq1*para.shape.L/para.shape.L_a/para.phy.cw;                             % simulation ka
  elseif length(para.simu.freq0) == 1
    kLmax=max(para.simu.k*para.shape.L)*1e-3*1.3;		% 1+ 3 std(dL/L)
   	para.simu.ni=max(para.simu.min_ni,ceil(kLmax*ns/(2*pi)));		% integration points along z-axis
 %  para.simu.ka0=2*pi*para.simu.freq0*para.shape.L/para.shape.L_a/para.phy.cw;
 % 	para.simu.ka1=linspace(min(para.simu.ka0)*(1-3.1*para.shape.Lstd), ...                      
 %										max(para.simu.ka0)*(1+3.1*para.shape.Lstd),length(para.shape.L))';   % simulation ka
	para.simu.ka1=linspace(min(para.simu.ka0)*(1-3.1*para.shape.Lstd), ...                      
										max(para.simu.ka0)*(1+3.1*para.shape.Lstd),length(para.simu.ka))';   % simulation ka
  else
     disp('the program cannot handle both frequency and length arrays!')
     f2=-1;
     return
  end
  para.simu.ka=para.simu.ka1;  
else
   ns=para.simu.min_ni;								% number of sample points per wave length
   para.simu.k=2*pi*para.simu.freq0(:)*1e3/para.phy.cw; 
   if para.simu.model_index == 2
      kLmax=max(para.simu.k)*max(para.shape.L)*1e-3*1.3;		% 1+ 3 std(dL/L)
   else
      kLmax=max(para.simu.k)*para.shape.L*1e-3*1.3;		% 1+ 3 std(dL/L)
   end
   para.simu.ni=max(para.simu.min_ni,ceil(kLmax*ns/(2*pi)));		% integration points along z-axis
   para.simu.ka0=2*pi*para.simu.freq0*para.shape.L/para.shape.L_a/para.phy.cw;
   para.simu.ka1=para.simu.ka0;
end
switch  para.simu.model
  case 1						% weakly scattering fluid elongated objects
    if para.simu.aveA_flag == 1
	   para.simu.ang=linspace(para.shape.ang-3.1*para.shape.dang, ...
							  para.shape.ang+3.1*para.shape.dang,para.simu.nA);
    else
	   para.simu.ang=para.shape.ang;
    end
    switch para.simu.model_index
    case 1           % elongated bent euphausiids
	   [ka0, ang, f]=DWBAbscat1(para);
    case 2           % umbrella shaped jelly fish
	   [ka0, ang, f]=DWBAbscat3(para);
    case 3           % Ellipsodal shape copepods 
	   [ka0, ang, f]=DWBA_ellipsoid_fun(para);
    case 4           % Analytical approximation solution from echo average paper (Stanton et al, 1994)
	   [ka0, ang, f2]=DWBAbscat4(para);   
       return
    case 5           % Stanton's high freq. Kirchhoff approximation
	   [ka0, ang, f2]=DWBAbscat5(para);
       return
    end
    if para.simu.aveA_flag == 1
        orient_ave_para=[para.shape.ang para.shape.dang];
        f1=orient_ave(ang,f,para.simu.aveA_PDF,orient_ave_para);
    else
        f1=sqrt(f.*conj(f));
    end
  case 2						% scattering by bubbles
	 [ka0, f1]=Bubblescat1(para);    
  case 3						% scattering by shelled objects
    switch para.simu.model_index
    case 1
	   [ka0, f1]=Shellscat1(para);           
    case 2           % Stanton et al. 1994, analytical approximation
	   [ka0, ang, f2]=DWBAbscat6(para); 
       return
    end
  case 4						% weakly scattering fluid objects with spherical shell 
    if para.simu.aveA_flag == 1
	   para.simu.ang=linspace(para.shape.ang-3.1*para.shape.dang, ...
							  para.shape.ang+3.1*para.shape.dang,para.simu.nA);
    else
	   para.simu.ang=para.shape.ang;
    end
	 [ka0, ang, f]=DWBAbscat2(para);
    if para.simu.aveA_flag == 1
        orient_ave_para=[para.shape.ang para.shape.dang];
        f1=orient_ave(ang,f,para.simu.aveA_PDF,orient_ave_para);
    else
        f1=sqrt(f.*conj(f));
    end
 end
 if para.simu.aveL_flag == 1
  len_ave_para=[para.simu.nL para.shape.Lstd];
  f2=length_ave(para.simu.ka1,reshape(para.simu.ka0,1,length(para.simu.ka0)),f1,para.simu.aveL_PDF,len_ave_para);
  f2=f2.*f2;
 else
  f2=f1.*f1;
 end
f2=reshape(f2,1,length(f2));