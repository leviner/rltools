%fishmodel.m

clear
global para
tmp=path;
i=findstr(tmp,'e:\revelle96\sildata\matfiles');
if isempty(i)
addpath(path,'e:\revelle96\sildata\matfiles')
end

tmp=path;
i=findstr(tmp,'e:\chu\zoomodel\wiebe');
if isempty(i)
   addpath(path,'e:\chu\zoomodel\wiebe')
   addpath(path,'e:\chu\zoomodel\wiebe\lib')
end

% Fluid-like Model (Fish using DWBA Bent Cylinder Model)
      %bd = 4;       % ratio of length to width of animal
      %d=0.25.*fish;           % width of fish body in size class assuming =.25*length
      para=set_para;				% obtain default parameters for Chu model computations
      load fishlen.dat
% Fluid-like Model (Fish using DWBA Bent Cylinder Model)
      %bd = 4;       % ratio of length to width of animal
      %d=0.25.*fish;           % width of fish body in size class assuming =.25*length
clear sigma_bs sigmabs
L=fishlen;  % (m) lengths (mm) of individual amphipods in 16 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            % lengths of different animal groups in mm
            L=reshape(L,length(L)*16,1); % Make data file into one row
            para.shape.L=max(L);

para.simu.model=1;
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=2;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.058;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.058;		% compressional sound speed contrast (outer layer)
para.shape.L_a =2*4.0; %2*ratio of length to width of animal 
para.shape.rho_L=3; 
para.simu.df=2;
nan_indx=find(isnan(L));
good_indx=find(~isnan(L));
L(nan_indx)=[];
para.simu.ka=para.simu.k.*L*1e-3./para.shape.L_a;

[sigma_bs_L, para]=zoo_bscat(para);
sigma_bs(1:length(para.simu.ka))=sigma_bs_L.*(L'*1e-3).^2;

sigmabs(nan_indx)=ones(size(nan_indx))*.0;
sigmabs(good_indx)=sigma_bs;
fishbs=sigmabs;

i=find(fishbs>0);
%h3=plot(para.simu.ka, 10*log10(fishbs(1,i)),'rs');
%title('Fish')

eisizebs_rev96salp
close all
figure
fishlen=reshape(fishlen,length(fishlen)*16,1);
fishlen=fishlen';
j=find(isnan(fishlen)==0);
h3=plot(fishlen(1,j), 10*log10(fishbs(1,j)),'rs');
hold on
fishbss=reshape(fishbss,length(fishbss)*16,1);
fishbss=fishbss';
k=find(fishbss>0);

plot(fishlen(1,k),10*log10(fishbss(1,k)),'bo');
title('Fish')
