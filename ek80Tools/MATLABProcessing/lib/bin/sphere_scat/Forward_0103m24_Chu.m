%Eisizebs_antarct_DWBA.m
% mfile to compute estimated volume backscattering using the DWBA model for fluid like animals 
% and the 94 model for pteropods and the data from the MOCNESS net tows 
% taken on NBP0103 krill patch study  (silhouette data

clear all;close all

global para
tic
tmp=path;
i=findstr(tmp,'c:\projects\Forwardproblem');
if isempty(i)
   addpath(path,'c:\bin\object_scat')
   addpath(path,'c:\bin\matlab')
   addpath(path,'c:\projects\Forwardproblem\zoomodel\wiebe')
   addpath(path,'c:\projects\Forwardproblem\zoomodel\wiebe\lib')
end

loadfils_ant_0103moc24;  % load the length data files to compute expected sv

expect=[1e-12 1e-12;1e-4 1e-4];

    
para=set_para;				% obtain default parameters for Chu model computations

%-------------------------------------------------------------------------------
% Fluid-like Model (Small copepod using DWBA Bent Cylinder Model)
%copslen=[copslen(:,1:4),copslen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.
L=copslen;  % (m) lengths (mm) of individual copepods in 9 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            %  d=0.3922.*copepod;      % width of copepod body assuming =.3*length
            % lengths of different animal groups in mm
[r,c]=size(L);            
L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.model_index = 3;  % scat. model: 1-> BCM, 2->Medusae, 3->ellipsod
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.02;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.058;		% compressional sound speed contrast (outer layer)
%%% ellipsoidal shape 
para.shape.L_a =2*2.5497; %2*ratio of length to width of animal 
para.shape.L_b = para.shape.L_a; % ellipsoid => prolate spheroid if L_a = L_b
para.shape.phi=0;    % azimuthal angle in radians: 0 => ka, pi/2 => kb

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
copsbs=reshape(sigmabs,r,c);


% compute the small copepod contribution to volume backscattering for each sample

for i=1:8   % do samples 1 to 10
   nbp0103m24totcopsbs(i)=sum(copsbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(5,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered     
   
end
figure(1)

h3=subplot(3,3,1);
axis([1e-12 1e-5 1e-12 1e-5])
%set(sp1,'YTick',10^-10:10^-4.3:10^-4,'yticklabel',10^-10:10^-4.3:10^-4)
loglog(nbp0103m24totcopsbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Small Copepods')
hold on
% kiwi97_acc(2:10,3)=10*log10(2:10,3)-14.5;
% kiwi97_acc(2:10,3)=10^(kiwi97_acc(2:10,3)/10);
%xy=[totcopbs(1,:)',kiwi97_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
set(h3,'XTickLabel',[]);
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
toc

%--------------------------
%-------------------------------------------------------------------------------
% Fluid-like Model (large copepod using DWBA Bent Cylinder Model)
clear sigma_bs sigmabs
para=set_para;				% obtain default parameters for Chu model computations
%copllen=[copllen(:,1:4),copllen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

L=copllen;  % (m) lengths (mm) of individual copepods in 8 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            %  d=0.3922.*copepod;      % width of copepod body assuming =.3*length
            % lengths of different animal groups in mm
[r,c]=size(L);            

L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.model_index = 3; % scat. model: 1-> BCM, 2->Medusae, 3->ellipsod
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.02;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.058;		% compressional sound speed contrast (outer layer)
%%% ellipsoidal shape 
para.shape.L_a =2*2.5497; %2*ratio of length to width of animal 
para.shape.L_b = para.shape.L_a; % ellipsoid => prolate spheroid if L_a = L_b
para.shape.phi=0;    % azimuthal angle in radians: 0 => ka, pi/2 => kb

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
coplbs=reshape(sigmabs,r,c);



% compute the large copepod contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totcoplbs(i)=sum(coplbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(6,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered     
end
figure(1)
h3=subplot(3,3,2);loglog(nbp0103m24totcoplbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Large Copepods')
hold on
%xy=[totlgcopbs(1,:)',kiwi97_acc(:,3)]
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
set(h3,'XTickLabel',[]);
text(log10(nbp0103m24totcoplbs(1,:)'),log10(nbp0103m24_acc(:,3)),num2str(nbp0103m24_acc(:,4)));
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%------------------------------------------------------------------------------
toc

%-------------------------------------------------------------------------------
% Fluid-like Model (Large Euphausiid using DWBA Bent Cylinder Model)
clear sigma_bs sigmabs
para=set_para;				% obtain default parameters for Chu model computations
%eudelen=[eudelen(:,1:4),eudelen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

L=eudelen;  % (m) lengths (mm) of individual copepods in 8 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            %  d=0.18665.*euphausiid    % width of euphausiid body 
            % lengths of different animal groups in mm
[r,c]=size(L);            
L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.model_index = 1;
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=5.485e-4*nanmean(L)+1.002;		% density contrast: rho1/rhow (outer layer) Calculated from Chu's Antarctic g/h vs. TOTAL L regressions
para.phy.hL=5.942e-4*nanmean(L)+1.004;		% compressional sound speed contrast (outer layer)
para.shape.L_a =2*5.3576; %2*ratio of length to width of animal 
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
eudebs=reshape(sigmabs,r,c);



% compute the large euphausiid and decapod contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24toteudebs(i)=sum(eudebs(:,i))*(2^aliqnbp0103m24(3,i))*(aliqnbp0103m24(25,i)/aliqnbp0103m24(7,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered
end

h3=subplot(3,3,3); loglog(nbp0103m24toteudebs(1,:),nbp0103m24_acc(:,3),'ro'); 
title('Large euphausiids/decapods')
hold on
%xy=[toteuphbs(1,:)',kiwi97_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
title('Large euphausiids')
axis([1e-12 1e-5 1e-12 1e-5])
set(h3,'XTickLabel',[]);

%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-----------------------------------------------------
toc

% Fluid-like Model (Euphausiid using DWBA Bent Cylinder Model)
clear sigma_bs sigmabs
para=set_para;				% obtain default parameters for Chu model computations
%eudelen=[eudelen(:,1:4),eudelen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

L=euphslen;  % (m) lengths (mm) of individual copepods in 8 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            %  d=0.18665.*euphausiid    % width of euphausiid body 
            % lengths of different animal groups in mm
[r,c]=size(L);            
L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.model_index = 1;
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
para.phy.g1=1.01;		% density contrast: rho1/rhow (outer layer) This is Chu's best estimate for animals smaller than he sampled in the Antarctic
para.phy.hL=1.01;		% compressional sound speed contrast (outer layer)
para.shape.L_a =2*5.3576; %2*ratio of length to width of animal 
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
euphsbs=reshape(sigmabs,r,c);



% compute the small krill contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24toteuphsbs(i)=sum(euphsbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(8,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of crust larv in for each sample
   % splits           #squares/total    volume filtered     
end
figure(1)
h3=subplot(3,3,4);loglog(nbp0103m24toteuphsbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Small euphausiids')
hold on
%xy=[totclarbs(1,:)',kiwi97_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-14 1e-5 1e-14 1e-5])
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%--------------------------
toc
%-------------------------------------------------------------------------------




%-------------------------------------------------------------------------------
% Fluid-like Model (amphipod using DWBA Bent Cylinder Model)
       %amph=amphlen.*1e-3;            % (m) #s/length class of amphipod in 16 samples
       %d=0.33315.*amph;           % width of amphipod body in size class assuming =.3*length
clear sigma_bs sigmabs
para=set_para;				% obtain default parameters for Chu model computations
%amphlen=[amphlen(:,1:4),amphlen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

L=amphlen;  % (m) lengths (mm) of individual amphipods in 16 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            % lengths of different animal groups in mm
[r,c]=size(L);            
L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.model_index = 1;
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.058;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.058;		% compressional sound speed contrast (outer layer)
para.shape.L_a =2*3.0021; %2*ratio of length to width of animal 
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
amphbs=reshape(sigmabs,r,c);



% compute the amphipod contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totamphbs(i)=sum(amphbs(:,i))*(2^aliqnbp0103m24(3,i))*(aliqnbp0103m24(25,i)/aliqnbp0103m24(10,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered
end

h3=subplot(3,3,7);loglog(nbp0103m24totamphbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Amphipods')
hold on
%xy=[totamphbs(1,:)',kiwi97_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
set(h3,'XTickLabel',[]);
ylabel('Measured Acoustic volume backscattering (sv)')
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-----------------------------------------------------
toc
%-------------------------------------------------------------------------------
% Fluid-like Model (Fish using DWBA Bent Cylinder Model)
      %bd = 4;       % ratio of length to width of animal
      %d=0.25.*fish;           % width of fish body in size class assuming =.25*length
clear sigma_bs sigmabs
%fishlen=[fishlen(:,1:4),fishlen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

L=fishlen;  % (m) lengths (mm) of individual amphipods in 8 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            % lengths of different animal groups in mm
[r,c]=size(L);            
L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.03;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.03;		% compressional sound speed contrast (outer layer)
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
fishbs=reshape(sigmabs,r,c);


% compute the fish contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totfishbs(i)=sum(fishbs(:,i))*(2^aliqnbp0103m24(4,i))*(aliqnbp0103m24(25,i)/aliqnbp0103m24(16,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits          #squares/total    volume filtered
end

h3=subplot(3,3,5);loglog(nbp0103m24totfishbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Fish')
hold on
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
set(h3,'XTickLabel',[]);
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-------------------------------------------------------------------------------

toc
%-------------------------------------------------------------------------------
% Fluid-like Model (chaetognath using DWBA Bent Cylinder Model)
				%bd = 17.151;       % ratio of length to width of animal
clear sigma_bs sigmabs
%chaelen=[chaelen(:,1:4),chaelen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

L=chaelen;  % (m) lengths (mm) of individual chaetognaths in 8 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            % lengths of different animal groups in mm
[r,c]=size(L);            
L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.model_index = 1;
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.058;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.058;		% compressional sound speed contrast (outer layer)
para.shape.L_a =2*17.151; %2*ratio of length to width of animal 
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
chaebs=reshape(sigmabs,r,c);


% compute the chaetognath contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totchaebs(i)=sum(chaebs(:,i))*(2^aliqnbp0103m24(3,i))*(aliqnbp0103m24(25,i)/aliqnbp0103m24(11,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered
end

h3=subplot(3,3,6);loglog(nbp0103m24totchaebs(1,:),nbp0103m24_acc(:,3),'ro');
title('Chaetognaths')
hold on
%xy=[totchaebs(1,:)',nbp0103m24_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
set(h3,'XTickLabel',[]);
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-------------------------------------------------------------------------------

toc

%-------------------------------------------------------------------------------
% Fluid-like Model (ostracod using DWBA Bent Cylinder Model)
       %bd = 2.5497;       % ratio of length to width of animal
clear sigma_bs sigmabs
%ostrlen=[ostrlen(:,1:4),ostrlen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.
     
L=ostrlen;  % (m) lengths (mm) of individual larval crustacean in 8 samples (nans fill out short columns)
            % *1e-3 to convert lengths to meters not needed for this
            % lengths of different animal groups in mm
[r,c]=size(L);            
L=reshape(L,r*c,1); % Make data file into one column
para.shape.L=max(L);            
para.simu.model=1;
para.simu.model_index = 1;
para.simu.aveA_flag=1;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=30;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.001;    	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=3;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.03;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.03;		% compressional sound speed contrast (outer layer)
para.shape.L_a =2*2.5497; %2*ratio of length to width of animal 
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
ostrabs=reshape(sigmabs,r,c);



% compute the ostracod contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totostrbs(i)=sum(ostrabs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(9,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered     
   
end
figure(1)

sp1=subplot(3,3,8);
axis([1e-12 1e-5 1e-12 1e-5])
%set(sp1,'YTick',10^-10:10^-4.3:10^-4,'yticklabel',10^-10:10^-4.3:10^-4)
loglog(nbp0103m24totostrbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Ostracods')
xlabel('Estimated Acoustic volume backscattering (sv)')   
hold on
%xy=[totostrabs(1,:)',nbp0103m24_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-------------------------------------------------------------------------------
toc
%-----------------------------------------------------
% Spherical Elastic Shell Model (Pteropod) from Stanton et al 1994

f  = 120000;  % (Hz)
c  = 1500;    % (m/s)
r  = 0.5;     % reflection coef. from fit of curves to data
%r=(g*h-1)/(g*h+1)
%	d  = .0009;   % (m) equivalent spherical diameter of body
%limalen=[limalen(:,1:4),limalen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.
pter=limalen.*1e-3;            % (m) lengths of individual pteropods in 10 samples (nans fill out short columns)

% compute sigmabs for each individual pteropod using tim's Elastic Shell model

step0=pter.^4;
step1 = (ones(size(pter))+(25/10)*pi^4*f^4*c^(-4).*step0).^(-1);
step2 = ((25/144)*pi^4.*pter.^6*f^4*r^2*c^(-4)).*step1;
pterobs=step2;         % backscattering cross-section for each individual pteropod
[r,c]=size(pterobs);
pterobs=reshape(pterobs,r*c,1); % look for nans and replace with zeros
i=find(isnan(pterobs)); 
pterobs(i)=ones(size(i))*.0;
pterobs=reshape(pterobs,r,c);

% compute the pteropod contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totlimabs(i)=sum(pterobs(:,i))*(2^aliqnbp0103m24(3,i))*(aliqnbp0103m24(25,i)/aliqnbp0103m24(13,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered 
end
figure(1)
subplot(3,3,9);loglog(nbp0103m24totlimabs(1,:),nbp0103m24_acc(:,3),'ro');
title('Limacina Pteropods')
hold on
%xy=[totpterbs(1,:)',nbp0103m24_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-----------------------------------------------------
toc

%--------------------------
% Spherical elastic Shell Model (Eggs --change g and h to copepod values))

f  = 120000;  % (Hz)
c  = 1500;    % (m/s)
r  = 0.058;     % reflection coef. from fit of curves to data
%r=(g*h-1)/(g*h+1)
%	d  = .0009;   % (m) equivalent spherical diameter of body

egg=sphrlen.*1e-3;            % (m) lengths of individuals in 8 samples (nans fill out short columns)

% compute sigmabs for each individual blob using tim's Elastic Shell model  (probably not accurate)

step0=egg.^4;
step1 = (ones(size(egg))+(25/10)*pi^4*f^4*c^(-4).*step0).^(-1);
step2 = ((25/144)*pi^4.*egg.^6*f^4*r^2*c^(-4)).*step1;
eggbs=step2;         % backscattering cross-section for each individual bivalve
[r,c]=size(eggbs);
eggbs=reshape(eggbs,r*c,1); % look for nans and replace with zeros
i=find(isnan(eggbs)); 
eggbs(i)=ones(size(i))*.0;
eggbs=reshape(eggbs,r,c);

% compute the egg contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totsphrbs(i)=sum(eggbs(:,i))*(2^aliqnbp0103m24(3,i))*(aliqnbp0103m24(25,i)/aliqnbp0103m24(20,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered 
end
figure(2)
subplot(3,3,1);loglog(nbp0103m24totsphrbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Spheres')
hold on
%xy=[toteggbs(1,:)',nbp0103m24_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-------------------------------------------------------------------------------

toc
%----------------------------------------------------------------
% Gelatinous Model -- Siphonophore nectophore parts ( using Stanton pers. comm. model for jelly-like stuff)

f  = 120000;  % (Hz)
c  = 1500;    % (m/s)
%	l  = .00218;   % (m) ave length of body (main portion)
%	bd = 2.5497;       % ratio of length to width of animal
r  = 0.028;   %  0.028 is reflection coef. from Monger et al. in press for aequorea
%r=(g*h-1)/(g*h+1)
%	d  = .00028;   % (m) diamter of bell or length of bract or nectophore
%	s  = 0.0;    % relative standard deviation of length (no deviation for single individuals
%snectlen=[snectlen(:,1:4),snectlen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

necto=snectlen.*1e-3;  % (m) lengths of individuals in 8 samples (nans fill out short columns)
% *1e-3 to convert lengths to meters
d=necto;    % diameter of medusa bell
ra=d/2;	     % radius of medusa bell
% compute sigmabs for each individual gelatinous plankter using Gelatinous model
% step one--calculate radii
r1=d.*.25;
r2=d.*.75;	
%step two
necbs=.25*((r1.^2)+(r2.^2))*(r.^2);  % backscattering cross-section for each individual
%according to Stanton (personal contact--10/23/1997) for 
%animals in a narrow size bin at high 
%frequencies--120kHz is ok.
[r,c] = size(necbs);
necbs=reshape(necbs,r*c,1); % look for nans and replace with zeros
i=find(isnan(necbs));
necbs(i)=ones(size(i))*.0;
necbs=reshape(necbs,r,c);

% compute the nectphores contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totsiphnbs(i)=sum(necbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(18,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered      
end
figure(2)

sp1=subplot(3,3,2);
axis([1e-12 1e-5 1e-12 1e-5])
%set(sp1,'YTick',10^-10:10^-4.3:10^-4,'yticklabel',10^-10:10^-4.3:10^-4)
loglog(nbp0103m24totsiphnbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Nectphores')
hold on
%xy=[totnecbs(1,:)',nbp0103m24_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-------------------------------------------------------------------

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%%% Contribution of target strength of Pneumatophore at 120 kHz
% Start new way using lookup table
f  = 120000;  % (Hz)
c  = 1500;    % (m/s)
k=2*pi*f/c;
a=spneulen*0.5/1000;  %1 mm divided by 1000 to put into meters
[r,col]=size(a);            
a=reshape(a,r*col,1); % Make data file into one column

ka=k*a;
load C:\Projects\ForwardProblem\mfiles\siphon_tbl.dat;
RTS=interp1(siphon_tbl(:,1),siphon_tbl(:,3),ka);
TS=RTS + 10*log10(pi*a.^2); %a2 is the radius of bubble in meters
pneumbs=10.^(TS/10);

i=find(isnan(pneumbs)); % look for nans and replace with zeros
pneumbs(i)=ones(size(i))*.0;

pneumbs=reshape(pneumbs,r,col);

%%% compute the contribution of pneumatophores to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totpneumbs(i)=sum(pneumbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(15,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered      
end

figure(2)
sp1=subplot(3,3,3);

axis([1e-12 1e-5 1e-12 1e-5])
%set(sp1,'YTick',10^-10:10^-4.3:10^-4,'yticklabel',10^-10:10^-4.3:10^-4)
loglog(nbp0103m24totpneumbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Pneumatophores')
hold on
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
set(h3,'XTickLabel',[]);
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx%---------------



%---------------------------------------------------------------
% Gelatinous Model -- Medusae ( using Chu's jellyfish_DWBA Bent Cylinder model)
clear sigma_bs sigmabs

para=set_para;				% obtain default parameters for Chu model computations

L=medlen;  % lengths (mm) of individual medusae in 9 samples (nans fill out short columns)

[m n]=size(L);
L=reshape(L,m*n,1); % Make data file into one row
para.simu.model=1;
para.simu.model_index=2;% 1: elongated,  2:jellyfish
para.simu.aveA_flag=0;	% average over angle
para.shape.ang=0;			% mean incident angle in degree
para.shape.dang=20;		% standard deviation std(ang)
para.simu.nA=60;			% number of discrete angles to be averaged over
para.simu.aveA_PDF=1;	% PDF model: 1-uniform, 2-Gaussian, 
para.simu.aveL_flag=0;	% average over length
para.shape.Lstd=0.04;  % relative length standard deviation
para.simu.nL=60;			% number of discrete angles to be averaged over
para.simu.f0=120;       % start frequency ;
para.simu.fe=120;       % end frequency - in this case only one frequency used;
para.phy.g1=1.02;		% density contrast: rho1/rhow (outer layer)
para.phy.hL=1.02;		% compressional sound speed contrast (outer layer)
para.shape.L_a =2; 
para.shape.b1_a1=1.31;  % ratio of a1 to a2:  a1/a2, 
								% a1 is the outer radius and a2 is the inner radius                   
								% b1 is the semi-major(minor) axis of the outer
                        % spheroid, and b2 is the one of the inner spheroid.
								% The spheroid is revolved around b1 or b2 axis.
								%  
para.shape.b2_a2=0.9; % chu suggests this value. Monger came up with 0.268; % b2/a2
para.shape.a2_a1=1;		% a2/a1
para.shape.a_a1=1;		% a is the radius of the mean circular opening
para.shape.a_a2=1;
para.shape.rho_L=3; 
para.simu.df=2;
nan_indx=find(isnan(L));
good_indx=find(~isnan(L));
L(nan_indx)=[];
para.simu.k=2*pi*para.simu.fe*1e3/para.phy.cw;
para.simu.ka=para.simu.k.*L*1e-3./para.shape.L_a;
para.shape.L=L;						% in mm

[sigma_bs_L, para]=zoo_bscat(para);
sigma_bs(1:length(para.simu.ka))=sigma_bs_L;%./(pi*(L'*1e-3/2).^2); Commented this out because DWBAscat3 returns TS not RTS

sigmabs(nan_indx)=ones(size(nan_indx))*.0;
sigmabs(good_indx)=sigma_bs;
medbs=reshape(sigmabs,m,n);



%%% compute the medusae contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 9
nbp0103m24totmedbs(i)=sum(medbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(23,i))/aliqnbp0103m24(2,i);% splits           #squares/total    volume filtered    

end
figure(2)

sp1=subplot(3,3,4);
axis([1e-12 1e-5 1e-12 1e-5])
%%%set(sp1,'YTick',10^-10:10^-4.3:10^-4,'yticklabel',10^-10:10^-4.3:10^-4)
loglog(nbp0103m24totmedbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Medusae')
hold on
%%%xy=[totmedbs(1,4:8)',oc334_acc_200(:,3)];
%%%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
%%%stringr=['r=' num2str(r)];
%%%text(10^(-10),10^(-11),stringr,'fontsize',8);
%%%-------------------------------------------------------------------------------




%-------------------------------------------------------------------
% Gelatinous Model (Siphonophore bracts using Stanton pers. comm. model)

f  = 120000;  % (Hz)
c  = 1500;    % (m/s)
%	l  = .00218;   % (m) ave length of body (main portion)
%	bd = 2.5497;       % ratio of length to width of animal
r  = 0.028;   % 0.028 is reflection coef. from Monger et al. in press for Pleurobrachia
%r=(g*h-1)/(g*h+1)
%	d  = .00028;   % (m) diamter of bell or length of bract or nectophore
%	s  = 0.0;    % relative standard deviation of length (no deviation for single individuals
%sbractlen=[sbractlen(:,1:4),sbractlen(:,6:9)]; % get rid of moc2 net 8 and moc3 net 8 data - not used.

bract=sbractlen.*1e-3;  % (m) lengths of individual medusae in samples (nans fill out short columns)
% *1e-3 to convert lengths to meters
d=bract;    % diameter of medusa bell
%ra=d/2;	     % radius of medusa bell
% compute sigmabs for each individual gelatinous plankter using tim's model
% step one--calculate radii
r1=d.*.25;
r2=d.*.75;	
%step two
brcbs=.25*((r1.^2)+(r2.^2))*(r.^2);  % backscattering cross-section for each individual
%according to Stanton (personal contact--10/23/1997) for 
%animals in a narrow size bin at high 
%frequencies--120kHz is ok.
[r,c] = size(brcbs);
brcbs=reshape(brcbs,r*c,1); % look for nans and replace with zeros
i=find(isnan(brcbs));
brcbs(i)=ones(size(i))*.0;
brcbs=reshape(brcbs,r,c);

% compute the siphnophore bract contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totsiphbbs(i)=sum(brcbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(25,i)/aliqnbp0103m24(17,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered    
end
figure(2)

sp1=subplot(3,3,5);
axis([1e-12 1e-5 1e-12 1e-5])
%set(sp1,'YTick',10^-10:10^-4.3:10^-4,'yticklabel',10^-10:10^-4.3:10^-4)
loglog(nbp0103m24totsiphbbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Bracts')
xlabel('Estimated Acoustic volume backscattering (sv)') 
hold on
%xy=[totbrcbs(1,:)',nbp0103m24_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])
%stringr=['r=' num2str(r)];
%text(10^(-10),10^(-11),stringr,'fontsize',8);
%-------------------------------------------------------------------
toc

%-------------------------------------------------------------------------------
%%% Fluid-like Model (Polychaete using Bent Cylinder Model)

	f  = 120000;  % (Hz)
	c  = 1500;    % (m/s)
%%%	l  = .00218;   % (m) ave length of body (main portion)
	bd = 17.151;       % ratio of length to width of animal
	r  = 0.058;   % reflection coef. from fit of curves to data
%%%	d  = .00028;   % (m) width of body
	s  = 0.0;    % relative standard deviation of length

pol=polylen.*1e-3;            % (m) #s/length class of polychaete in 9 samples
d=0.0583.*pol;           % width of polychaete body in size class assuming =.1*length
%%% compute sigmabs for each polychaete size class using tim's Bent Cylinder Model

  step1 = cos(pi*f.*d.*c^(-1).*(4*(ones(size(pol)))-0.5*pi*(pi*f.*d.*c^(-1)+0.4*(ones(size(pol)))).^(-1)));
  step2 = (1-exp(-8*pi^2*f^2.*d.^2.*s^2*c^(-2)).*step1);
  polbs= (0.08*r^2.*pol.^2*bd^(-1)).*step2;  % backscattering cross-section for each size class

  [r,c] = size(polbs);
  polbs=reshape(polbs,r*c,1); % look for nans and replace with zeros
  i=find(isnan(polbs));
  polbs(i)=ones(size(i))*.0;
  polbs=reshape(polbs,r,c);

%%% compute the polychaete contribution to volume backscattering for each sample  
for i=1:8   % do samples 1 to 10
   nbp0103m24totpolbs(i)=sum(polbs(:,i))*2^aliqnbp0103m24(3,i)*(aliqnbp0103m24(end,i)/aliqnbp0103m24(12,i))/aliqnbp0103m24(2,i);        %*convertion to volume= total sv of copepods in for each sample
   % splits           #squares/total    volume filtered    
end
figure(2);
subplot(3,3,6);
axis([1e-12 1e-5 1e-12 1e-5])
%set(sp1,'YTick',10^-10:10^-4.3:10^-4,'yticklabel',10^-10:10^-4.3:10^-4)
loglog(nbp0103m24totpolbs(1,:),nbp0103m24_acc(:,3),'ro');
title('Polychaetes')
hold on
%xy=[totbrcbs(1,:)',nbp0103m24_acc(:,3)];
%fregres3
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-12 1e-5 1e-12 1e-5])


%-------------------------------------------------------------------------------
nbp0103m24totalbs = nbp0103m24totamphbs(1,:)+...
    nbp0103m24totchaebs(1,:)+...
    nbp0103m24totcoplbs(1,:)+...
    nbp0103m24totcopsbs(1,:)+...
    nbp0103m24toteuphsbs(1,:)+...
    nbp0103m24toteudebs(1,:)+...
    nbp0103m24totfishbs(1,:)+...
    nbp0103m24totlimabs(1,:)+...
    nbp0103m24totmedbs(1,:)+...
    nbp0103m24totostrbs(1,:)+...
    nbp0103m24totpolbs(1,:)+...
    nbp0103m24totsiphbbs(1,:)+...
    nbp0103m24totsiphnbs(1,:)+...
    nbp0103m24totsphrbs(1,:)+...    
    nbp0103m24totpneumbs(1,:);

figure(3)
loglog(nbp0103m24totalbs(1,:),nbp0103m24_acc(:,3),'ok');
%title('log total')
hold on
expect=[1e-11 1e-11;1e-4 1e-4];
loglog(expect(:,1),expect(:,2));  %expected relationship=[1e-12 1e-12;1e-5 1e-5]
axis([1e-11 1e-4 1e-11 1e-4])
xlabel('Model estimated sv')
ylabel('Acoustic volume backscattering (sv)')

%-------------------------------------------------------------------------------
figure(4)
expect1 =[-90 -90; -65 -65];
nbp0103m24xy=[10*log10(nbp0103m24totalbs(1,:))',10*log10(nbp0103m24_acc(:,3))];
z=nbp0103m24_acc(:,4);
[v,u,r,emat,fit]=fregres(nbp0103m24xy);

hold on
plot(expect1(:,1),expect1(:,2),'k');
axis([-90 -65 -90 -65])
%text(-82.5,-72.5,'r=0.55');
xlabel('Model estimated SV (in dB)')
ylabel('Acoustic volume backscattering (SV in dB)')
title('Comparison of Estimated SV and Measured SV')

% calculate percent contribution to backscattering cross-section
nbp0103m24pcentsv=[(nbp0103m24totcoplbs./nbp0103m24totalbs)*100;
         (nbp0103m24totcopsbs./nbp0103m24totalbs)*100;
         (nbp0103m24toteuphsbs./nbp0103m24totalbs)*100;
         (nbp0103m24toteudebs./nbp0103m24totalbs)*100;
         (nbp0103m24totfishbs./nbp0103m24totalbs)*100;
         (nbp0103m24totlimabs./nbp0103m24totalbs)*100;
         (nbp0103m24totsiphbbs./nbp0103m24totalbs)*100;
         (nbp0103m24totsiphnbs./nbp0103m24totalbs)*100;
         (nbp0103m24totmedbs./nbp0103m24totalbs)*100;
         (nbp0103m24totpolbs./nbp0103m24totalbs)*100;
         (nbp0103m24totchaebs./nbp0103m24totalbs)*100;
         (nbp0103m24totostrbs./nbp0103m24totalbs)*100;
         (nbp0103m24totpneumbs./nbp0103m24totalbs)*100;
         (nbp0103m24totamphbs./nbp0103m24totalbs)*100;
         (nbp0103m24totsphrbs./nbp0103m24totalbs)*100;];
 

 
nbp0103m24pcentsv_taxa =['copl';
    'cops';
    'eups';
    'eude';
    'fish';
    'lima';
    'sbra';
    'snec';
    'medu';
    'poly';
    'chae';
    'ostr'
    'pneu'
    'amph'
    'sphr']; %vector of taxon names in order.

figure(5)
pcolor(nbp0103m24pcentsv)
colorbar
xlabel('NET NUMBER')
ylabel('ANIMAL INDEX')
title('PERCENTAGE CONTRIBUTION')

nbp0103m24pcentsva=nbp0103m24pcentsv';
nbp0103m24pcentsv=[nbp0103m24_acc(:,4),10*log10(nbp0103m24_acc(:,3)),10*log10(nbp0103m24totalbs)',nbp0103m24pcentsv'];
figure(6)
plot(nbp0103m24pcentsv(:,2)-nbp0103m24pcentsv(:,3),'o-r');
hold on
plot([1 8],(mean(nbp0103m24pcentsv(:,2)-nbp0103m24pcentsv(:,3)))*[1 1],'--k');
title('Difference between Measured SV and Estimated SV')
xlabel('Sample Number')
ylabel('Measured SV - Estimated SV')
toc
%this figure plots up only those net estimates that are close to measured. 
net_indx=[1:7];
figure(7)
expect1 =[-90 -90; -65 -65];
msv=10*log10(nbp0103m24_acc(:,3));
estsv=nbp0103m24totalbs';
nbp0103m242xy_ts=[10*log10(estsv(net_indx,1)),  msv(net_indx,1)];

z=nbp0103m24_acc(:,4);
[v1,u1,r1,emat,fit]=fregres(nbp0103m242xy_ts);
hold on
plot(expect1(:,1),expect1(:,2),'k');
axis([-90 -65 -90 -65])
%text(-82.5,-72.5,'r=0.55');
xlabel('Model estimated SV (in dB)')
ylabel('Acoustic volume backscattering (SV in dB)')
title('Comparison of Estimated SV and Measured SV')

toc

%eisizebs_rev96dwba2 %do computation for Revelle96 data set 
%revkiwixy=[revxy;kiwi2xy]; %combine "good" data from both kiwi and revelle96 MOCNESS tows.
%[v1,u1,r1,emat,fit]=fregres(revkiwixy); %this calculates both sets of data that seem to have estimated values close to measured.

%mean(revkiwixy(:,2)-revkiwixy(:,1))
%measuredSV-estimated SV (for samples included in regression)
