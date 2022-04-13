clear
%close all
clc

sys_indx=4;

[TS_tmp,f_tmp]=WC_TS(38.1,0);
[Psi,BP,freq]=Beam_Pattern_ftheta(sys_indx,0,0);

[min_ij,ij1]=min(abs(f_tmp-freq(1)));clear min_ij
[min_ij,ij2]=min(abs(f_tmp-freq(end)));clear min_ij

f=f_tmp(ij1:ij2);
TS=TS_tmp(ij1:ij2);
clear ij1 ij2 TS_tmp f_tmp

figure
set(gcf,'color','w','units','inches','position',[2, 2, 4, 4])
axes('units','inches','position',[0.5 0.4 3.2 3.5])

plot(f/1000,TS,'k','linewidth',2)
xlim([160 260])
set(gca,'fontsize',12,'linewidth',2)
%xlabel('Frequency (kHz)')
%ylabel('Target Strength (dB)')

ylabel('Target Strength (dB)','fontweight','bold')
xlabel('Frequelcy (kHz)','fontweight','bold')

axis([160 260 -65 -35])

[Psi,BP,freq]=Beam_Pattern_ftheta(sys_indx,1,0);
BP_interp=interp1(freq,BP,f);
hold on
%plot(f/1000,TS+BP_interp-nanmean(BP_interp(1:2)),'b')
plot(f/1000,TS+BP_interp,'g','linewidth',2)
hold off

[Psi,BP,freq]=Beam_Pattern_ftheta(sys_indx,2,0);
BP_interp=interp1(freq,BP,f);
hold on
%plot(f/1000,TS+BP_interp-nanmean(BP_interp(1:2)),'r')
plot(f/1000,TS+BP_interp,'r','linewidth',2)
hold off

[Psi,BP,freq]=Beam_Pattern_ftheta(sys_indx,3,0);
BP_interp=interp1(freq,BP,f);
hold on
%plot(f/1000,TS+BP_interp-nanmean(BP_interp(1:2)),'g')
plot(f/1000,TS+BP_interp,'m','linewidth',2)
hold off

[Psi,BP,freq]=Beam_Pattern_ftheta(sys_indx,4,0);
BP_interp=interp1(freq,BP,f);
hold on
%plot(f/1000,TS+BP_interp-nanmean(BP_interp(1:2)),'m')
plot(f/1000,TS+BP_interp,'y','linewidth',2)
hold off
l1=legend('0^o','1^o','2^o','3^o','4^o','location','southwest');
set(l1,'fontsize',10);

grid on

return
[Psi,BP,freq]=Beam_Pattern_ftheta(sys_indx,5,0);
BP_interp=interp1(freq,BP,f);
hold on
%plot(f/1000,TS+BP_interp-nanmean(BP_interp(1:2)),'c')
plot(f/1000,TS+BP_interp,'c')
hold off
