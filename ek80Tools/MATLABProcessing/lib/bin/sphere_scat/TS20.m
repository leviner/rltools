clear 
close all

addpath c:\bin\matlab\toolboxes\seawater

ka_min=10;
ka_max=50;
R=20e-3./2;
proc_flag=1;%ka
scale_flag=1;%linear
out_flag=3;%scaled scat amplitude
cw=sw_svel(1,20,1)
g=14900./1000;
hR=6853/cw;
hT=4171/cw;
para_flag=[10000,ka_min,ka_max,g,hR,hT,180];

[outx,outy]=elastic_fs(proc_flag,scale_flag,out_flag,para_flag);
f=cw.*outx./2./pi./R;


figure(1)
plot(outx,10.*log10(outy./pi./2),'k')
set(gca,'linewidth',[2],'fontsize',[12])
xlabel('ka','fontsize',[12])
ylabel('REDUCED TARGET STRENGTH (dB)','fontsize',[12])
ylim([-40 0])


figure(2)
plot(f/1000,10.*log10(outy./2.*pi.*R.^2),'k','linewidth',2)
set(gca,'linewidth',[2],'fontsize',[12])
xlabel('FREQUENCY (kHz)','fontsize',[12])
ylabel('TARGET STRENGTH (dB)','fontsize',[12])
ylim([-70 -30])
grid on
axis([950 1050 -55 -35])