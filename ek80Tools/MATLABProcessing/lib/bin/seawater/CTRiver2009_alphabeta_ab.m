clear
close all
clc

S=[2:1:30];T=0.1055*S+9.11;
P=5;%checked that over 1-20 range very little dependence
for ii=1:length(S)
    alph(ii,:) = sw_alpha(S(ii)*ones(1,length(T)),T,P);
    bet(ii,:)  = sw_beta(S(ii)*ones(1,length(T)),T,P);
    %cw(ii,:)   = sw_svel(S(ii)*ones(1,length(T)),T,P);
end
[a,b]=sw_aandb(S,T,P);
alph_m=mean(mean(alph))
bet_m=mean(mean(bet))
a_m=mean(mean(a))
b_m=mean(mean(b))

figure(1)
plot(S,alph./(alph(:,1)*ones(1,length(T))))
xlim([0 35])
%ylim([0 0.2])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('\alpha (T,S)','fontsize',16)

figure(1)
plot(S,alph.*1e3)
hold on
plot(15,alph_m.*1e3,'kd')
hold off
xlim([0 35])
ylim([0 0.2])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('\alpha (T,S)','fontsize',16)

figure(2)
plot(S,bet*1e3)
hold on
plot(15,bet_m.*1e3,'kd')
hold off
xlim([0 35])
ylim([0.75 0.8])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('\beta (T,S)','fontsize',16)

figure(3)
plot(S,a*1e3)
hold on
plot(15,a_m.*1e3,'kd')
hold off
xlim([0 35])
ylim([2.3 3.1])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('a(T,S)','fontsize',16)

figure(4)
plot(S,b*1e3)
hold on
plot(15,b_m.*1e3,'kd')
hold off
xlim([0 35])
ylim([0.8 0.87])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('b(T,S)','fontsize',16)

figure(5)
A2=(a+alph).^2;
plot(S,A2*1e6)
xlim([0 35])
ylim([4.5 8.5])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('A^2(T,S)','fontsize',16)

figure(6)
B2=(b+bet).^2;
plot(S,B2*1e6)
xlim([0 35])
ylim([2.44 2.62])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('B^2(T,S)','fontsize',16)

figure(7)
C=(A2*0.1055^2+B2+2*0.1055*sqrt(A2.*B2));
C_new=((a_m+alph_m).^2*0.1055^2+(b_m+bet_m).^2+2*0.1055*sqrt((a_m+alph_m).^2.*(b_m+bet_m).^2));
plot(S,C*1e6)
hold on
plot(15,C_new*1e6,'rd')
hold off
xlim([0 35])
%ylim([2.44 2.62])
set(gca,'linewidth',2,'fontsize',12)
xlabel('Salinity (psu)','fontsize',12)
ylabel('C','fontsize',16)

A2p=(3.19-0.13).^2.*1e-6;
B2p=(0.96+0.8).^2.*1e-6;
Cp=(A2p.*0.1055^2+2.*0.1055.*sqrt(A2p.*B2p)+B2p);%(5./6./(2).^(8/3))*

hold on
plot(33,Cp*1e6,'kd','markerfacecolor','k')
plot(15,3.49,'kd','linewidth',[1],'markerfacecolor','k')
hold off
