clc;        clear all;      close all;      format long

%% MoSi
x = [1.127461510506305  14.245125499784351   5.780781253337324...
        19.649609489208387   4.873159703434938  14.681656921944104   3.144758491240632   0.338559239864453   1.011683444359494];
wv1=linspace(0.2,25,2000);           omega=1.2398./wv1;
Drude=x(2)^2./(omega.^2 + 1i*omega*x(3));
Lorentz=x(4)^2./(x(5)^2 - omega.^2 - 1i*omega*x(6));
Lorentz1=x(7)^2./(x(8)^2 - omega.^2 - 1i*omega*x(9));
ncal2=sqrt(x(1)-Drude+Lorentz+Lorentz1);
figure(1);      subplot(1,3,[2 3]);
plot(wv1,real(ncal2),'Color',[0 0.4470 0.7410],'linewidth',1.3);     hold on;      
xlim([0.2 25]);     xticks([0.2 5:5:25]);  
ylim([1 22]);       yticks(1:3:22);
xlabel('Wavelength (μm)','FontSize',16);               
ylabel('Refractive index, n','FontSize',16);      
set(gca,'LineWidth',1.1,'fontsize',16);   

figure(2);      subplot(1,3,[2 3]);
plot(wv1,imag(ncal2),'Color',[0.8500, 0.3250, 0.0980],'linewidth',1.3);       hold on;           
xlim([0.2 25]);     xticks([0.2 5:5:25]); 
ylim([1 22]);       yticks(1:3:22);
xlabel('Wavelength (μm)','FontSize',16);               
ylabel('Extinction coefficient, k','FontSize',16);
set(gca,'LineWidth',1.1,'fontsize',16);   

% Python nk simulation @4K
load('MoSi_nk_4K.mat');
figure(1);      subplot(1,3,[2 3]);
plot(wv, real(nk_result),'--','Color',[0 0.4470 0.7410],'linewidth',1.3);     hold on;
legend({'FTIR model, T = 293 K','simulated, T = 4K'},'Location','best','FontSize',14);          
legend boxoff;      
figure(2);      subplot(1,3,[2 3]);
plot(wv, imag(nk_result),'--','Color',[0.8500, 0.3250, 0.0980],'linewidth',1.3);     hold on;
legend({'FTIR model, T = 293 K','simulated, T = 4K'},'Location','best','FontSize',14);          
legend boxoff;
