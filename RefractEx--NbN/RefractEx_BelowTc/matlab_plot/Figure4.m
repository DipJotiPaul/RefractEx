clc;        clear all;      close all;      format long

%% NbN
x=[6.906980816237389   6.716385517669770   1.599998412858263...   
            17.847903843066636   1.935406412562817  19.999981207736077];
wv1=linspace(0.2,25,2000);           omega=1.2398./wv1;
Drude=x(2)^2./(omega.^2 + 1i*omega*x(3));
Lorentz=x(4)^2./(x(5)^2 - omega.^2 - 1i*omega*x(6));
ncal1=sqrt(x(1)-Drude+Lorentz);
figure(1);      subplot(1,3,[2 3]);
plot(wv1,real(ncal1),'Color',[0 0.4470 0.7410],'linewidth',1.3);     hold on;      
xlim([0.2 25]);         xticks([0.2 5:5:25]); 
ylim([1 19]);           yticks(1:3:21);
xlabel('Wavelength (μm)','FontSize',16);               
ylabel('Refractive index, n','FontSize',16);
set(gca,'LineWidth',1.1,'fontsize',16);   

figure(2);      subplot(1,3,[2 3]);
plot(wv1,imag(ncal1),'Color',[0.8500, 0.3250, 0.0980],'linewidth',1.3);       hold on;           
xlim([0.2 25]);         xticks([0.2 5:5:25]);   
ylim([1 19]);           yticks(1:3:21);
xlabel('Wavelength (μm)','FontSize',16);               
ylabel('Extinction coefficient, k','FontSize',16);     
set(gca,'LineWidth',1.1,'fontsize',16);   

% Python nk simulation @4K
load('NbN_nk_4K.mat');
figure(1);      subplot(1,3,[2 3]);
plot(wv, real(nk_result),'--','Color',[0 0.4470 0.7410],'linewidth',1.3);     hold on;
legend({'FTIR model, T = 293 K','simulated, T = 4K'},'Location','best','FontSize',14);          
legend boxoff;  
figure(2);      subplot(1,3,[2 3]);
plot(wv, imag(nk_result),'--','Color',[0.8500, 0.3250, 0.0980],'linewidth',1.3);     hold on;
legend({'FTIR model, T = 293 K','simulated, T = 4K'},'Location','best','FontSize',14);          
legend boxoff;
