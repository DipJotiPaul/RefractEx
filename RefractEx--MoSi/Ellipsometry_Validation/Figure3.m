clc;        clear all;      close all;      

%% MoSi
x = [1.127461510506305  14.245125499784351   5.780781253337324...
        19.649609489208387   4.873159703434938  14.681656921944104   3.144758491240632   0.338559239864453   1.011683444359494]; 
% ---> RMSE: 0.15836, R²: 0.97118*
% x = [1.09928      15.0579      5.48423      20.2915      4.71809       15.005      3.16627     0.362083      0.95123]; ---> RMSE: 0.4461, R²: 0.77128
% x = [1.12134      15.1643      5.56605      23.4202      4.65586      13.1823      3.13649     0.370438     0.929163];---> RMSE: 0.59775, R²: 0.58934
% x = [4.9382      7.3189      1.3339      3.6639     0.34858     0.95085];---> RMSE: 1.6491, R²: -2.1257
% x = [1.1119      7.3053       1.331      3.6454     0.34447     0.96744];---> RMSE: 1.7287, R²: -2.4348

wv1=linspace(0.2,1.7,1000);           
omega=1.2398./wv1;
Drude=x(2)^2./(omega.^2 + 1i*omega*x(3));
Lorentz1=x(4)^2./(x(5)^2 - omega.^2 - 1i*omega*x(6));
Lorentz2=x(7)^2./(x(8)^2 - omega.^2 - 1i*omega*x(9));
ncal2=sqrt(x(1)-Drude+Lorentz1+Lorentz2);

figure(1);      subplot(1,3,[2 3]);
yyaxis left;    plot(wv1,real(ncal2),'linewidth',1.3);     hold on;       
xlim([0.2 1.7]);  xticks(0.2:0.3:1.7);  
ylim([1 6]);      yticks(1:1:6);
xlabel('Wavelength (μm)','FontSize',16);          
ylabel('Refractive index, n','FontSize',16);
set(gca,'LineWidth',1.1,'fontsize',16);   
yyaxis right;    plot(wv1,imag(ncal2),'linewidth',1.3);       hold on;           
ylim([1 6]);   yticks(1:1:6);
ylabel('Extinction coefficient, k','FontSize',16);          
set(gca,'LineWidth',1.1,'fontsize',16);       

%Semilab nearIR extracted data
omega2 = 1.2398./wv1;
Eps_inf = 1.08758;      E_p = 19.88504;       E_gamma = 9.68516;        
f = 8.7975;     E0 = 3.15189;   gamma = 4.8756;     
Drude2 = E_p^2./(omega2.^2+1i*omega2.*E_gamma);         
Lorentz2 = f*E0^2./(E0^2-omega2.^2-1i*omega2.*gamma); 
nk2=sqrt(Eps_inf - Drude2 + Lorentz2);          
real_n2 = real(nk2);            
imag_k2 = imag(nk2);
figure(1);      subplot(1,3,[2 3]);
yyaxis left;    plot(wv1,real_n2,'--','linewidth',1.3);     hold on;
yyaxis right;    plot(wv1,imag_k2,'--','linewidth',1.3);    hold on;    
legend({'FTIR model','Ellipsometry'},'Location','southeast','FontSize',16);          
legend boxoff;

% Add indicators for left y-axis
left_indicator_x = 1.1; 
left_indicator_y = 4.65;
rectangle('Position', [left_indicator_x-0.1, left_indicator_y-0.1, 0.1, 0.65], ...
                    'Curvature', [1, 1], 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1.3);
annotation('arrow', [0.69, 0.63], [0.77, 0.77], 'Color', [0 0.4470 0.7410], 'LineWidth', 1.3);

% Add indicators for right y-axis
right_indicator_x = 1.2; 
right_indicator_y = 4.1; 
rectangle('Position', [right_indicator_x-0.1, right_indicator_y-0.1, 0.1, 0.65], ...
                'Curvature', [1, 1], 'EdgeColor', [0.8500 0.3250 0.0980], 'LineWidth', 1.3);
annotation('arrow', [0.73, 0.79], [0.62, 0.62], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.3);

%% Calculate RMSE and R-squared (R²)
actual_data = [real_n2, imag_k2];                         
simulated_data = [real(ncal2), imag(ncal2)];  
RMSE = sqrt(mean((actual_data - simulated_data).^2));           
disp(['RMSE: ', num2str(RMSE)]);
SS_tot = sum((actual_data - mean(actual_data)).^2);          % Total sum of squares
SS_res = sum((actual_data - simulated_data).^2);                % Residual sum of squares
R2 = 1 - (SS_res / SS_tot);  disp(['R²: ', num2str(R2)]);          % R-squared formula