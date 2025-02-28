clc;        clear all;      close all;      

%% NbN
x=[6.906980816237389   6.716385517669770   1.599998412858263...   
                    17.847903843066636   1.935406412562817  19.999981207736077]; 
% ---> RMSE: 0.29988, R²: 0.87879*
% x = [5.71166      6.76862      1.50626      16.6193      1.80117      19.3804];  ---> RMSE: 0.38637, R²: 0.79879
% x = [6.14949      6.69554      1.50756      16.6898      1.69575      19.8207];  ---> RMSE: 0.35762, R²: 0.82762
% x = [0.970595      6.75674        1.521      16.4046      1.72057       18.711];  ---> RMSE: 0.68622, R²: 0.3653
% x = [1.25145      6.74509      1.51497      16.2943      1.71901      19.4698];  ---> RMSE: 0.66828, R²: 0.39806
% x = [6.96956      6.83388      1.53481      18.1111      1.97749      20.2155];  ---> RMSE: 0.36931, R²: 0.81617
% x = [6.94437      6.88063      1.54415      17.9138      2.00503      20.1402];  ---> RMSE: 0.36839, R²: 0.81708
% x = [6.95739      6.88985      1.54701      18.0673      2.01179      20.6105];  ---> RMSE: 0.36818, R²: 0.81729

wv1=linspace(0.24,1.69,1000);           
omega=1.2398./wv1;
Drude=x(2)^2./(omega.^2 + 1i*omega*x(3));
Lorentz=x(4)^2./(x(5)^2 - omega.^2 - 1i*omega*x(6));
ncal1=sqrt(x(1)-Drude+Lorentz);

figure(1);      subplot(1,3,[2 3]);
yyaxis left;    plot(wv1,real(ncal1),'linewidth',1.3);     hold on;       
ylim([1 5]);        yticks(1:1:5);
xlim([0.2 1.7]);    xticks(0.2:0.3:1.7);                        
set(gca,'LineWidth',1.1,'fontsize',16);   
xlabel('Wavelength (μm)','FontSize',16);          
ylabel('Refractive index, n','FontSize',16);

yyaxis right;    plot(wv1,imag(ncal1),'linewidth',1.3);       hold on;           
ylim([1 6]);   yticks(1:1:6);
set(gca,'LineWidth',1.1,'fontsize',16); 
ylabel('Extinction coefficient, k','FontSize',16);          

%Semilab nearIR extracted data
wv31 = xlsread('NbN nk.csv',1,'A2:A1167')*1e-3;        
real_n31=xlsread('NbN nk.csv',1,'B2:B1167');     
imag_k31=xlsread('NbN nk.csv',1,'D2:D1167');   
real_n2 = interp1(wv31,real_n31,wv1);       
imag_k2 = interp1(wv31,imag_k31,wv1);

figure(1);      subplot(1,3,[2 3]);
yyaxis left;     plot(wv1,real_n2,'--','linewidth',1.3);     hold on;
yyaxis right;    plot(wv1,imag_k2,'--','linewidth',1.3);       hold on;    
legend({'FTIR model','Ellipsometry'},'Location','southeast','FontSize',16);          
legend boxoff;

% Add indicators for left y-axis
left_indicator_x = 0.5; 
left_indicator_y = 2.65;
rectangle('Position', [left_indicator_x-0.1, left_indicator_y-0.1, 0.1, 0.65], ...
                    'Curvature', [1, 1], 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1.3);
annotation('arrow', [0.49, 0.43], [0.45, 0.45], 'Color', [0 0.4470 0.7410], 'LineWidth', 1.3);

% Add indicators for right y-axis
right_indicator_x = 0.5; 
right_indicator_y = 1.5; 
rectangle('Position', [right_indicator_x-0.1, right_indicator_y-0.1, 0.1, 0.65], ...
                'Curvature', [1, 1], 'EdgeColor', [0.8500 0.3250 0.0980], 'LineWidth', 1.3);
annotation('arrow', [0.5, 0.56], [0.2, 0.2], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.3);
% ax = gca;   exportgraphics(ax,'Fig3a.pdf','ContentType','vector');


%% Calculate RMSE and R-squared (R²)
actual_data = [real_n2, imag_k2];                         
simulated_data = [real(ncal1), imag(ncal1)];  
RMSE = sqrt(mean((actual_data - simulated_data).^2));           
disp(['RMSE: ', num2str(RMSE)]);
SS_tot = sum((actual_data - mean(actual_data)).^2);          % Total sum of squares
SS_res = sum((actual_data - simulated_data).^2);                % Residual sum of squares
R2 = 1 - (SS_res / SS_tot);  disp(['R²: ', num2str(R2)]);          % R-squared formula