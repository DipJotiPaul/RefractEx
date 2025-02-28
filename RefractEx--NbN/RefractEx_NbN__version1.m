clc;        clear;      close all;      
root = cd;             format long   

%% Load FTIR data
load('FTIR_NbN_Spectra.mat');  
wv1=linspace(wv(1),wv(end),300); 
incidence_angle=[0 12 30 35 40 45 50 55 60];
polarization=45;
thickness=[9.93e-3 500];           % thickness in um
nCaF2=interp1(nk_CaF2(1,:),real(nk_CaF2(3,:)),wv1)+1i*interp1(nk_CaF2(1,:),imag(nk_CaF2(3,:)),wv1); 
T0=interp1(wv,NbN_Tran,wv1);       R12=interp1(wv,NbN_R12,wv1); 
R30=interp1(wv,NbN_R30,wv1);       R35=interp1(wv,NbN_R35,wv1);       
R40=interp1(wv,NbN_R40,wv1);       R45=interp1(wv,NbN_R45,wv1);  
R50=interp1(wv,NbN_R50,wv1);       R55=interp1(wv,NbN_R55,wv1);  
R60=interp1(wv,NbN_R60,wv1);       

%% Perform Nonlinear Optimization with Termination Criterion
% E = hc/lambda, E (eV) = 1239.8 / lambda (nm)
err=@(x)0;       err1=@(x)0;         err2=@(x)0;         
err3=@(x)0;      err4=@(x)0;         err5=@(x)0;         
nk=cell(length(wv1),1);              N = length(wv1);

for k1=1:N
    omega = 1.2398/wv1(k1);    
    Drude = @(x) x(2)^2/(omega^2 + 1i*omega*x(3));
    Lorentz = @(x) x(4)^2/(x(5)^2 - omega^2 - 1i*omega*x(6));
    eps = @(x) x(1)-Drude(x)+Lorentz(x);
    nk{k1}=@(x) sqrt(eps(x));
    T0_fit=@(x) tmm_thinfilm(wv1(k1),incidence_angle(1),polarization,thickness,nk{k1}(x),nCaF2(k1),1);            
    err1=@(x) err1(x)+(T0(k1)-T0_fit(x))^2; 
    R12_fit=@(x) tmm_thinfilm(wv1(k1),incidence_angle(2),polarization,thickness,nk{k1}(x),nCaF2(k1),2);           
    err2=@(x) err2(x)+(R12(k1)-R12_fit(x))^2;
    R30_fit=@(x) tmm_thinfilm(wv1(k1),incidence_angle(3),polarization,thickness,nk{k1}(x),nCaF2(k1),2);           
    err3=@(x) err3(x)+(R30(k1)-R30_fit(x))^2;
    R45_fit=@(x) tmm_thinfilm(wv1(k1),incidence_angle(6),polarization,thickness,nk{k1}(x),nCaF2(k1),2);           
    err4=@(x) err4(x)+(R45(k1)-R45_fit(x))^2;
    R60_fit=@(x) tmm_thinfilm(wv1(k1),incidence_angle(9),polarization,thickness,nk{k1}(x),nCaF2(k1),2);           
    err5=@(x) err5(x)+(R60(k1)-R60_fit(x))^2;                         
end
MSE=@(x) (err1(x)+err2(x)+err3(x)+err4(x)+err5(x))/(5*N);
RMSE=@(x) sqrt(MSE(x)); 
x0_initial=[1.5   6.5   1.5  17.5   1.5  19.5];
options=optimset('Display','iter','PlotFcns','optimplotfval','TolX',1e-10,'MaxFunEvals',9e7,'MaxIter',500); 
[x,err_val,exitflag,output]=fminsearch(RMSE,x0_initial,options); 
disp(['x: ', num2str(x)]);

% x0_initial=[6.5   6.5   1.5  17.5   1.5  19.5];
% x = [5.71166      6.76862      1.50626      16.6193      1.80117      19.3804];  ---> iterations 150, RMSE: 1.7611, R²: 0.96115
% x = [6.14949      6.69554      1.50756      16.6898      1.69575      19.8207];  ---> iterations 500, RMSE: 1.7813, R²: 0.96027
% x0_initial=[1   6.5   1.5  17.5   1.5  19.5];
% x = [0.970595      6.75674        1.521      16.4046      1.72057       18.711];  ---> iterations 500, RMSE: 1.7676, R²: 0.96088
% x0_initial=[1.5   6.5   1.5  17.5   1.5  19.5];
% x = [1.25145      6.74509      1.51497      16.2943      1.71901      19.4698];  ---> iterations 500, RMSE: 1.7767, R²: 0.96047
% x0=[6.906980816237389   6.716385517669770   1.599998412858263...   
%                                   17.847903843066636   1.935406412562817  19.999981207736077];

%% Calculate NbN Refractive Indices
omega=1.2398./wv1;  
Drude=x(2)^2./(omega.^2 + 1i*omega*x(3));
Lorentz=x(4)^2./(x(5)^2 - omega.^2 - 1i*omega*x(6));
ncal=sqrt(x(1)-Drude+Lorentz);

%% Comparing FTIR data with the Simulated Data
nk_layer=[ones(1,length(wv1));ncal;nCaF2;ones(1,length(wv1))];
[T0_sim,R0_sim,~]=transfer_matrix(wv1,incidence_angle(1),polarization,thickness,nk_layer);
[T12_sim,R12_sim,~]=transfer_matrix(wv1,incidence_angle(2),polarization,thickness,nk_layer);
[T30_sim,R30_sim,~]=transfer_matrix(wv1,incidence_angle(3),polarization,thickness,nk_layer);
[T45_sim,R45_sim,~]=transfer_matrix(wv1,incidence_angle(6),polarization,thickness,nk_layer);
[T60_sim,R60_sim,~]=transfer_matrix(wv1,incidence_angle(9),polarization,thickness,nk_layer);
figure(1);        plot(wv1,real(ncal),'-b',wv1,imag(ncal),'-r','linewidth',1.1);      title('nk'); 
figure(2);        plot(wv1,T0,'-r',wv1,T0_sim,'--r','linewidth',1.1);                 title('Tran');     
figure(3);        plot(wv1,R12,'-r',wv1,R12_sim,'--r','linewidth',1.1);               title('R12');      
figure(4);        plot(wv1,R30,'-r',wv1,R30_sim,'--r','linewidth',1.1);               title('R30');           
figure(5);        plot(wv1,R45,'-r',wv1,R45_sim,'--r','linewidth',1.1);               title('R45');      
figure(6);        plot(wv1,R60,'-r',wv1,R60_sim,'--r','linewidth',1.1);               title('R60');      

%% Calculate RMSE and R-squared (R²)
actual_data = [T0, R12, R30, R45, R60];                         % Actual values for T0, R12, R30, R45, R60
simulated_data = [T0_sim, R12_sim, R30_sim, R45_sim, R60_sim];  % Simulated values for T0, R12, R30, R45, R60
RMSE = sqrt(mean((actual_data - simulated_data).^2));           
disp(['RMSE: ', num2str(RMSE)]);
SS_tot = sum((actual_data - mean(actual_data)).^2);             % Total sum of squares
SS_res = sum((actual_data - simulated_data).^2);                % Residual sum of squares
R2 = 1 - (SS_res / SS_tot);  disp(['R²: ', num2str(R2)]);       % R-squared formula
