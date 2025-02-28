clc;        clear;      close all;      
root = cd;             format long   

%% Load FTIR data
load('FTIR_MoSi_Spectra.mat');  
wv1=linspace(MoSi15_caf2_IR.wv(1),MoSi15_caf2_IR.wv(end),300); 
incidence_angle=[0 12 30 35 40 45 50 55 60];
polarization=45;
thickness=[14.981e-3 500];      % thickness in um
nCaF2=interp1(nk_CaF2(1,:),real(nk_CaF2(3,:)),wv1)+1i*interp1(nk_CaF2(1,:),imag(nk_CaF2(3,:)),wv1); 
T0=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.tran,wv1);          
R12=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref12,wv1); 
R30=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref30,wv1);     
R35=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref35,wv1);       
R40=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref40,wv1);     
R45=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref45,wv1);  
R50=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref50,wv1);     
R55=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref55,wv1);  
R60=interp1(MoSi15_caf2_IR.wv,MoSi15_caf2_IR.ref60,wv1);      

%% Perform Nonlinear Optimization with Termination Criterion
% E = hc/lambda, E (eV) = 1239.8 / lambda (nm)
err=@(x)0;       err1=@(x)0;         err2=@(x)0;         
err3=@(x)0;      err4=@(x)0;         err5=@(x)0;         
nk=cell(length(wv1),1);              N = length(wv1);

for k1=1:N
    omega = 1.2398/wv1(k1);    
    Drude = @(x) x(2)^2/(omega^2 + 1i*omega*x(3));
    Lorentz1 = @(x) x(4)^2/(x(5)^2 - omega^2 - 1i*omega*x(6));
    Lorentz2 = @(x) x(7)^2/(x(8)^2 - omega^2 - 1i*omega*x(9));
    eps = @(x) x(1)-Drude(x)+Lorentz1(x)+Lorentz2(x);
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
x0_initial=[1.127461510506305  14.245125499784351   5.780781253337324...
        19.649609489208387   4.873159703434938  14.681656921944104   3.144758491240632   0.338559239864453   1.011683444359494];
options=optimset('Display','iter','PlotFcns','optimplotfval','TolX',1e-10,'MaxFunEvals',9e7,'MaxIter',500); 
[x,err_val,exitflag,output]=fminsearch(RMSE,x0_initial,options); 
disp(['x: ', num2str(x)]);

% x0_initial=[1.1 14.2 5.78 19.65 4.87 14.6 3.14 0.34 1];
% x = [1.09928      15.0579      5.48423      20.2915      4.71809       15.005      3.16627     0.362083      0.95123];  ---> iterations 150, RMSE: 1.2898, R²: 0.99247
% x0 = [1.127461510506305  14.245125499784351   5.780781253337324...
%         19.649609489208387   4.873159703434938  14.681656921944104   3.144758491240632   0.338559239864453   1.011683444359494];
% x = [1.12134      15.1643      5.56605      23.4202      4.65586      13.1823      3.13649     0.370438     0.929163];  ---> iterations 500, RMSE: 1.2787, R²: 0.9926

%% Calculate MoSi Refractive Indices
omega=1.2398./wv1;  
Drude=x(2)^2./(omega.^2 + 1i*omega*x(3));
Lorentz1=x(4)^2./(x(5)^2 - omega.^2 - 1i*omega*x(6));
Lorentz2=x(7)^2./(x(8)^2 - omega.^2 - 1i*omega*x(9));
ncal=sqrt(x(1)-Drude+Lorentz1+Lorentz2);

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
