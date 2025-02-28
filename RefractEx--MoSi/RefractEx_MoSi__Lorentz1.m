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
    eps = @(x) x(1)-Drude(x)+Lorentz1(x);
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
x0_initial=[1.1 6.93 1.38 3.64 0.34 0.98];
options=optimset('Display','iter','PlotFcns','optimplotfval','TolX',1e-10,'MaxFunEvals',9e7,'MaxIter',500); 
[x,err_val,exitflag,output]=fminsearch(RMSE,x0_initial,options); 
disp(['x: ', num2str(x)]);

% x0_initial=[4.84 6.93 1.38 3.64 0.34 0.98];
% x = [4.9382      7.3189      1.3339      3.6639     0.34858     0.95085];  ---> iterations 150, RMSE: 1.2989, R²: 0.99236
% x0_initial=[1.1 6.93 1.38 3.64 0.34 0.98];
% x = [1.1119      7.3053       1.331      3.6454     0.34447     0.96744];  ---> iterations 500, RMSE: 1.2987, R²: 0.99236

%% Calculate MoSi Refractive Indices
omega=1.2398./wv1;  
Drude=x(2)^2./(omega.^2 + 1i*omega*x(3));
Lorentz1=x(4)^2./(x(5)^2 - omega.^2 - 1i*omega*x(6));
ncal=sqrt(x(1)-Drude+Lorentz1);

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
