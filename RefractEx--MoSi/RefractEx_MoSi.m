clc;        clear;      close all;      base_dir = cd;      format long 

%% Load FTIR data
cd('FTIR_MoSi_Spectra/Transmission');      
[wvlen1,Spec1,~]=LoadSpectra('MoSi.spa',1);    % in percentage
cd(base_dir);   cd('FTIR_MoSi_Spectra/Veemax/R30');      
[wvlen2,Spec2,~]=LoadSpectra('15nm MoSi CaF2.spa',1);    %Spec2=Spec2*1e-2;
cd(base_dir);   cd('FTIR_MoSi_Spectra/Veemax/R45');      
[wvlen3,Spec3,~]=LoadSpectra('15nm MoSi CaF2.spa',1);    %Spec3=Spec3*1e-2;
cd(base_dir);   cd('FTIR_MoSi_Spectra/Veemax/R60');      
[wvlen4,Spec4,~]=LoadSpectra('15nm MoSi CaF2.spa',1);    %Spec4=Spec4*1e-2;

%% Measurement Parameters
cd(base_dir);                            wv1 = linspace(wvlen1(1), wvlen1(end), 600);    
omega = 1.2398 ./ wv1;          angles = [0 30 45 60];                   
polarization = 45;                     thickness_nom = [14.981e-3, 500];    % in microns
sigma_thickness = [0.65e-3, 50];     sigma_angle = 1;                   % ±1° perturbation
R_measured = [interp1(wvlen1, Spec1, wv1); interp1(wvlen2, Spec2, wv1);
              interp1(wvlen3, Spec3, wv1); interp1(wvlen4, Spec4, wv1)];
sigma_R = 2;        % in percentage

load('FTIR_MoSi_Spectra.mat')
nCaF2 = interp1(nk_CaF2(1,:), nk_CaF2(3,:), wv1, 'linear', 'extrap');
x = [2.1339    0.079068    0.024601   0.0048217];
Lorentz=x(2)^2./(x(3)^2 - omega.^2 - 1i*omega*x(4));
nsub = real(sqrt(x(1)+Lorentz));      % ksub = imag(sqrt(x(1)+Lorentz)); 
ksub = smoothdata(imag(nCaF2),'gaussian',100);    ncal = nsub + 1i * ksub;

%% Optimization Parameters
n_trials = 10;
x_results = zeros(n_trials, 9);
err_results = zeros(1,n_trials);
x0=[1.127461510506305  14.245125499784351   5.780781253337324...
       19.649609489208387   4.873159703434938  14.681656921944104   3.144758491240632   0.338559239864453   1.011683444359494];
options = optimset('Display','iter','TolX',1e-4,'TolFun', 1e-4,'MaxFunEvals',1e4,'MaxIter',300);
rng('default');

%% Monte Carlo Simulation
for trial = 1:n_trials
    % Perturb thickness and reflectance
    thickness = thickness_nom + sigma_thickness .* randn(1, 2);
    R_perturbed = R_measured + sigma_R .* randn(size(R_measured));
    angle_perturbed = angles + sigma_angle .* (2*rand(size(angles)) - 1);  
    % Define reduced chi-squared objective
    disp(['iteration ' num2str(trial)]);
    objfun = @(x) compute_chi2_fast(x, omega, wv1, R_perturbed, ncal, thickness, angle_perturbed, polarization, sigma_R);
    [x_fit, err_val] = fminsearch(objfun, x0, options);
    x_results(trial, :) = x_fit;
    err_results(trial) = err_val;
    x0 = x_fit;
end

%% Post-Processing
x_mean = mean(x_results, 1);
x_std = std(x_results, 0, 1);
disp('Fitted parameters (mean ± std):');
for i = 1:length(x_mean)
    fprintf('x(%d): %.5f ± %.5f\n', i, x_mean(i), x_std(i));
end

%% Calculate MoSi Refractive Indices
% N = length(omega);                      % number of spectral points
% n_trials = size(x_results, 1);          % number of Monte Carlo trials
% ncal_all = zeros(n_trials, N);          % store each computed ncal
% 
% for i = 1:n_trials
%     x = x_results(i, :);
%     Drude  = x(2)^2 ./ (omega.^2 + 1i * omega * x(3));
%     Lorentz1 = x(4)^2 ./ (x(5)^2 - omega.^2 - 1i * omega * x(6));
%     Lorentz2 = x(7)^2 ./ (x(8)^2 - omega.^2 - 1i * omega * x(9));
%     ncal_all(i, :) = sqrt(x(1) - Drude + Lorentz1 + Lorentz2);
% end
% 
% % Compute mean and std over the trials
% ncal_mean = mean(ncal_all, 1);
% ncal_std  = std(ncal_all, 0, 1);
% 
% figure(1); hold on;
% plot(wv1, real(ncal_mean), '-b', 'LineWidth', 1.2);
% plot(wv1, imag(ncal_mean), '-r', 'LineWidth', 1.2);
% fill([wv1 fliplr(wv1)], [real(ncal_mean + ncal_std) fliplr(real(ncal_mean - ncal_std))], ...
%      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill([wv1 fliplr(wv1)], [imag(ncal_mean + ncal_std) fliplr(imag(ncal_mean - ncal_std))], ...
%      'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% title('n and k with Uncertainty from x\_std'); xlabel('Wavelength (µm)'); ylabel('n, k');
% legend('n', 'k', 'n ± σ', 'k ± σ');

%% Use ncal_mean ± ncal_std to generate uncertainty in T_sim and R_sim
% mode_list = [1, 2, 2, 2];           % 1 = transmission, others = reflection
% figure(2); 
% 
% for j = 1:length(angles)
%     angle = angles(j);
%     mode = mode_list(j);
% 
%     % Mean simulation
%     nk_mean=[ones(1,length(wv1));ncal_mean;nCaF2;ones(1,length(wv1))];
%     [T_mean, R_mean, ~] = transfer_matrix(wv1, angle, polarization, thickness_nom, nk_mean);
%     nk_stdp=[ones(1,length(wv1));ncal_mean + ncal_std;nCaF2;ones(1,length(wv1))];
%     [T_upper, R_upper, ~] = transfer_matrix(wv1, angle, polarization, thickness_nom, nk_stdp);
%     nk_stdn=[ones(1,length(wv1));ncal_mean - ncal_std;nCaF2;ones(1,length(wv1))];
%     [T_lower, R_lower, ~] = transfer_matrix(wv1, angle, polarization, thickness_nom, nk_stdn);
% 
%     % Choose T or R
%     if mode == 1
%         sim = T_mean;
%         sim_upper = T_upper;
%         sim_lower = T_lower;
%         label = 'T (%)';
%     else
%         sim = R_mean;
%         sim_upper = R_upper;
%         sim_lower = R_lower;
%         label = 'R (%)';
%     end
% 
%     % Plot
%     subplot(length(angles), 1, j); hold on;
%     fill([wv1 fliplr(wv1)], [sim_upper fliplr(sim_lower)], ...
%         [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     plot(wv1, R_measured(j, :), 'k', 'LineWidth', 1.2);
%     plot(wv1, sim, 'r--', 'LineWidth', 1.3);
% 
%     xlim([min(wv1), max(wv1)]); grid on;
%     ylabel(label); title(sprintf('\\theta = %d^{\\circ}', angle));
%     legend('Model ± \sigma', 'Measured', 'Model', 'Location', 'best');
% end
% 
% xlabel('Wavelength (\mum)');
% sgtitle('Measured vs Simulated Spectra using n ± \sigma');

%% Chi-Squared Function
function chi2_val = compute_chi2_fast(x, omega, wv1, R_data, n_sub, thickness, angles, pol, sigma_R)
    N = length(wv1);
    n_angles = length(angles);
    chi2 = 0;
    n_k = sqrt(x(1) - (x(2)^2) ./ (omega.^2 + 1i * omega * x(3)) + ...
                     (x(4)^2) ./ (x(5)^2 - omega.^2 - 1i * omega * x(6)) + ...
                     (x(7)^2) ./ (x(8)^2 - omega.^2 - 1i * omega * x(9)));
    for j = 1:n_angles
        mode = (j ~= 1) + 1;                % 1 = transmission, 2 = reflection
        R_model = arrayfun(@(i) tmm_thinfilm(wv1(i), angles(j), pol, thickness, n_k(i), n_sub(i), mode), 1:N);
        chi2 = chi2 + sum(((R_data(j, :) - R_model) / sigma_R).^2);
    end
    chi2_val = chi2 / (N * n_angles-9);  % reduced chi-squared
end