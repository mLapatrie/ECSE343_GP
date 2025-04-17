clear; clc; close all;
%% Variable definition

% Nominal Component Values
R1_nom = 30e3; % R1 nominal value in ohms
R2_nom = 18e3; % R2 nominal value in ohms
C1_nom = 0.01e-6; % C1 nominal value in farads
C2_nom = 0.0047e-6; % C2 nominal value in farads

% Tolerance
tol = 0.05;

% Bounds for component based on tolerance
R1_lower = R1_nom * (1 - tol);   R1_upper = R1_nom * (1 + tol);
R2_lower = R2_nom * (1 - tol);   R2_upper = R2_nom * (1 + tol);
C1_lower = C1_nom * (1 - tol);   C1_upper = C1_nom * (1 + tol);
C2_lower = C2_nom * (1 - tol);   C2_upper = C2_nom * (1 + tol);

% Simulation parameters
N = 10000; % Number of Monte Carlo samples
f = 1000; % Frequency of the input signal in Hz
Vin_amp = 1; % Input amplitude in Volts

% Array for storing |Vout|
Vout_samples = zeros(N,1);

%% Section 2: Monte Carlo Simulation Loop
for i = 1:N
    % Randomly sample components from uniform distribution
    R1_sample = R1_lower + (R1_upper - R1_lower) * rand();
    R2_sample = R2_lower + (R2_upper - R2_lower) * rand();
    C1_sample = C1_lower + (C1_upper - C1_lower) * rand();
    C2_sample = C2_lower + (C2_upper - C2_lower) * rand();
    
    % Evaluate output using Sallen-Key simulation function
    Vout = simulate_Sallenkey(R1_sample, R2_sample, C1_sample, C2_sample, Vin_amp, f);
    
    % Store the magnitude of Vout
    Vout_samples(i) = abs(Vout);
end

%% Section 3: Statistical analysis

% Basic information 
mean_Vout = mean(Vout_samples);
std_Vout  = std(Vout_samples);
variance_Vout = var(Vout_samples);
median_Vout = median(Vout_samples);
min_Vout = min(Vout_samples);
max_Vout = max(Vout_samples);

% Percentiles (5th, 25th, 50th, 75th, 95th)
percentiles = prctile(Vout_samples, [5 25 50 75 95]);

% Skewness
skewness_Vout = skewness(Vout_samples);

% Kurtosis
kurtosis_Vout = kurtosis(Vout_samples);

% 95% Confidence Interval for the Mean
alpha = 0.05;
z = norminv(1 - alpha/2);  % Two-tailed (1.96 for 95% CI)
CI_lower = mean_Vout - z * std_Vout/sqrt(N);
CI_upper = mean_Vout + z * std_Vout/sqrt(N);

% Errors
StandardError = std_Vout / sqrt(N);
RelativeErrorPercent = (StandardError / mean_Vout) * 100;

%% Display the computed statistics
fprintf('\nMonte Carlo Simulation Results (N = %d samples):\n', N);
fprintf('---------------------------------------------\n');
fprintf('Mean of |Vout|                 : %.4f V\n', mean_Vout);
fprintf('Standard Deviation of |Vout|    : %.4f V\n', std_Vout);
fprintf('Variance of |Vout|              : %.4f V^2\n', variance_Vout);
fprintf('Median of |Vout|                : %.4f V\n', median_Vout);
fprintf('Minimum of |Vout|              : %.4f V\n', min_Vout);
fprintf('Maximum of |Vout|              : %.4f V\n', max_Vout);
fprintf('5th Percentile                 : %.4f V\n', percentiles(1));
fprintf('25th Percentile                : %.4f V\n', percentiles(2));
fprintf('50th Percentile (Median)       : %.4f V\n', percentiles(3));
fprintf('75th Percentile                : %.4f V\n', percentiles(4));
fprintf('95th Percentile                : %.4f V\n', percentiles(5));
fprintf('Skewness                       : %.4f\n', skewness_Vout);
fprintf('Kurtosis                       : %.4f\n', kurtosis_Vout);
fprintf('95%% Confidence Interval for Mean: [%.4f, %.4f] V\n', CI_lower, CI_upper);
fprintf('Standard Error of the Mean     : %.4f V\n', StandardError);
fprintf('Relative Error (%%)             : %.2f%%\n', RelativeErrorPercent);

%% Section 4: Plots

% Estimated PDF
nbins = 50;
[counts, binCenters] = hist(Vout_samples, nbins);
pdf_Vout = counts / trapz(binCenters, counts);

figure;
bar(binCenters, pdf_Vout, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
xlabel('|Vout| (V)');
ylabel('Probability Density');
title('Estimated PDF of |Vout| via Monte Carlo Simulation');
grid on;

% Empirical Cumulative Distribution Function (ECDF)
figure;
[f_ecdf, x_ecdf] = ecdf(Vout_samples);
stairs(x_ecdf, f_ecdf, 'LineWidth', 2);
xlabel('|Vout| (V)');
ylabel('Cumulative Probability');
title('Empirical Cumulative Distribution Function of |Vout|');
grid on;

% Boxplot for data spread
figure;
boxplot(Vout_samples, 'Notch','on', 'Labels', {'|Vout|'});
ylabel('|Vout| (V)');
title('Boxplot of |Vout| Values');
grid on;

%% Helper Function simulate_Sallenkey

% Computes Vout of a Sallen-Key filter circuit
% using transfer function for unity-gain low-pass filter:
% H(s) = 1 / (R1*R2*C1*C2*s^2 + ((R1+R2)*C2 + R1*C1)*s + 1)
% s = j*omega with omega = 2*pi*f, and Vout is H(jω)*Vin_amp
function Vout = simulate_Sallenkey(R1, R2, C1, C2, Vin_amp, f)
    omega = 2*pi*f; % Angular frequency
    s = 1i * omega; % Complex frequency variable
    
    % Compute denominator of the transfer function
    denom = R1*R2*C1*C2 * s^2 + ((R1 + R2)*C2 + R1*C1) * s + 1;
    H = 1 / denom; % Filter transfer function
    
    % Calculate circuit output voltage
    Vout = H * Vin_amp;
end
