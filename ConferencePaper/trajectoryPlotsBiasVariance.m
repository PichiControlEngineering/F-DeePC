%% Generate two informative plots with the same RMSE
% NOTE: THIS VERSION WAS EVENTUALLY UNUSED IN THE PAPER PUBLISHED ON ARXIV 
% Date: 20-02-2025
% Author: Gert Vankan
% This script creates an illustrative plot to display the difference
% between structural error and noise error
close all

T = 1;         %length of training-data signal
Ts = 0.01;
x_init_data_gen = [Ts];
t_plot = 0:Ts:T-Ts;

A = [-0.9];
A2 = [-0.6];

B = [0.0609];
C = [0];
D = 0;

sys_disc = ss(A,B,C,D,Ts);
u_e = randn(T,1);
u_sin = @(x,u,t) 100*sin((2*pi/100)*t);
[x_d, ~] = systemDataGen(A, B, C, D, Ts, T, x_init_data_gen, u_sin);
[x_2, ~] = systemDataGen(A2, B, C, D, Ts, T, x_init_data_gen, u_sin);

%% Generating a third signal, w_1_tilde, with the same rms as w_1 - w_2 
noise_power = rms(x_d - x_2)^2;  % Power of difference
signal_power = rms(x_2)^2;    
SNR_i = mag2db(signal_power/noise_power); % Add noise such that the RMSE of both signals (disturbed and noised) is equal
x_noisy = addNoise(x_d, SNR_i);
%% Plotting standards % % 
% Defining font size for plots
standard_fontsize = 14;
set(0,'defaulttextinterpreter','latex')
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
%setting font sizes
set(0, 'DefaultAxesFontSize', standard_fontsize);
set(0, 'DefaultLegendFontSize', standard_fontsize*3);

%% plotting trajectories - w_1 and w_2
fig1 = figure(1);
plot(t_plot, x_d,'black', LineWidth=2); grid on; hold on
legend('$y_0$', "FontSize", 16)

fig_structuralBias = figure(2);
plot(t_plot, x_d,'black', LineWidth=2); grid on; hold on
plot(t_plot, x_2,'b', LineWidth=2)

%Creating filled space
filled_x = [t_plot, fliplr(t_plot)];
filled_space = [x_d, fliplr(x_2)];
fill(filled_x, filled_space, 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

%settings
xlim([0, T-Ts])
% fig.Position = [100, 100, 1200, 400]+100;
legend('$y_0$', '$y_1$',"", "FontSize", 16)

% plotting trajectories - w_1 and w_1_tilde
% subplot(1,2,2)
fig_noiseVariance = figure(3);
plot(t_plot, x_d, 'black', LineWidth=2); grid on; hold on
plot(t_plot, x_noisy, 'r', LineWidth=2)

%Creating filled space
filled_x = [t_plot, fliplr(t_plot)];
filled_space = [x_d, fliplr(x_noisy)];
fill(filled_x, filled_space, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

%settings
xlim([0, T-Ts])
legend('$y_0$', '$\tilde{y}_0$',"", "FontSize", 16)
fig_structuralBias.Position = [100, 100, 400, 300]
fig_noiseVariance.Position = [100, 100, 400, 300]
exportgraphics(fig_structuralBias,strcat(pwd,"\Figures\ExampleTrajectories\structuralBias.pdf"),'ContentType','vector')
exportgraphics(fig_noiseVariance,strcat(pwd,"\Figures\ExampleTrajectories\noiseVariance.pdf"),'ContentType','vector')

% sgtitle("The noisy signal $\tilde{y}_0$, and the 'federated' signal $y_2$", "Interpreter", "Latex")