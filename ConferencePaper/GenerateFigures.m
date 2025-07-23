%% Test out different types of systems in the F-DPC framework
% Version: 30-03-2025
% Author: Gert Vankan
% We generate all figures which are used in the conference paper
clear
close

%Options
regenerate_data = 0;     % This option restarts monte carlo simulation for the parameters below
seperate_subfigures = 1; % this option seperates and saves the subfigures individually;
filter_inactive_controllers = 0; % The controller have a tendency to be over-regularized and not do anything, this option filters those controllers from the plots

% Plot settings (general)
standard_fontsize = 14; %set standard fontsize
% figure_directory = "\Figures\OldSysOldPert_regularization";
% figure_directory = "\Figures\OldSysDiagPert";
figure_directory = "\Figures\OldSysOldPert";

%%% EXISTING DATASETS %%%
% filePath = "\DataSets\OldSystemOldPert_SlackVariable.mat";
% filePath = "\DataSets\OldSystemNewPert.mat";
filePath = "\DataSets\OldSystemOldPert.mat";      
% filePath = "\DataSets\DifferentSystem.mat";

% Where files are saved
current_location = convertCharsToStrings(pwd);

% % % PARAMETERS % % % - these only have effect when setting
K = 150;
N_gridsize = 15; %amount of values of lambda_g in the logspace
binNo = 25;      % number of bins for the histogram

beta = 0.1; 
k_spread = 0.05;
M = 55; 
M2 = 21; 
M3 = 11;
M_list = [M,M2,M3];

SNR = 20;

%plotting and data
T_sim = 50; %time that the system generates data and plots
x_init = [1; 1]; %initial position of the system in simulation

% DeePC parameters
T_ini = 3;      %Initialization time
N = 3;          %prediction horizon
T = 50;         %length of training-data signal

% optimization function for DeePC - Cost function terms
Q = 1;
R = 0.01;
lambda_y = []; % '[]' turns slack variable off

% noise
sigma_e = 0;       %Noise on the measured output

% System settings
Ts = 1;                       %sampling rate of the system: DONT CHANGE
x_init_data_gen = [-1; 1]; %random initial start for data generating step
t = 0:Ts:T_sim-Ts;

% SS matrices from the "R_EDDPC" paper (Breschi, 2023)
A = [0.7326, -0.0891; 0.1722, 0.9909];
B = [0.0609; 0.0064];
C = [0,1];
D = 0;

% SS matrices in the "DifferentSystem.mat
% A = [0, 0.9; -1, 0.9];
% B = [1.1; 1];
% C = [0,1];
% D = 0;
sys_disc = ss(A,B,C,D,Ts);

num_controllers = 1+length(M_list); %we have 1 regular trajectory, and 3 federated ones
%% REFERENCE SIGNALS
 %set a reference signal - Must be a function of t - below some examples
r = @(t) 0*t + 0;
% r = @(t) 2*sin(2*pi*t/(20)) ; %set a reference signal - Must be a function of t

%% script prerequisites 
filePath_Figures = current_location + figure_directory;
filePath_DataSets = current_location + filePath;

if regenerate_data == 1
    u_diff_grid = zeros(num_controllers,N_gridsize,K);
    y_diff_grid = zeros(num_controllers,N_gridsize,K);
    u_rms_grid = zeros(num_controllers,N_gridsize,K);
    y_rms_grid = zeros(num_controllers,N_gridsize,K);
    y_rms_optimal_grid = zeros(K,1);
    u_rms_optimal_grid = zeros(K,1);

    lambda_g_range = logspace(-3,3,N_gridsize);

    %setting input function for data-generation: make sure P.O.E. is sufficient
    disp(['regenerating the data for the experiments, K = ', num2str(K)])
    close all

    for i_randn = 1:K
        %regenerate a new noise iteration each time
        rng('shuffle')
        u_e = randn(T,1);
        input_func_noise = @(x,u,t) u_e(t);

        if i_randn == 1
            tic()
        else
            disp(['Monte Carlo run ', num2str(i_randn), '/', num2str(K)])
        end
        %% HANKEL MATRICES
        % % Generate data for the federated hankel matrices
       
        % Initialize Hankel for oracle
        [x_d, u_d] = systemDataGen(A,B,C,D, [], T, x_init_data_gen, input_func_noise);
        y_d = C*x_d;
        y_d_noisy = addNoise(y_d, SNR);
        w_base = [u_d; y_d];
    
        [U_p_noiseless, Y_p_noiseless, U_f_noiseless, Y_f_noiseless] = createHankel(u_d, y_d, T_ini, N);
        [U_p_noisy, Y_p_noisy, U_f_noisy, Y_f_noisy] = createHankel(u_d, y_d_noisy, T_ini, N);
    
        % % (2) Generate data for the non-federated, noiseless situation 
        H_noiseless = [U_p_noiseless; Y_p_noiseless; U_f_noiseless; Y_f_noiseless];
        H_noisy = [U_p_noisy; Y_p_noisy; U_f_noisy; Y_f_noisy];
    
        % Pre-allocate matrices for the Hankelmats and weightings alpha
        H_fed_cell = cell(size(M_list));
        federated_dataMats = cell(size(M_list));
        for i = 1:length(M_list)
            H_fed_cell{i} = generateFederatedHankel(M_list(i),k_spread,beta,H_noisy,A,B,C,T,T_ini,N,x_init_data_gen,input_func_noise,SNR);
            federated_dataMats{i} = createMatrices(H_fed_cell{i},T_ini,N,Q,R);
        end
        
        %% SIMULATING SYSTEMS WITH DEFINED CONTROLLER
        % The order is: 1. Regular (noisy), 
        %               2. Federated matrices
        %               3. oracle controller (noiseless)
        regular_dataMats = createMatrices(H_noisy,T_ini,N,Q,R);
        H_lists = [{regular_dataMats}, federated_dataMats]; %concatenating the Hankel matrices into one big list
        
        OptimalDataMats = createMatrices(H_noiseless,T_ini,N,Q,R);
        oracle_ctrl = @(x, u, t) DeePCcontroller(sys_disc.C*x, u, t, T_ini, N, OptimalDataMats, r, 0,[]);
        [x_traj_oracle, u_traj_oracle,y_pred_oracle,g_opt,~] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, Ts, sigma_e, oracle_ctrl);
        
        %setting range to test lambda
        [lambda_g_opts, u_diff_grid_i, y_diff_grid_i, u_rmse_i, y_rmse_i, optval_diff_grid] = regularizationTuning(lambda_g_range, x_traj_oracle, u_traj_oracle, T_ini, T_sim, N, H_lists, sys_disc, x_init);

        % distance from oracle trajectory (in RMSE sense)
        y_diff_grid(:,:,i_randn) = y_diff_grid_i;
        u_diff_grid(:,:,i_randn) = u_diff_grid_i;
    
        %absolute RMS of the F-DeePC & Dee-PC controller
        y_rms_grid(:,:,i_randn) = y_rmse_i;
        u_rms_grid(:,:,i_randn) = u_rmse_i;
    
        %RMSE & RMSU of the oracle trajectory
        y_rms_optimal_grid(i_randn) = rms(x_traj_oracle(2,:));
        u_rms_optimal_grid(i_randn) = rms(u_traj_oracle);
    
        %%
        if i_randn == 1
            runtime = toc();
            disp(['Approximate Monte-Carlo sim time is: ', num2str(runtime*K, 3), ' s'])
            pause(1)
        end
    end

    %saving nescessary system data
    save(filePath_DataSets)
else
    close all
    disp('Loading file:')

    workspace_data = load(filePath_DataSets); % Load everything into a struct
    workspace_data = rmfield(workspace_data, {'seperate_subfigures', 'filter_inactive_controllers','filePath_Figures', 'current_location', 'filePath_DataSets'}); % Remove the unwanted variable
    fn = fieldnames(workspace_data); % Get remaining variable names
    % Assign each variable to the workspace
    for i = 1:length(fn)
        assignin('base', fn{i}, workspace_data.(fn{i}));
    end
    regenerate_data = 0;

    % Loading oracle controller
    OptimalDataMats = createMatrices(H_noiseless,T_ini,N,Q,R);
    oracle_ctrl = @(x, u, t) DeePCcontroller(sys_disc.C*x, u, t, T_ini, N, OptimalDataMats, r, 0,[]);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining font size for plots
set(0,'defaulttextinterpreter','latex')
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
%setting font sizes
set(0, 'DefaultAxesFontSize', standard_fontsize);
set(0, 'DefaultLegendFontSize', standard_fontsize*3);

%%
%we also want to check on 2 special conditions: 
% 1. Instable controller, 2. Controller not doing anything 
if filter_inactive_controllers
    u_grid_inactive = u_rms_grid<=rms(u_traj_oracle)/20; %we say our controller is "inactive" if it does not have a fraction of the input energy of the oracle traj 
    y_diff_grid(u_grid_inactive) = NaN;
    u_diff_grid(u_grid_inactive) = NaN;
end
% averaging monte carlo results
y_diff_grid_m = nanmean(y_diff_grid(:,:,:), 3);
u_diff_grid_m = nanmean(u_diff_grid(:,:,:), 3);
y_rms_grid_m = nanmean(y_rms_grid(:,:,:), 3);
u_rms_grid_m = nanmean(u_rms_grid(:,:,:),  3);

% Taking the median 
y_diff_grid_median = nanmedian(y_diff_grid(:,:,:), 3);
u_diff_grid_median = nanmedian(u_diff_grid(:,:,:), 3);
% the Lambda_g's for which the median is minimal
lambda_g_min_median = lambda_g_range*(nanmedian(y_diff_grid(:,:,:),3) == min(nanmedian(y_diff_grid(:,:,:),3),[],2))';

% and Standard dev
y_diff_grid_std = std(y_diff_grid(:,:,:),[],3);
u_diff_grid_std = std(u_diff_grid(:,:,:), [],3);

%computing mean RMS for RMS comparison plot
y_rms_optimal_m = nanmean(y_rms_optimal_grid);
u_rms_optimal_m = nanmean(u_rms_optimal_grid);


%% Looking at histograms over Monte carlo sim creating Histogram data for every monte-carlo iteration

num_histograms = size(y_diff_grid, 2); % Number of rows in y_analysis

%we determine the size of the bin edges between -0.1 and the upper 95th
%prctile (to prevent shaping the histogram around extreme outliers)

%Initialize bin edges
bin_edges_y = linspace(-0.005, prctile(y_diff_grid,95,"all"), binNo+1);
bin_edges_y_RMS  = linspace(0, prctile(y_rms_grid,95,"all"), binNo+1);
bin_edges_u = linspace(-0.005, prctile(u_diff_grid,95,"all"), binNo+1);

mean_edge_val_y = (bin_edges_y(1:end-1) + bin_edges_y(2:end)) / 2;
mean_edge_val_y_RMS = (bin_edges_y(1:end-1) + bin_edges_y(2:end)) / 2;
mean_edge_val_u = (bin_edges_u(1:end-1) + bin_edges_u(2:end)) / 2;
mean_edge_val_u(1) = 0; mean_edge_val_y(1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We generate an num_histograms x nBins cell containing the histogram data,
%and save these in a cell for each seperate controller
hist_vals_y_cell = cell(num_controllers,1);
hist_vals_y_RMS_cell = cell(num_controllers,1);
hist_vals_u_cell = cell(num_controllers,1);

for i_controller_type = 1:num_controllers
    y_analysis = squeeze(y_diff_grid(i_controller_type,:,:)); 
    y_analysis_RMS = squeeze(y_rms_grid(i_controller_type,:,:));
    u_analysis = squeeze(u_diff_grid(i_controller_type,:,:));
    
    % Pre-allocate
    hist_vals_y = zeros(num_histograms, binNo);
    hist_vals_u = zeros(num_histograms, binNo);
    hist_vals_y_RMS = zeros(num_histograms, binNo);

    % 'FOR' loop to get the histogram values for each lambda_g value
    for i=1:num_histograms
        figure(93);
        h_y = histogram(y_analysis(i,:), bin_edges_y);
        figure(94)
        h_u = histogram(u_analysis(i,:), bin_edges_u);
        figure(95);
        h_y_rms = histogram(y_analysis_RMS(i,:), bin_edges_y_RMS);

        hist_vals_y(i,:) = [h_y.Values];
        hist_vals_u(i,:) = [h_u.Values];
        hist_vals_y_RMS(i,:) = [h_y_rms.Values];
    end

    hist_vals_y_cell{i_controller_type} = hist_vals_y;
    hist_vals_y_RMS_cell{i_controller_type} = hist_vals_y_RMS;
    hist_vals_u_cell{i_controller_type} = hist_vals_u;
    close(93); close(94); close(95)
end
%% Generating a 3d-Histogram
%set titles
hist_title_y = ["${RMSE}_y$", "${RMSE}_y^{fed}$"];
hist_title_u = ["${RMSE}_u$", "${RMSE}_u^{fed}$"];
plot_sizes = [200, 200, 500, 400];

for i_plot_controller = [1,2]
    % Prepare a 3D bar plot
    lambda_log_spaced = log10(lambda_g_range);
    [X_mesh_y, Y_mesh_y] = meshgrid(lambda_log_spaced, mean_edge_val_y);
    [X_mesh_u, Y_mesh_u] = meshgrid(lambda_log_spaced, mean_edge_val_u);
   
    % First plot showing RMSE_y
    if seperate_subfigures
        fig_y = figure(9+i_plot_controller);
        fig_y.Position = plot_sizes;
    else
        histogram_3d = figure(10);
        histogram_3d.Position = [200, 100, 1200, 700];
        subplot(2,2,i_plot_controller)
    end

    %plot histogram
    s = surf(X_mesh_y, Y_mesh_y, hist_vals_y_cell{i_plot_controller}','HandleVisibility','off'); hold on; 
    %plot mean line
    plot3(lambda_log_spaced, y_diff_grid_median(i_plot_controller,:),ones(1,num_histograms).*K+1,'red--', 'LineWidth',0.8,'DisplayName','$M_y$')
    plot3(lambda_log_spaced, y_diff_grid_m(i_plot_controller,:),ones(1,num_histograms).*K+1,'black--', 'LineWidth',1,'DisplayName','$\mu_y$')
    
    ylim([0,0.8]) %for the old system comparison plots
    ylabel(hist_title_y(i_plot_controller));
    set(s,'linestyle','none')

    %Second plot displaying RMSE_u
    if seperate_subfigures
        fig_u = figure(11+i_plot_controller);
        fig_u.Position = plot_sizes;
    else
        subplot(2,2,i_plot_controller+2)
    end
    s = surf(X_mesh_u, Y_mesh_u, hist_vals_u_cell{i_plot_controller}','HandleVisibility','off'); hold on
    set(s,'linestyle','none')
    
    % Plotting mean and Median
    plot3(lambda_log_spaced, u_diff_grid_m(i_plot_controller,:),ones(1,num_histograms).*max(hist_vals_u+1,[],'all'),'black--', 'LineWidth',1,'DisplayName','$M_u$')
    plot3(lambda_log_spaced, u_diff_grid_median(i_plot_controller,:),ones(1,num_histograms).*max(hist_vals_u+1,[],'all'),'red--', 'LineWidth',1,'DisplayName','$\mu_u$')
    ylim([0, max(Y_mesh_u,[],"all")])
    ylabel(hist_title_u(i_plot_controller));
end
% Setting plot settings for all plots
for k_2 = 1:4
    view(2); box on;
    if seperate_subfigures
        figure(9+k_2)  
    else
        subplot(2,2,k_2)
    end
    
    if seperate_subfigures && any(k_2 == [2,4])
        cbar = colorbar;
        cbar.Label.String = 'counts';
        a = figure(9+k_2);
        a.Position = plot_sizes + [0, 0, 70, 0];
    end
    cmap = parula(256); % 256 is the number of colors
    cmap(1, :) = [1, 1, 1]; % RGB for white
    colormap(cmap);
    c_map_max = K; %round up the max
    clim([0 c_map_max])
    
    xlim([log10(min(lambda_g_range)), log10(max(lambda_g_range))])
    xticks([-3, -2, -1, 0, 1, 2, 3]); % Set tick locations
    xticklabels({'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'}); % Custom tick labels in LaTeX format
    
    xlabel('$$\lambda_g$$', "Interpreter","latex");
    zlabel('Counts');
    % legend('Location', 'northeast')
end

if seperate_subfigures
    for k_2 = 1:4
    filename_traj_plot = filePath_Figures+"\3d_histogram_"+convertCharsToStrings(num2str(k_2))+".pdf";
    exportgraphics(figure(9+k_2),filename_traj_plot, 'ContentType', 'vector')
    end
end

% Determining optimal Lambda: where on average the u oracle bias (and consequently
% y oracle bias are minimal. DO NOT use lambda_g_opts, since it is flawed
% (still figuring out what's going wrong there)
lambda_g_mean_opts = (u_diff_grid_m==min(u_diff_grid_m,[],2))*lambda_g_range';
%% Plotting absolute errors
figure(3)
for i = 1:K
semilogx(lambda_g_range,y_rms_grid(1,:,i), 'b.'); hold on; grid on
semilogx(lambda_g_range,y_rms_grid(2,:,i), 'r.'); 
end
ylabel('$RMS(y_i)$')
xlabel('$\lambda_g$')
ylim([0,max(y_rms_grid,[],"all")])
close(3) %unused in final report
%% Plotting an example trajectory
regular_controller_func = @(x, u, t) DeePCcontroller(sys_disc.C*x, u, t, T_ini, N, regular_dataMats, r, lambda_g_mean_opts(1), lambda_y);
[x_traj(1:2,:), u_traj(1,:)] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, Ts, sigma_e, regular_controller_func);

federated_controller_func = @(x, u, t) DeePCcontroller(sys_disc.C*x, u, t, T_ini, N, federated_dataMats{2}, r, lambda_g_mean_opts(2), lambda_y);
[x_traj(3:4,:), u_traj(2,:)] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, Ts, sigma_e, federated_controller_func);

[x_traj(5:6,:), u_traj(3,:)] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, Ts, sigma_e, oracle_ctrl);

%plotting results as an example
leg = {['reg'], "", ['fed'], "", "orac", ""}; 

plotDataCell = cell(3,3);
for i1 = 1:3
    plotDataCell{1,i1} = x_traj((2*i1)-1:(2*i1),:); 
    plotDataCell{2,i1} = u_traj(i1,:); 
end

trajectory_plot = PlotFunction(t,plotDataCell,T_ini-1,T_sim,r,[],leg,seperate_subfigures);
if seperate_subfigures
    for i_trajectory_plots = 1:3
        trajectory_plot(i_trajectory_plots).Position = 0.75*[200, 100, 800, 250];
        filename_traj_plot = strcat(filePath_Figures,"\trajectory_",num2str(i_trajectory_plots),".pdf");
        exportgraphics(trajectory_plot(i_trajectory_plots),filename_traj_plot, 'ContentType', 'vector')
    end
else
    trajectory_plot.Position = [200, 100, 800, 700];
end
%% Histogram of the RMSE of the controllers at optimal values for the federated and regular trajectories
legend_text = {['reg, $\lambda_g = ', num2str(lambda_g_mean_opts(1),2),'.0$']};

RMSE_histogram = figure(1);
for i_histogram = 1:2
    hist_vals = hist_vals_y_cell{i_histogram};
    histvals_optimal_lambda = hist_vals(lambda_g_range==lambda_g_mean_opts(i_histogram+1),:);  hold on; grid on; box on;
    histogram("BinEdges",bin_edges_y,"BinCounts",histvals_optimal_lambda);
    
    leg_i = ['fed, $\lambda_g = ', num2str(lambda_g_mean_opts(i_histogram+1),2),'$'];
    legend_text = [legend_text, {leg_i}];
end
legend(legend_text,'FontSize',standard_fontsize)
xlabel('${RMSE}_y, {RMSE}_y^{fed}$')
ylabel('counts')
% save figure
if seperate_subfigures
    filename_histogram = strcat(filePath_Figures,"\histogram_optimal_lambda.pdf");
    RMSE_histogram.Position = 0.7*[500, 250, 600, 500];
    exportgraphics(RMSE_histogram,filename_histogram, 'ContentType', 'vector')
end
%% Histogram of the RMS histogram at optimal values for the federated and regular trajectories
legend_text = {['reg']};

RMS_histogram = figure(3);
for i_histogram = 1:2
    hist_vals = hist_vals_y_RMS_cell{i_histogram};
    histvals_optimal_lambda = hist_vals(lambda_g_range==lambda_g_mean_opts(i_histogram+1),:);  hold on; box on;
    histogram("BinEdges",bin_edges_y,"BinCounts",histvals_optimal_lambda);
    
    % leg_i = ['fed, $\lambda_g = ', num2str(lambda_g_mean_opts(i_histogram+1),2),'$'];
    leg_i = ['fed'];
    legend_text = [legend_text, {leg_i}];
end
legend(legend_text,'FontSize',standard_fontsize,"Location","best")
xlabel('${RMS}_y, {RMS}_y^{fed}$')
ylabel('counts')

% save figure
if seperate_subfigures
    filename_histogram = strcat(filePath_Figures,"\RMS_histogram_optimal_lambda.pdf");
    RMS_histogram.Position = 0.7*[500, 250, 600, 500];
    exportgraphics(figure(3),filename_histogram, 'ContentType', 'vector')
end
%% Box whisker plot
%first formatting all the data for the boxplot
RMS_optimal_boxplot = zeros(size(y_diff_grid,3),num_controllers);
boxplot_labels = cell(1,num_controllers);
boxplot_order = [1,num_controllers:-1:2]; %we have to swap around the order, such that it is [reg, M=low, M=med, M=hi] 
index_counter = 1;
for i_boxplot = boxplot_order
    RMS_optimal_boxplot(:,index_counter) = y_diff_grid(i_boxplot,lambda_g_range==lambda_g_mean_opts(i_boxplot),:);

    %automatically generating labels
    if i_boxplot == 1
        boxplot_labels{index_counter} = 'reg';
    else
        boxplot_labels{index_counter} = ['M = ', num2str(M_list(i_boxplot-1))];
    end
    index_counter = index_counter + 1;
end
%Plotting + settings
bxplot_fig = figure(2);
bxplot_fig.Position = 0.7*[500, 250, 600, 500];
boxplot(RMS_optimal_boxplot,'Labels',boxplot_labels);

ylabel("$RMSE_y, RMSE_y^{fed}$")
ylim([0,1.1*prctile(RMS_optimal_boxplot,100,"all")])

% save figure
if seperate_subfigures
    filename_boxplot = strcat(filePath_Figures,"\boxplot_lambda.pdf");
    exportgraphics(figure(2),filename_boxplot, 'ContentType', 'vector')
end

%%%%%%%%%%%%%%%%
%% Functions %%%
%%%%%%%%%%%%%%%%
function H_federated = generateFederatedHankel(M,k_spread,beta,H0,A,B,C,T,T_ini,N,x_init_data_gen,input_func_noise,SNR)
    alpha = zeros(M, 1);
    H_cell_fed = cell(M,1);
    
    k_federated = linspace(-k_spread,k_spread,M); %define spread around distribution
    if M == 1
        k_federated = 0;
    end

    for i = 1:M
        % A_i = A + ([0, 1; -1, 0]*k_federated(i));
        A_i = A + (eye(2)*k_federated(i));

        if abs(eig(A_i)) > 1
            % disp(eig(A_i)); 
            error('instable system')
        end
        
        [x_d_i, u_d_i] = systemDataGen(A_i,B,C,0, [], T, x_init_data_gen, input_func_noise);
        y_d_fed_i = addNoise(C*x_d_i, SNR);
        [U_p_fed, Y_p_fed, U_f_fed, Y_f_fed] = createHankel(u_d_i, y_d_fed_i, T_ini, N);
        H_cell_fed{i} = [U_p_fed; Y_p_fed; U_f_fed; Y_f_fed];
        
        % Compute the weight alpha(i)
        alpha(i) = exp(-beta * norm(H0 - H_cell_fed{i}, 2));
        
        % Compute the weighted sum of H matrices using vectorized operation
    end
    alpha = alpha./sum(alpha);
    H_federated = sum(cat(3, H_cell_fed{:}) .* reshape(alpha, 1, 1, []), 3);
end

function matrixStruct = createMatrices(H,T_ini,N,Q,R)
    %This function creates the optimization parameters used in the DeePCcontroller function
    %Mainly used for clarity
    % Extract submatrices from H_sum
    nu = 1; ny = 1; %we assume a SISO setting
    width = size(H,2);

    U_p = H(1:(T_ini * nu), 1:width);
    Y_p = H((T_ini * nu + 1):(T_ini * nu + T_ini * ny), 1:width);
    U_f = H((T_ini * nu + T_ini * ny + 1):((T_ini + N) * nu + T_ini * ny), 1:width);
    Y_f = H((T_ini + N) * nu + T_ini * ny + 1:end, 1:width);
    
    %constructing quadratic problem matrices
    Hessian = (Y_f' * Y_f * Q) + (U_f' * U_f * R);
    f_pre_r = -2 * Q * Y_f;
    
    % Equality constraints
    H_p = [U_p; Y_p];
    
    matrixStruct = struct('H',Hessian,'f',f_pre_r,'H_p',H_p,'U_f',U_f,'Y_f',Y_f,'Q',Q,'R',R);
end

function fig = PlotFunction(t,PlotDataCell,T_ini,T_sim,r, tit, leg, seperate_subfigures)
    %This function plots results from a simulation step
    % To be given in the PlotDataCell format
    if seperate_subfigures
        fig = [figure(20),figure(21),figure(22)];
    else
        fig = figure(20);
    end
    % Defining Colors: (blue, red, yellow, ?orange?)
    colorlist = {'blue','r',"#77AC30",'#BF40BF'};

    %determine maxima:
    fullmat = cell2mat(PlotDataCell);
    x1_max = max(fullmat(1,:),[],"all"); x1_min = -max(-fullmat(1,:),[],"all");
    x2_max = max(fullmat(2,:),[],"all"); x2_min = -max(-fullmat(2,:),[],"all");
    u_max = max(fullmat(3,:),[],"all"); u_min = -max(-fullmat(3,:),[],"all");

    % Plot
    for i = 1:size(PlotDataCell,2)
        col = colorlist{i};
        linesize = 0.8;
    
        x_traj = PlotDataCell{1,i};
        u_traj = PlotDataCell{2,i};
    
        %%%%%%%%%%% x_1 PLOT %%%%%%%%%%%%%%
        if seperate_subfigures
            figure(20)
        else
        subplot(3,1,1)
        end
        stairs(t, x_traj(1,:), 'Color', col,"LineWidth",linesize); 
        hold on; grid on; box on;
        xline(T_ini, 'black:','LineWidth', 1, 'Label','');
        ylabel('$x_1$')
        
        ylim([-3 1] + [x1_min, x1_max]);
        xlim([0, T_sim])
        xlabel('$t (s)$')

        %%%%%%%%%% x_2 PLOT %%%%%%%%%%%%%%
        if seperate_subfigures
            figure(21)
        else
            subplot(3,1,2)
        end
        stairs(t, x_traj(2,:), 'Color', col,"LineWidth",linesize); 
        hold on; grid on; box on;;
        % stairs(t+Ts, y_pred_traj(1:end), 'o', 'Color', col); % we are PREDICTING the next output, so we should move all y_
        xline(T_ini, 'r:', 'LineWidth', 1, 'Label','');
        
        ylim([-3 1] + [x2_min, x2_max]);
        xlim([0, T_sim])
        
        xlabel('$t (s)$')
        ylabel('$x_2$')
        
        %%%%%%%%%%%%%% INPUTS PLOT %%%%%%%%%%
        if seperate_subfigures
            figure(22)
        else
            subplot(3,1,3)
        end
        stairs(t, u_traj, 'Color', col,"LineWidth",linesize); 
        hold on; grid on; box on;;
        xline(T_ini, 'r:', 'LineWidth', 1, 'Label','');
        
        ylim([-6.5,1]);
        xlim([0, T_sim])
        
        xlabel('$t (s)$')
        ylabel('$u$')
    end

    if exist("r","var") && ~all(r(t) == 0)
        if seperate_subfigures
            figure(22)
        else
            subplot(3,1,2)
        end
        stairs(t, r(t),'black--', LineWidth=0.75)
    end

    if exist("tit","var")
        sgtitle(tit);
    end

    %add legend to all subplots
    if exist("leg","var")
        for i = 1:3
            if seperate_subfigures
                figure(19+i)
            else
                subplot(3,1,i)
            end
            legend(leg, 'location', 'southeast')
        end
    end
end