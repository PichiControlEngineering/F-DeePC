%% Test out different types of systems in the F-DPC framework
% NOTE: THIS VERSION WAS EVENTUALLY UNUSED IN THE PAPER PUBLISHED ON ARXIV 
% Version: 19-01-2024
% Author: Gert Vankan
% We compare between a noise-disturbed hankel matrix (like regular), and a
% system-difference based distrubed hankel matrix (like federated). We
% compare the trajectories and parameters in these plots. Waterfall plots
% used in report
% % 
% % % TURN THIS ON IN YOU WANT TO WAIT A LONG TIME % % %

% % % PARAMETERS % % %

% We compare: 
% 1.  regular -> Optimal
% 2.  federated
clear
beta = 0; 
k_spread = 0.03; 
SNR = [];

%plotting and data
T_sim = 30; %time that the system generates data and plots
x_init = [10; -1]; %initial position of the system in simulation

% DeePC parameters
T_ini = 2;      %Initialization time
N = 2;          %prediction horizon
T = 30;         %length of training-data signal

% optimization function for DeePC - Cost function terms
Q = 1;
R = 1;
lambda_y = []; %'[]' turns slack variable off

% noise
SNR_1 = SNR;        %SNR on collected data - set to '[]' to turn off
SNR_2 = SNR;        %SNR for the blue line (x_2) - M=1
sigma_e = 0;       %Noise on the measured output

% % System settings: assumes mass-spring system (can easily be changed by altering system_dataGen2.m)
d = 0;                        %damping on the system (both data-gen and simulation)
Ts = 1;                       %sampling rate of the system: DONT CHANGE
x_init_data_gen = [5;-2]; %random initial start for data generating step

% The different simulated system from the "R_EDDPC" paper
A = [0.7326, -0.0891; 0.1722, 0.9909];
B = [0.0609; 0.0064];
C = [0,1];
D = 0;
sys_disc = ss(A,B,C,D,Ts);
    
%% REFERENCE SIGNALS
 %set a reference signal - Must be a function of t - below some examples
r = @(t) 0*t + 0;
% r = @(t) 2*sin(2*pi*t/(20)) ; %set a reference signal - Must be a function of t


%%setting input function for data-generation: make sure P.O.E. is sufficient
rng(1,'twister')
u_e = randn(T,1);
input_func_noise = @(x,t) u_e(t);


%% HANKEL MATRICES
% % Generate data for the federated hankel matrices
width = T - N - T_ini + 1; %Set width of the hankel matrices

% Initialize Hankel for oracle
[x_d, u_d] = system_dataGen(A,B,C,D, [], T, x_init_data_gen, input_func_noise);
y_d = C*x_d;
w_base = [u_d; y_d];

[U_p_noiseless, Y_p_noiseless, U_f_noiseless, Y_f_noiseless] = createHankel(u_d, (C*x_d), T_ini, N);

% % (2) Generate data for the non-federated, noiseless situation 
H_noiseless = [U_p_noiseless; Y_p_noiseless; U_f_noiseless; Y_f_noiseless];


% Add a single disturbance to A_i

A_i = A + [1 0; 0 1].*(k_spread);

if abs(eig(A_i)) > 1
    disp(eig(A_i))
    error('instable system')
end

[x_d_i, u_d_i] = system_dataGen(A_i,B,C,D, [], T, x_init_data_gen, input_func_noise);
y_d_fed = C*x_d_i;
w_traj_federated = [u_d_i; y_d_fed];

[U_p_fed, Y_p_fed, U_f_fed, Y_f_fed] = createHankel(u_d_i, y_d_fed, T_ini, N);
H_federated = [U_p_fed; Y_p_fed; U_f_fed; Y_f_fed];
figure()
plot(x_d_i(2,:)); grid on; hold on
plot(x_d(2,:)); grid on; hold on

%% Creating a second hankel matrix which is only noise, no structural error due to system differences
MSE_fed = rms(y_d_fed - y_d)^2;  % Power of difference
signal_power = rms(y_d)^2;    
SNR_i = mag2db(signal_power/MSE_fed); % Add noise such that the RMSE of both signals (disturbed and noised) is equal
y_d_noisy = addNoise(y_d, SNR_i);
w_traj_noise_only = [u_d_i; y_d_noisy];

MSE_reg = rms(y_d_noisy-y_d).^2;

[U_p2, Y_p2, U_f2, Y_f2] = createHankel(u_d_i, y_d_noisy, T_ini, N);
H_noisy = [U_p2; Y_p2; U_f2; Y_f2];

    
disp(['RMS_regular = ', num2str(sqrt(MSE_reg)), '; RMS_fed = ', num2str(sqrt(MSE_fed))])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_traj = zeros(6,T_sim);
u_traj = zeros(3,T_sim);
y_pred_traj = zeros(3,T_sim);

%% Plotting an example trajectory
lambda_g_plot = 1e-1;
regular_dataMats = createMatrices(H_noisy,T_ini,N,Q,R,width);
regular_controller_func = @(x, u, t) DeePC_controller(sys_disc.C*x, u, t, T_ini, N, regular_dataMats, r, 1, lambda_y);
[x_traj(1:2,:), u_traj(1,:), y_pred_traj(1,:), g_reg, ~, ~] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, sys_disc.C, Ts, sigma_e, T_ini, regular_controller_func);

federated_dataMats = createMatrices(H_federated,T_ini,N,Q,R,width);
federated_controller_func = @(x, u, t) DeePC_controller(sys_disc.C*x, u, t, T_ini, N, federated_dataMats, r, 1, lambda_y);
[x_traj(3:4,:), u_traj(2,:), y_pred_traj(2,:), g_fed, ~, ~] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, sys_disc.C, Ts, sigma_e, T_ini, federated_controller_func);

OptimalDataMats = createMatrices(H_noiseless,T_ini,N,Q,R,width);
oracle_ctrl = @(x, u, t) DeePC_controller(sys_disc.C*x, u, t, T_ini, N, OptimalDataMats, r, 0,[]);
[x_traj(5:6,:), u_traj(3,:), y_pred_traj(3,:), g_opt,~] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, sys_disc.C, Ts, sigma_e, T_ini, oracle_ctrl);

%plotting result
tit = ["Example system for \lambda_g = ", num2str(lambda_g_plot)];
leg = {"regular", "", "federated", "","oracle",""}; 
t = 0:Ts:T_sim-Ts;

plotDataCell = cell(3,3);
for i1 = 1:size(plotDataCell,2)
    plotDataCell{1,i1} = x_traj((2*i1)-1:2*i1,:); 
    plotDataCell{2,i1} = u_traj(i1,:); 
    plotDataCell{3,i1} = y_pred_traj(i1,:); 
end
fig = PlotFunction(t,plotDataCell,T_ini,T_sim,r,[],leg);

%% looking at g differences
z_lims = 1.1.*[0,max(abs(cell2mat(g_fed(1,T_ini+1:end))),[],"all")]; 
[XX, YY] = meshgrid((T_ini:Ts:T_sim-Ts), 1:1:T-T_ini-N+1);

g_norm_plot = figure();

g_norm_plot.Position = [150, 100, 1450, 600];
subplot(1,2,1)
waterfall(YY, XX, abs(cell2mat(g_reg(1,T_ini+1:end))))
ylabel("T (s)"); xlabel("i");
title('Regularized')
view(55+90, 30)
zlim(z_lims)
% figure()
subplot(1,2,2)
waterfall(YY, XX, abs(cell2mat(g_fed(1,T_ini+1:end))))
ylabel("T (s)"); xlabel("i");
zlim(z_lims)
sgtitle('structure of the parameter vector $|g(k)|$', "interpreter", "latex")
title('Federated')
view(55+90, 30)

%% Functions
function matrixStruct = createMatrices(H,T_ini,N,Q,R,width)
%This function creates the optimization parameters used in the DeePC_controller function
%Mainly used for clarity
% Extract submatrices from H_sum
nu = 1; ny = 1; %we assume a SISO setting
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

function [u, y, g, quadprog_optval,tabrow] = DeePC_controller(y, u, t, T_ini, N, H, r, lambda_g, ~)
% The DeePC controller employs an online optimization based on data from a
% "data-gathering" step. 
% To be used as an implicit function in the 
    f_pre_r = H.f;
    H_p = H.H_p;
    U_f = H.U_f;
    Y_f = H.Y_f;
    R = H.R;
    Q = H.Q;
    H = H.H;
    
    l = size(U_f, 2);
    if t < T_ini
        % Before DeePC is triggered
        u = 0;
        y = 0;
        g=zeros(l,1);
        quadprog_optval = NaN;
        tabrow = 0;
    else
        % disp(['t = ', num2str(t)])
        %Constructing r matrix
        d_t = 0:1:N-1;
        r_mat = r(t+d_t)';

        f = (r_mat'*f_pre_r)';
        
        % Implement a special procedure on t = T_ini, where we append 1
        % zero to the l_ini vector
        if t == T_ini
            y_ini = [0,y(:, t-T_ini+1:t-1)];
            u_ini = [0,u(:, t-T_ini+1:t-1)];
        else
            %regular procedure
            y_ini = y(:, t-T_ini:t-1);
            u_ini = u(:, t-T_ini:t-1);
        end
        l_ini = [u_ini'; y_ini'];
        W = 2*(H + eye(size(H)).*lambda_g);
        
        % Start optimalization
        quadprogopts = optimoptions("quadprog", "MaxIterations", 2.5e2,"ConstraintTolerance", 1e-7, 'Display','off');
        [g,quadprog_optval,exitflag] = quadprog(W,f,[],[],H_p,l_ini,[],[],[],quadprogopts);
        
        if exitflag ~= 1
            % g = W \ (H_p'*inv(H_p * inv(W) * H_p')*(l_ini + H_p / W * f)-f);
            % warning('Innacurate problem, impossible to draw conclusions')
            quadprog_optval = NaN;
        end

        % Predict control input and output
        u_predicted = U_f * g;
        u = u_predicted(1);

        y_predicted = Y_f * g;
        y = y_predicted(1);
        
        % OPTIONAL: Report on analytic-theoretical difference
        % optval_other = R.*norm(u_predicted).^2 + Q.*norm(y_predicted - r_mat).^2 + lambda_g*norm(g).^2 - r_mat'*Q*r_mat;
        % difference = quadprog_optval - optval_other;
        % disp(['difference between c_minimizer and c_analytic: ', num2str(difference)])

        tabrow = [R.*norm(u_predicted).^2, Q.*norm(y_predicted - r_mat).^2, lambda_g*norm(g).^2, quadprog_optval];
    end
end 


function [lambda_g_opt,u_bias, y_bias, u_rms, y_rms, optval_diff_grid] = regularizationTuning(lambda_g_range, x_optimal, u_optimal, T_ini, T_sim, N, H_list, sys_disc, x_init,Q,R) %this function looks at the initial state to determine an optimal \lambda_g tuning
    %We take 2 samples taken from the same data-set, and compare the cost
    %function values. Do they converge at some value of lambda_g?
    u_bias = zeros(length(H_list),length(lambda_g_range));
    u_rms = zeros(length(H_list),length(lambda_g_range));
    y_bias = zeros(length(H_list),length(lambda_g_range));
    y_rms = zeros(length(H_list),length(lambda_g_range));

    optval_diff_grid = zeros(length(H_list),length(lambda_g_range));

    for j = 1:length(H_list)
        for i = 1:length(lambda_g_range)

            deePC_ctrl = @(x, u, t) DeePC_controller(sys_disc.C*x, u, t, T_ini, N, H_list{j}, @(t) 0*t, lambda_g_range(i), []);
            [x_traj_i, u_traj_i, ~, ~, optval_difference] = stepStates(x_init, T_sim, sys_disc.A, sys_disc.B, sys_disc.C, 1, 0, T_ini, deePC_ctrl);
            u_bias(j,i) = rmse(u_traj_i, u_optimal);
            u_rms(j,i) = rms(u_traj_i);

            y_bias(j,i) = rmse(sys_disc.C*x_traj_i, sys_disc.C*x_optimal);
            y_rms(j,i) = rms(sys_disc.C*x_traj_i);
            optval_diff_grid(j,i) = norm(optval_difference);
        end
    end
    cost_matrix = u_bias.^2;
    [~,lambda_index_grid] = min(cost_matrix,[],2);
    for j = 1:length(H_list)
        lambda_g_opt(j) = lambda_g_range(lambda_index_grid(j));
    end
end

function fig = PlotFunction(t,PlotDataCell,T_ini,T_sim,r, tit, leg)
    %This function plots results from a simulation step
    % To be given in the PlotDataCell format
    fig = figure();
    Ts = t(2)-t(1);
    
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
        y_pred_traj = PlotDataCell{3,i};
    
        % % plotting x_1
        subplot(3,1,1)
        stairs(t, x_traj(1,:), 'Color', col,"LineWidth",linesize); 
        hold on; grid on
        xline(T_ini, 'black:','LineWidth', 1, 'Label','');
        ylabel('x_1')
        
        ylim([-1 1] + [x1_min, x1_max]);
        xlim([0, T_sim])
        
        % % x_2
        subplot(3,1,2)
        stairs(t, x_traj(2,:), 'Color', col,"LineWidth",linesize); 
        hold on; grid on;
        % stairs(t+Ts, y_pred_traj(1:end), 'o', 'Color', col); % we are PREDICTING the next output, so we should move all y_
        xline(T_ini, 'r:', 'LineWidth', 1, 'Label','');
        
        ylim([-1 1] + [x2_min, x2_max]);
        xlim([0, T_sim])
        
        xlabel('t (s)')
        ylabel('x_2')
        
        subplot(3,1,3)
        stairs(t, u_traj, 'Color', col,"LineWidth",linesize); 
        hold on; grid on;
        xline(T_ini, 'r:', 'LineWidth', 1, 'Label','');
        
        ylim(0.2*[-1 1] + [u_min,u_max]);
        xlim([0, T_sim])
        
        xlabel('t (s)')
        ylabel('u')
    end

    if exist("r") && ~all(r(t) == 0)
        subplot(3,1,2)
        stairs(t, r(t),'black--', LineWidth=0.75)
    end

    if exist("tit")
        sgtitle(tit);
    end

    %add legend to all subplots
    if exist("leg")
    for i = 1:3
        subplot(3,1,i)
        legend(leg, 'location', 'best')
    end
    end
end

%simulate
function [x_traj, u_traj,y_pred_traj,g_cell,cvx_optval,tabrow] = stepStates(x_init, T, A_disc, B_disc, C_disc, Ts, sigma_e, T_ini, input_func)
    %this function works to simulate any controller 'input_func' on the
    %discrete system defined by the state_space matrices A_disc, B_disc
    x_traj = zeros(size(x_init, 1), T) ; % Saving past states for plotting
    u_traj = zeros(size(B_disc, 2), T);
    y_pred_traj = zeros(1, T); % Add this line to store predicted outputs
    x = x_init;
    g_cell = cell(1, T/Ts);
    cvx_optval = zeros(1,T/Ts);
    for i = 1:Ts:T
        % Save current state
        x_traj(:, i) = x; 

        % OPTIONAL: Add noise
        if i == T_ini+1
            [u, y_pred, g, cvx_optval_i,tabrow] = input_func(x_traj+randn(size(x_init))*sigma_e, u_traj, i);
        else
            [u, y_pred, g, cvx_optval_i,~] = input_func(x_traj+randn(size(x_init))*sigma_e, u_traj, i);
        end
        y_pred_traj(:, i) = y_pred;
        

        % Update state
        x = A_disc * x + B_disc .* u;
        u_traj(:, i) = u;
        g_cell{1, i} = g;
        cvx_optval(1,i) = cvx_optval_i;
    end
end

function [u_LQR,y_pred,g,cvx_optval] = LQR_control(x, u, t, r, K, T_ini)
    %LQR controller
    if t <= T_ini
        % Before DeePC is triggered
        u_LQR = 0;
    else
        e = r(t) - x(:, t);
        u_LQR = K*e;
    end
    y_pred = NaN;
    g = NaN;
    cvx_optval=NaN;
    i = NaN;
end