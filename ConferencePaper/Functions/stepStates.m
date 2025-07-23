function [x_traj, u_traj, varargout] = stepStates(x_init, T, A_disc, B_disc, Ts, ~, input_func, varargin)
% This function works to simulate any controller 'input_func' on the discrete
% system defined by a state-space model
% Inputs:
% - x_init:             Initial state
% - T:                  Total simulation time
% - A_disc, B_disc:     Discrete-time state-space matrices
% - Ts:                 Sampling time
% - sigma_e:            Measurement noise standard deviation
% - input_func(x, u, t):The controller function handle
% - varargin:
%       "Timer" - prints time step outputs
%       "holdState" - fixes the state for varargin{2} steps, then releases
%       state

% Outputs:
% - x_traj:             State trajectory
% - u_traj:             Input trajectory
% - varargout:
%   - g_cell:           Cell of g values for each time step
%   - cvx_optval:       Optimization cost at each time step
%   - alpha_cell:       alpha selected by the controller at each time step

% Initialization
steps = round(T / Ts);
x_traj = zeros(size(x_init, 1), steps);
u_traj = zeros(size(B_disc, 2), steps);
computation_times = zeros(1, steps);
n_outputs = nargout(input_func);

if nargout > 2
    g_cell = cell(1, steps);
    cvx_optval = zeros(1, steps);
    if nargout > 3
        alpha_cell = cell(1, steps);
    end
end

timer_on_flag = any(strcmpi(varargin, "timer"));
holdState_flag = any(strcmpi(varargin, "holdState"));

x = x_init;

for i = 1:steps
    t = i * Ts;

    % Display timer if requested
    if timer_on_flag
        tic
    end

    % Save current state
    x_traj(:, i) = x;

    % Check for "holdState" and how long to hold
    if holdState_flag
        hold_idx = find(strcmpi(varargin, "holdState")) + 1;
        if hold_idx <= length(varargin)
            hold_duration = varargin{hold_idx};
        else
            error('"holdState" must be followed by the number of steps to hold the state.');
        end
    else
        hold_duration = 0;
    end
    
    % Determine control action
    if holdState_flag && i < hold_duration
        u = zeros(size(B_disc, 2), 1); % zero input during hold
        x = x_init; % frozen state
    else
        if abs(n_outputs) == 1
            u = input_func(x_traj, u_traj, t);
        elseif nargout > 2
            [u, ~, g, cvx_optval_i] = input_func(noisy_x, u_traj, t);
            if iscell(g)
                alpha_i = g{2};
                g_i = g{1};
            else
                g_i = g;
                alpha_i = [];
            end
        end
        % Update state
        x = A_disc * x + B_disc * u;

        if timer_on_flag
            computation_times(i) = toc;
        end
    end
    
    % Store input
    u_traj(:, i) = u;

     if timer_on_flag
        varargout{1} = computation_times;
     end

    % Store optional outputs
    if nargout > 2 && ~timer_on_flag
        g_cell{1, i} = g_i;
        cvx_optval(1, i) = cvx_optval_i;
        if nargout > 3
            alpha_cell{1, i} = alpha_i;
        end
    end
end
end

