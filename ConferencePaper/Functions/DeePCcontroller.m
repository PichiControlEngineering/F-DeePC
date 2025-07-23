function [u, y, g, quadprog_optval] = DeePCcontroller(y, u, t, T_ini, N, H, r, lambda_g, lambda_y, u_bounds, y_bounds, uInactive)
% The DeePC controller employs an online optimization based on data from a
% "data-gathering" step. 
rng('shuffle')

%optional inputs
if ~exist("r", "var")
r = @(t) 0;
end

if ~exist("lambda_g","var")
lambda_g = 0;
end

if ~exist("lambda_y","var")
lambda_y = [];
end

if ~exist("y_bounds","var")
    y_bounds = [];
end

if ~exist("u_bounds","var")
    u_bounds = [];
end

if exist("uInactive","var")
    if strcmpi(uInactive,'Fixed')
        rng(t)
    else
        error("set variable to uInactive 'PRBSinputs', or leave empty")
    end
end
% Extract the data properties from H = createMatrices()
    f_pre_r = H.f;
    H_p = H.H_p;
    U_f = H.U_f;
    Y_f = H.Y_f;
    if isfield(H, 'Ts')
        Ts = H.Ts;
    else 
        Ts = 1;
    end

    H = H.H;
    
    l = size(U_f, 2);
    if t < T_ini
        % Before DeePC is triggered
        % u = 0;
        u = 3*(rand(1)-0.5);
        y = 0;
        g=zeros(l,1);
        quadprog_optval = NaN;
    else

        % Constructing r matrix
        d_t = 0:Ts:N-Ts;
        r_mat = r(t+d_t)';

        f = (r_mat'*f_pre_r)';
        
        % Implement a special procedure on t = T_ini, where we append 1
        % zero to the l_ini vector
        if t == T_ini
            y_ini = [y(:, t-T_ini+1:t)];
            u_ini = [0, u(:, t-T_ini+1:t-1)];
        else
            %regular procedure
            y_ini = y(:, t-T_ini+1:t);
            u_ini = u(:, t-T_ini:t-1);
        end
        l_ini = [u_ini'; y_ini'];
        
        W = 2*(H + eye(size(H)).*lambda_g);
        quadprogopts = optimoptions("quadprog", "MaxIterations", 1e3,"ConstraintTolerance", 1e-8, 'Display','off');

        % OPTIONAL: ADD SAFETY CONSTRAINTS
        A_ineq = [];
        b_ineq = [];
        
        if ~isempty(u_bounds)
            A_ineq = [-U_f; U_f];
            b_ineq = repmat([-u_bounds(1); u_bounds(2)], N,1);
        end
        
        if ~isempty(y_bounds)
            A_ineq = [A_ineq; -Y_f; Y_f]; 
            b_ineq = [b_ineq; repmat([-y_bounds(1); y_bounds(2)], N,1)];
        end

        %If the
        % slack variable is given, define slack variable problem
        if isempty(lambda_y)
            % Start regular optimalization
            [g,quadprog_optval,exitflag] = quadprog(W,f,A_ineq,b_ineq,H_p,l_ini,[],[],[],quadprogopts);
        else
            % QP with slack variable
            slack_var_mat = eye(T_ini).*lambda_y;
            W_slack = blkdiag(W,slack_var_mat);
            H_slack = [H_p, [zeros(size(H_p,1) - T_ini, T_ini); eye(T_ini)]];
            f_slack = [f; zeros(size(slack_var_mat,1),1)]; 
            A_ineq_slack = [A_ineq, zeros(size(A_ineq,1),size(slack_var_mat,1))];
            % Start slack variable optimalization

            [g_slack,quadprog_optval,exitflag] = quadprog(W_slack,f_slack,A_ineq_slack,b_ineq, H_slack,l_ini,[],[],[],quadprogopts);
            g = g_slack(1:end-T_ini); % Extract g
        end
        
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
    end
end 