function [x_traj, u_traj] = systemDataGen(A,B,C,D,Ts,T,x_init,u_func)
    %% Defining & discretizing the system
    if ~isempty(Ts) 
        % Optional discretizetion in the function, will skip if Ts is not provided
        ss_cont = ss(A,B,C,D);
        ss_disc = c2d(ss_cont,Ts,'tustin');
        
        A = ss_disc.A; B = ss_disc.B;
    else
        % If the system is already given discretized, Ts is set to 1
        Ts = 1;
    end
    
    %If x_init is accidentally given as row vector, swap its columns
    if size(x_init) == [1, size(A,1)]
        warning('Initial condition is given as a row vector, transposing x_init to column vector')
        x_init = x_init';
    end
    x_traj = zeros(size(A,1),T/Ts);
    u_traj = zeros(size(B,2),T);
    
    x = x_init;
    %Start Step_states basic
    for i = 1:T/Ts
        t = (i-1)*Ts; %Time goes from 0 to T
        % Save current state
        x_traj(:, i) = x; 

        % OPTIONAL: Add noise; given by sigma_e
        u = u_func(x_traj, u_traj, i);
        
        
        % Update state
        x = A * x + B .* u;
        u_traj(:, i) = u;
    end
end
