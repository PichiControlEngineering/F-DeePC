function [varargout] = traj2Hankel(w,L,varargin)

% Inputs: w (trajectory), L (number of hankel rows)
% varargin: "split_hankels", this enables the optional splitting of the
% Hankel matrices, given U_p, Y_p, U_f, Y_f, with the former being
% length T_ini, and the latter being N
% Varargin{1}: must be "split" if you want to split hankel
% Varargin{2}: will be T_ini, and N will be L - T_ini

[n_trajectories, T] = size(w);

hankel_rows = T*n_trajectories; % Number of rows in the Hankel matrix
hankel_cols = T;     % Number of columns in the Hankel matrix

% Initialize Hankel matrix with zeros
H = zeros(hankel_rows, hankel_cols);

% Fill the Hankel matrix
for i = 1:hankel_cols
    for j = 1:(hankel_cols - i + 1)
        H(n_trajectories*(j-1)+1:n_trajectories*j, i) = w(:, j+i-1);
    end
end
trunc_rows = round(L*n_trajectories); % Number of rows in the Hankel matrix
trunc_cols = round((T - L) + 1);     % Number of columns in the Hankel matrix

if ~isempty(varargin) 
    if varargin{1} == "split"
        T_ini = varargin{2};
    
        H_U = H(1:2:trunc_rows, 1:trunc_cols);
        H_Y = H(2:2:trunc_rows, 1:trunc_cols);
    
        varargout{1} = H_U(1:T_ini, :); % U_p
        varargout{2} = H_Y(1:T_ini, :); %Y_p
        varargout{3} = H_U(T_ini+1:end, :); % U_f
        varargout{4} = H_Y(T_ini+1:end, :); % Y_f
    end
else
    varargout{1} = H(1:trunc_rows, 1:trunc_cols);
end
end