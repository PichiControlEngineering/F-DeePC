function [Up, Yp, Uf, Yf] = createHankel(u_d,y_d,Tini,N)
    nu = size(u_d,1);
    ny = size(y_d,1);
    
    % Ensure U and Y are of equal size
    if size(u_d,2) ~= size(y_d,2)
        error('U and Y do not have same length')
    end
    
    % Ensure Willems fundemental lemma holds
    if size(u_d,2) < nu*(Tini+N+1) 
        error('T >= n_u*(T_ini + N + n); lower T_ini or N, or create a longer data sample')
    end

    T = size(u_d,2) - N - Tini + 1;
    Lmax = floor(T/3);

    Up = zeros(Tini*nu, T); Uf = zeros(N*nu, T);
    Yp = zeros(Tini*ny, T); Yf = zeros(N*ny, T);

    for i = 1:Tini
        Up((i-1)*nu+1:i*nu, :) = u_d(:, i :i+T-1);
        Yp((i-1)*ny+1:i*ny, :) = y_d(:, i:i+T-1);
    end

    for i = 1:N
        Uf((i-1)*nu+1:i*nu, :) = u_d(:,i+Tini:i+Tini+T-1);
        Yf((i-1)*ny+1:i*ny, :) = y_d(:, i+Tini:i+Tini+T-1);
    end
end