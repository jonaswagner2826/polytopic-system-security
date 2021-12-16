function PolytopicSysSim(A_all, L, ALPHA_real, ALPHA_hat, X_0, X_hat_0, N, B, C, D, K)
    
    arguments
        A_all % Plan Sys Matrix
        L_all % Observer Gain Matrix
        ALPHA_real % Plant Sys Params
        ALPHA_hat % Observer Sys Params
        X_0 % Plant Sys ICs
    % Optional Parameters
        X_hat_0 = X_0 % Observer Sys ICs (assume = X_0 => E_0 = 0)
        N = 100 % Num simu steps
        B_all = zeros(size(A_all,1),1) % B matrix (assume zeros(n,p=1))
        C_all = eye(size(A_all,1)) % C matrix (assume eye(n))
        D_all = zeros(size(C,1),size(A_all,1)) % D matrix (assume zeros(q,n))
        K_all = zeros(size(B,2),size(A_all,1)) % K matrix (assume zeros(p,n))
    end
    
    % System Dimensions
    n = size(A_all,1);
    m = size(A_all,3);
    p = size(B,2);
    q = size(C,1);
    
    num_x_0 = size(X_0,2);
    num_real = size(ALPHA_real,2);
    num_hat = size(ALPHA_hat,2);
    
    % Save Data Arrays
    X_data = zeros(N, n, num_x_0, num_real, num_hat);
    X_hat_data = zeros(N, n, num_x_0, num_real, num_hat);
    Y_data = zeros(N, q, num_x_0, num_real, num_hat);
%     E_data = zeros(N, n, num_x_0, num_alpha_real, num_alpha_hat);
%     R_data = zeros(N, q, num_x_0, num_alpha_real, num_alpha_hat);
    
    A_real_data = zeros(n,n,num_real,num_hat);
    A_hat_data = zeros(n,n,num_real,num_hat);
    
    
    
    
    % Initial Conditions loop
    for idx_x_0 = 1:num_x_0
%         x_0 = X_0(:,idx_x_0);
        x_hat_0 = X_hat_0(:,idx_x_0);
        % Alpha_real loop
        for idx_real = 1:num_real
            Alpha_real = ALPHA_real(:,idx_real);
            % Alpha_hat loop
            for idx_hat = 1:num_hat
                Alpha_hat = ALPHA_hat(:,idx_hat);
    
                % Polytopic Calculations
                A_real = zeros(n);
                A_hat = zeros(n);
                for i = 1:m
                    A_real = A_real + Alpha_real(i) * A_all(:,:,i);
                    A_hat = A_hat + Alpha_hat(i) * A_all(:,:,i);
                end
                A_real_data(:,:,idx_real,idx_hat) = A_real;
                A_hat_data(:,:,idx_real,idx_hat) = A_hat;
                Delta_A = A_real - A_hat;
                delta = norm(Delta_A);
    
%                 % Sim Setup
%                 X = zeros(N,n);
%                 X_hat = zeros(N,n);
%                 U = zeros(N,p);
%                 Y = zeros(N,q);
%                 R = zeros(N,q);
    
%                 % Error Values
%                 E = zeros(N,n);
    
                % Inital setup (k=0)
                x = x_0;
                x_hat = x_hat_0;
                e = x - x_hat;
                
                % Simulation
                for k = 1:N
                    % Control
                    u = 0; % No Input
                    
                    
                    % Update Eqs
                    %$y_{k-1} = C x_{k-1} + D u_{k-1}$
                    y = C * x + D * u;
                    %$x_k = A x_{k-1} + B u_{k-1}$
                    x = A_real * x + B * u;
                    %$x_hat_k = A_hat x_hat_{k-1}
                    %           + B u_{k-1} 
                    %           + L (y_{k-1} - C x_hat_{k-1})$
                    x_hat = A_hat * x_hat...
                            + B * u ...
                            + L * (y - C * x_hat);
                    %$e_k = x_k - x_hat_k$
                    e = x - x_hat;
                    %$r_k = C * e_k$
                    r = C * e;
    
                    % Save Data
%                     X(k,:) = x;
%                     X_hat(k,:) = x_hat;
%                     U(k,:) = u;
%                     Y(k,:) = y;
%                     E(k,:) = e;
%                     R(k,:) = r;
                    X_data(k,:,idx_x_0,idx_real,idx_hat) = x;
                    X_hat_data(k,:,idx_x_0,idx_real,idx_hat) = x_hat;
                    Y_data(k,:,idx_x_0,idx_real,idx_hat) = y;
%                     E_data(k,:,idx_x_0,idx_real,idx_hat) = e;
%                     R_data(k,:,idx_x_0,idx_real,idx_hat) = r;
                end
            end
        end
    end
    E_data = X_data - X_hat_data;
end