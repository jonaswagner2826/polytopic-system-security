function [sim_results] = PolytopicSysSimSingle(...
                    sys_real, sys_hat, K, L, x_0, x_hat_0, u_0, N)
    %POLYTOPICSYSSIMSINGLE simulates a single sim
    arguments
        sys_real
        sys_hat
        K
        L
        x_0
        x_hat_0
        u_0
        N
    end
    
    A = sys_real{1}{1};
    B = sys_real{1}{2};
    C = sys_real{1}{3};
    D = sys_real{1}{4};
    
    A_hat = sys_hat{1}{1};
    B_hat = sys_hat{1}{2};
    C_hat = sys_hat{1}{3};
    D_hat = sys_hat{1}{4};
    
    % System Dimensions
    n = size(A,1);
    p = size(B,2);
    q = size(C,1);
    
    % Data Arrays
    sim_results.X = zeros(n,N);
    sim_results.X_hat = zeros(n,N);
    sim_results.U = zeros(p,N);
    sim_results.Y = zeros(q,N);
    sim_results.Y_hat = zeros(q,N);
    
    % k = 0
    x = x_0;
    x_hat = x_hat_0;
    u = u_0;
    y = C * x + D * u;
    y_hat = C_hat * x_hat + D * u;
    
    
    for k = 1:N
        %$x_{k} = A x_{k-1} + B u_{k-1}$
        x = A * x + B * u;
        %$x_hat_k = A_hat x_hat_{k-1}
                   %+ B_hat u_{k-1} 
                   %+ L_hat (y_{k-1} - C x_hat_{k-1})$
        x_hat = A_hat * x_hat...
                + B_hat * u ...
                + L * (y - y_hat);
        %$u_{k} = K * x_hat_{k}$
        u = K * x_hat;
        
        %$y_{k} = C x_{k} + D u_{k}$
        y = C * x + D * u;
        %$y_hat_{k-1} = C x_hat_{k-1} + D u_{k-1}$
        y_hat = C_hat * x_hat + D_hat * u;
        
        % Store Data
        sim_results.X(:,k) = x;
        sim_results.X_hat(:,k) = x_hat;
        sim_results.U(:,k) = u;
        sim_results.Y(:,k) = y;
        sim_results.Y_hat(:,k) = y_hat;        
    end
    
    % Additional Calc
    sim_results.E = sim_results.X - sim_results.X_hat;
    sim_results.R = sim_results.Y - sim_results.Y_hat;
end

