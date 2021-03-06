function [sim_results] = PolytopicSysSimSingle(sys, sim_params)
    
    

% function [sim_results] = PolytopicSysSimSingle(...
%                     sys_real, sys_hat, K, L, x_0, x_hat_0, u_0, N)
%     %POLYTOPICSYSSIMSINGLE simulates a single sim
%     arguments
%         sys_real
%         sys_hat
%         K
%         L
%         x_0
%         x_hat_0
%         u_0
%         N
%     end
    
    % System Values from sys struct
    %sys_real
    A = sys.sys_real.A;
    B = sys.sys_real.B;
    C = sys.sys_real.C;
    D = sys.sys_real.D;
    %sys_hat
    A_hat = sys.sys_hat.A;
    B_hat = sys.sys_hat.B;
    C_hat = sys.sys_hat.C;
    D_hat = sys.sys_hat.D;
    %K
    K = sys.K;
    %L
    L = sys.L;
    
    %x_0
    x_0 = sys.x_0;
    %x_hat_0
    x_hat_0 = sys.x_hat_0;
    %u_0
    u_0 = sys.u_0;
    
    % Simulation Parameters
    %Standard
    N = sim_params.N;
    u_k = @(x) sim_params.u_k(x, K); %updated to be just of x
    SigmaInv = sim_params.SigmaInv;
    
    
    
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
    sim_results.Z = zeros(1,N);
    
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
        u = u_k(x_hat);%K * x_hat;
        
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
        sim_results.Z(:,k) = (y-y_hat)'*SigmaInv*(y-y_hat);
    end
    
    % Additional Calc
    sim_results.E = sim_results.X - sim_results.X_hat;
    sim_results.R = sim_results.Y - sim_results.Y_hat;
%     sim_results.Z = ...
%         sim_results.R.' * SigmaInv * sim_results.R;
end

