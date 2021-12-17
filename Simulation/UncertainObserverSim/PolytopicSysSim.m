% function [X_data, X_hat_data, varargout] = PolytopicSysSim(...
%                 A_all, L, ALPHA_real, ALPHA_hat, X_0, X_hat_0,...
%                 N, B_all, C_all, D_all, K, u_0)
%     
%     arguments
%         A_all (:,:,:) double % Plan Sys Matrix
%         L (:,:) double % Observer Gain Matrix
%         ALPHA_real (:,:) double % Plant Sys Params
%         ALPHA_hat (:,:) double% Observer Sys Params
%         X_0 (:,1) double% Plant Sys ICs
%     % Optional Parameters
%         X_hat_0 (:,1) double = X_0 % Observer ICs (assume X_0 => E_0 = 0)
%         N (1,1) uint16 = 100 % Num simu steps
%         B_all (:,:,:) double = zeros(size(A_all,1),1) % B matrix (assume zeros(n,p=1))
%         C_all (:,:,:) double = eye(size(A_all,1)) % C matrix (assume eye(n))
%         D_all (:,:,:) double = zeros(size(C_all,1),size(A_all,1)) % D matrix (assume zeros(q,n))
%         K (:,:) double = zeros(size(B_all,2),size(A_all,1)) % K matrix (assume zeros(p,n))
%         u_0 (:,1) = zeros(size(B_all,2)) % Initial u_0 value... assume 0
%         u_x function_handle = @(x,K) K * x % u(x) = K * x
%     end
%     
%     % Polytopic LPV System Dimensions
%     n = size(A_all,1);
%     m = size(A_all,3);
%     p = size(B,2);
%     q = size(C,1);
%     
%     % Simulation Parameter Dimensions
%     num_x_0 = size(X_0,2);
%     num_real = size(ALPHA_real,2);
%     num_hat = size(ALPHA_hat,2);
%     
%     % Save Data Arrays
%     X_data = zeros(N, n, num_x_0, num_real, num_hat);
%     X_hat_data = zeros(N, n, num_x_0, num_real, num_hat);
%     U_data = zeros(N, p, num_x_0, num_real, num_hat);
%     Y_data = zeros(N, q, num_x_0, num_real, num_hat);
% %     E_data = zeros(N, n, num_x_0, num_alpha_real, num_alpha_hat);
% %     R_data = zeros(N, q, num_x_0, num_alpha_real, num_alpha_hat);
%     
% %     A_real_data = zeros(n,n,num_real,num_hat);
% %     A_hat_data = zeros(n,n,num_real,num_hat);
%     
%     % Cell of sys w/ A, B, C, D
%     sys_all = {A_all, B_all, C_all, D_all};
% 
%     % ndgrid of system and parameters
%     [x_0_data, alpha_real_data, alpha_hat_data] = ndgrid(...
%         num2cell(X_0,1), num2cell(ALPHA_real,1), num2cell(ALPHA_hat,1));
% 
%     sysCalc = @(sys_all, alpha_data) cellfun(@(alpha) cellfun(@(A_all)...
%                 PolytopicSysCalc(A_all, alpha),...
%                     sys_all, 'UniformOutput',false),...
%                     alpha_data, 'UniformOutput',false);
% 
% 
%     sys_real_data = sysCalc(sys_all, alpha_real_data);
%     sys_hat_data = sysCalc(sys_all, alpha_real_data);
%     
%     sim_results_data = 
%     

    % Sim Results - X, X_hat, U, Y, E, R
%     sim_result_data = cell(6, num_x_0, num_real, num_hat)
    
    
    
%     % Initial Conditions loop
%     for idx_x_0 = 1:num_x_0
%         x_0 = X_0(:,idx_x_0);
%         x_hat_0 = X_hat_0(:,idx_x_0);
%         % Alpha_real loop
%         for idx_real = 1:num_real
%             Alpha_real = ALPHA_real(:,idx_real);
%             % Alpha_hat loop
%             for idx_hat = 1:num_hat
%                 Alpha_hat = ALPHA_hat(:,idx_hat);
%                 
%                 % Polytopic Calculations
%                 A_real = zeros(n);
%                 A_hat = zeros(n);
%                 for i = 1:m
%                     A_real = A_real + Alpha_real(i) * A_all(:,:,i);
%                     A_hat = A_hat + Alpha_hat(i) * A_all(:,:,i);
%                 end
%                 A_real_data(:,:,idx_real,idx_hat) = A_real;
%                 A_hat_data(:,:,idx_real,idx_hat) = A_hat;
%                 Delta_A = A_real - A_hat;
% %                 delta = norm(Delta_A);
%                 
% %                 % Sim Setup
% %                 X = zeros(N,n);
% %                 X_hat = zeros(N,n);
% %                 U = zeros(N,p);
% %                 Y = zeros(N,q);
% %                 R = zeros(N,q);
%     
% %                 % Error Values
% %                 E = zeros(N,n);
%                 % Inital setup (k=0)
%                 x = x_0;
%                 x_hat = x_hat_0;
%                 u = u_0;
%                 e = x - x_hat;
%                 %$y_{0} = C x_{0} + D u_{0}$
%                 y = C * x + D * u;
%                 %$y_hat_{0} = C x_hat_{0} + D u_{0}
%                 y_hat = C * x_hat + D * u;
%                 
%                 % Simulation
%                 for k = 1:N
%                     % Control
% %                     u = 0; % No Input
%                     u = u_x(x_k,K);
%                     
%                     
%                     % Update Eqs
%                     %$x_{k} = A x_{k-1} + B u_{k-1}$
%                     x = A_real * x + B * u;
%                     %$x_hat_k = A_hat x_hat_{k-1}
%                     %           + B_hat u_{k-1} 
%                     %           + L_hat (y_{k-1} - C x_hat_{k-1})$
%                     x_hat = A_hat * x_hat...
%                             + B * u ...
%                             + L * (y - y_hat);
%                     
%                     %$e_k = x_k - x_hat_k$
%                     e = x - x_hat;
%                     %$r_k = C * e_k$
%                     r = C * e;
%                     
%                     %$y_{k-1} = C x_{k-1} + D u_{k-1}$
%                     y = C * x + D * u;
%                     %$y_hat_{k-1} = C x_hat_{k-1} + D u_{k-1}
%                     y_hat = C * x_hat + D * u;
%     
%                     % Save Data
% %                     X(k,:) = x;
% %                     X_hat(k,:) = x_hat;
% %                     U(k,:) = u;
% %                     Y(k,:) = y;
% %                     E(k,:) = e;
% %                     R(k,:) = r;
%                     X_data(k,:,idx_x_0,idx_real,idx_hat) = x;
%                     X_hat_data(k,:,idx_x_0,idx_real,idx_hat) = x_hat;
% %                     Y_data(k,:,idx_x_0,idx_real,idx_hat) = y;
% %                     E_data(k,:,idx_x_0,idx_real,idx_hat) = e;
% %                     R_data(k,:,idx_x_0,idx_real,idx_hat) = r;
%                 end
%             end
%         end
%     end
%     E_data = X_data - X_hat_data;
%     Y_data = 
    
function [results,sim_params] = PolytopicSysSim(...
        A_all, C_all, L, alpha_real_all, alpha_hat_all, X_0, ...
        X_hat_0, N, B_all, D_all, K, U_0, u_k, SigmaInv...
        )
    arguments
        A_all (:,:,:) double % Plant Sys Matrix
        C_all (:,:,:) double % Plant Output Matrix
        L (:,:) double % Observer Gain Matrix
        alpha_real_all (:,:) double % Plant Sys Params
        alpha_hat_all (:,:) double% Observer Sys Params
        X_0 (:,:,:) double% Plant Sys ICs
    % Optional Parameters
        X_hat_0 = [] % Observer ICs (assume [] => X_hat_0 = X_0 => E_0 = 0)
        N (1,1) double = 100 % Num simu steps
        B_all (:,:,:) double = zeros(size(A_all,1),1) % B matrix (assume zeros(n,p=1))
        D_all (:,:,:) double = zeros(size(C_all,1),size(B_all,2)) % D matrix (assume zeros(q,n))
        K (:,:) double = zeros(size(B_all,2),size(A_all,1)) % K matrix (assume zeros(p,n))
        U_0 (:,:,:) double = zeros(size(B_all,2)) % Initial u_0 value... assume 0
        u_k function_handle = @(x,K) K * x % u(x) = K * x
        SigmaInv (:,:) = eye(size(C_all,1));
    end     
    %% Input Processing
    % System Matrices
    sys_poly = struct( ...
        'A', A_all, ...
        'B', B_all, ...
        'C', C_all, ...
        'D', D_all ...
        );
    K_all = {K};
    L_all = {L};
    % Initial Conditions
    X_0_all = num2cell(X_0,1);
    if numel(X_hat_0) == 0, X_hat_0_all = X_0_all(1);
    else, X_hat_0_all = num2cell(X_hat_0,1); end
    
    U_0_all = num2cell(U_0,1);
    % Polytopic Sys Calculations
    sys_real_all = PolytopicSysCalc(sys_poly, alpha_real_all);
    real_all = struct(  'alpha', num2cell(alpha_real_all,1),...
                        'sys', sys_real_all);
    sys_hat_all = PolytopicSysCalc(sys_poly, alpha_hat_all);
    hat_all = struct(  'alpha', num2cell(alpha_hat_all,1),...
                        'sys', sys_hat_all);
                    
    % Grid of Simulation Systems
    [real_grid, hat_grid, K_grid, L_grid,...
            X_0_grid, X_hat_0_grid, U_0_grid] = ...
        ndgrid(...
            real_all, hat_all, K_all, L_all,...
            X_0_all, X_hat_0_all, U_0_all...
        );
    
    if numel(X_hat_0) == 0, X_hat_0_grid = X_0_grid; end
    
    %% Simulation Systems Struct Def
    sim_sys = struct(...
        'alpha_real', reshape({real_grid.alpha}, size(real_grid)),...
        'alpha_hat', reshape({hat_grid.alpha}, size(hat_grid)),...
        'sys_real', reshape({real_grid.sys}, size(real_grid)),...
        'sys_hat', reshape({hat_grid.sys}, size(hat_grid)),...
        'K', K_grid,...
        'L', L_grid,...
        'x_0', X_0_grid,...
        'x_hat_0', X_hat_0_grid,...
        'u_0', U_0_grid...
        );
    
    %% Simulation Paramteres
    sim_params.N = N;
    sim_params.u_k = u_k;
    sim_params.SigmaInv = SigmaInv;
    
    %% Simulate System
    sim_results = arrayfun(...
        @(sys) PolytopicSysSimSingle(sys, sim_params),...
        sim_sys);
    
    %% Concatinate Results
    results =  cell2struct(...
        [struct2cell(sim_sys); 
         struct2cell(sim_results)],...
        [fieldnames(sim_sys);
         fieldnames(sim_results)]...
         );
end