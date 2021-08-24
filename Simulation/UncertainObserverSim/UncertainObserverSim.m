% Uncertain Observer Toy System Simulation
% Jonas Wagner
% 2021-08-11 @ 6:12PM

clear
close all

% Sim Parameters
N = 20;

% Optional Parts
runSim = false;
testCtrb = false;
testObsv = false;
testL = false;
plotAllResiduals = false;
plotResidualStatistics = false;
independentSim = true;


% % Toy System Def
% A(:,:,1) = [-0.80, 0.25; 0.25,-0.30];
% A(:,:,2) = [ 0.30, 0.70; 0.70, 0.00];
% A(:,:,3) = [-0.30, 0.65; 0.55, 0.10];
% A(:,:,4) = [ 0.55,-0.20;-0.40,-0.30];
% 
% B = [1.5; -0.5];
% C = [1, 0];
% D = 0;

% % Toy Sys Def 2
% V = [0.2785, 0.9575;  0.5469, 0.9649]; % generated randomly
% A(:,:,1) = V * diag([0.8,0.9]) * inv(V);
% A(:,:,2) = V * diag([0.85,-0.95]) * inv(V);
% A(:,:,3) = V * diag([-0.85,0.95]) * inv(V);
% A(:,:,4) = V * diag([-0.9,-0.8]) * inv(V);

% B = 0;
% C = [0.2, 0.5];
% D = 0;

% Toy Sys Def 3
V =[0.7780    0.3669    0.8727    0.6220;
    0.4267    0.7948    0.2858    0.0751;
    0.2800    0.0387    0.6568    0.9668;
    0.3335    0.7267    0.2319    0.6100]; % Generated Randomly
A(:,:,1) = V * diag([-0.410    0.9096    0.7019   -0.7042]) * inv(V);
A(:,:,2) = V * diag([0.7324   -0.7779   -0.5655    0.1609]) * inv(V);
A(:,:,3) = V * diag([0.1511    0.7571    0.5847    0.0015]) * inv(V);
A(:,:,4) = V * diag([0.8152   -0.4052    0.3449   -0.5620]) * inv(V);

B =[0.2779   -0.1673;
    0.7863    0.4796;
   -0.8786    0.7859;
   -0.6485   -0.9483];

C = [0.4745   -0.1475   -0.3936    0.0739;
    -0.6376   -0.8037    0.5602    0.5381];

D = 0;

% System Dimensions
n = size(A,1);
m = size(A,3);
p = size(B,2);
q = size(C,1);

% Test Controllability
if testCtrb
    for i = 1:m
        rank(ctrb(A(:,:,i),B))
    end
end

% Test Observability
if testObsv
    for i = 1:m
        rank(obsv(A(:,:,i),C))
    end
end


% Initial Conditions
%single set of I.C.s
if independentSim
randTransform = [0.4664    0.8965;    0.4564    0.4274;    0.1278    0.1758;    0.6399    0.2949];
x_0_mag = 1;
x_0_theta = pi/4;
x_0 = randTransform * x_0_mag * [cos(x_0_theta); sin(x_0_theta)];
x_hat_0 = x_0;
end

%multiple I.C.s
X_0_mag = [1];
X_0_theta = [0, pi/4, pi/2, 3*pi/4, pi];
X_0 = round(X_0_mag * [cos(X_0_theta); sin(X_0_theta)],2);
X_hat_0 = X_0;


% Parameter Definition
%single set of alphas...
if false
Alpha_real = [1,0,0,0];%[0.25, 0.25, 0.25, 0.25];
% Alpha_hat = [0,0,0,1];%[0.375, 0.125, 0.125, 0.125];
Alpha_hat = [0.999,0.001,0,0]
end

%Multiple Alphas
num_extra_alpha_hat = 50;
ALPHA_real = [eye(m), [0.25;0.25;0.25;0.25]];% normalize(rand(m),1,'norm',1)];
ALPHA_hat = [eye(m), [0.25;0.25;0.25;0.25],...
    normalize(rand(m, num_extra_alpha_hat),1,'norm',1)];


% CVX Designing
% % Control Design: u = K * x_hat
% tol = 1e-6;
% cvx_clear
% cvx_begin sdp
%     variable Q(p,n)
%     variable P(n,n) symmetric
%     maximize(trace(P))
%     subject to
%         P >= 0;
%         for i = 1:m
%             [P, (A(:,:,i)* P + B * Q)';
%             (A(:,:,i) * P + B * Q) , P] >= 0;
%         end
% cvx_end
% 
% K = Q * inv(P)


% Observer Design: x_hat = A_hat * x_hat + B * u + L * (y - y_hat)
cvx_begin sdp quiet
    variable X(n,q)
    variable Q(n,n) symmetric
    tol = 1e-6;
    minimize(trace(Q))
%     maximize(trace(Q));
    subject to
        for i = 1:m
            Q >= tol*eye(n);
            [Q, (Q*A(:,:,i) - X*C)';
             Q*A(:,:,i) - X*C, Q] >= tol*eye(2*n);
        end
cvx_end

L = inv(Q) * X

% Testing L
if testL
    disp('Testing L Design')
    for i = 1:m
        A(:,:,i)
        A_LC = A(:,:,i) - L * C
        eig_A_LC = eig(A_LC)
    end
    A_real
    A_LC = A_real - L*C
    eig_A_LC = eig(A_LC)
    A_hat
    A_LC = A_hat - L*C
    eig_A_LC = eig(A_LC)
end


% Save Data
% X_data = zeros(size(ALPHA_real,2), size(ALPHA_hat,2), N, n);
% X_hat_data = zeros(size(ALPHA_real,2), size(ALPHA_hat,2), N, n);
R_data = zeros(N, q, size(X_0,2), size(ALPHA_real,2), size(ALPHA_hat,2));

%% Simulation
if runSim
% Initial Conditions loop
for idx_x_0 = 1:size(X_0,2)
    x_0 = X_0(:,idx_x_0);
    x_hat_0 = X_hat_0(:,idx_x_0);
    % Alpha_real loop
    for idx_real = 1:size(ALPHA_real,2)
        Alpha_real = ALPHA_real(:,idx_real);
        % Alpha_hat loop
        for idx_hat = 1:size(ALPHA_hat,2)
            Alpha_hat = ALPHA_hat(:,idx_hat);

            % Polytopic Calculations
            A_real = zeros(n);
            A_hat = zeros(n);
            for i = 1:m
                A_real = A_real + Alpha_real(i) * A(:,:,i);
                A_hat = A_hat + Alpha_hat(i) * A(:,:,i);
            end
            Delta_A = A_real - A_hat;
            delta = norm(Delta_A);

            % Sim Setup
            X = zeros(N,n);
            X_hat = zeros(N,n);
            U = zeros(N,p);
            Y = zeros(N,q);
            R = zeros(N,q);

            % Error Values
            E = zeros(N,n);
            E_calc = zeros(N,n);
            E_calc_2 = zeros(N,n);
            E_norm = zeros(N,1);
            E_norm_calc = zeros(N,1);
            E_norm_calc_2 = zeros(N,1);
            E_norm_bound = zeros(N,1);

            % Inital setup (k=0)
            x = x_0;
            x_hat = x_hat_0;
            e = x - x_hat;
            e_calc_2 = x_0 - x_hat_0;
            e_norm_calc_2 = norm(e_calc_2);

            % Simulation
            for k = 1:N
                % Control
                u = 0; % No Input

                % Error Calc 2(recursive)
                %e_k = Delta_A * x_{k-1} + (A_hat - L*C) * e_{k-1};
                e_calc_2 = Delta_A * x_hat + (Delta_A + A_hat - L*C) * e_calc_2;
                E_calc_2(k,:) = e_calc_2;

                % Error norm cal 2 (recursive)
                %\norm{e_k} \leq \delta * \norm{\hat{x}_{k-1}} + (\delta +
                %\norm{\hat{A} - LC}) \norm{e_{k-1}}
                e_norm_calc_2 = delta * norm(x_hat) + (delta + norm(A_hat - L*C)) * e_norm_calc_2;
                E_norm_calc_2(k,:) = e_norm_calc_2;

                % Update Eqs
                %y_{k-1} = C x_{k-1} + D * u_{k-1}
                y = C * x + D * u;
                %x_k = A * x_{k-1} + B * u_{k-1}
                x = A_real * x + B * u;
                %x_hat_k = A_hat * x_hat_{k-1} + B * u_{k-1} + L * (y_{k-1} - C * x_hat_{k-1})
                x_hat = A_hat * x_hat + B * u + L * (y - C * x_hat);
                %e_k = x_k - x_hat_k
                e = x - x_hat;
                %r_k = C * e_k
                r = C * e;

                % Save Data
                X(k,:) = x;
                X_hat(k,:) = x_hat;
                U(k,:) = u;
                Y(k,:) = y;
                E(k,:) = e;
                R(k,:) = r;

                % Error Calc Estimate
                e_calc = x_0 - x_hat_0;
                for i = 0:(k-1)
                    e_calc = e_calc + (A_hat - L * C)^(i) * Delta_A * A_real^(k - 1 - i) * x_0;
                end
                E_calc(k,:) = e_calc;

                % Error Norm
                E_norm(k,1) = norm(x - x_hat);

                % Error Norm Calc
                e_norm_calc = 0;
                for i = 0:(k-1)
                    e_norm_calc = e_norm_calc +...
                        norm((A_hat - L * C)^(i) * Delta_A * A_real^(k -1 - i)) * norm(x_0);
                end
                E_norm_calc(k,:) = e_norm_calc;

                % Error Norm Bound
                E_norm_bound(k,:) = norm(Delta_A) * norm(x_0) * ...
                    (norm(A_hat - L*C)^k - norm(A_real)^k)/...
                    (norm(A_hat - L*C) - norm(A_real));
            end

            % Detailed Plot (of single alpha)
            if false
            figure('Position', [0,0,1500,800])
            sgtitle(num2str(x_0))

            subplot(2,3,1)
            plot(X(:,1));
            hold on
            plot(X_hat(:,1));
            title('State and Estimate (X1)')
            legend('Actual', 'Estimate')

            subplot(2,3,4)
            plot(X(:,2));
            hold on
            plot(X_hat(:,2));
            title('State and Estimate (X2)')
            legend('Actual', 'Estimate')

            subplot(2,3,2)
            plot(E(:,1))
            hold on
            plot(E_calc(:,1))
            plot(E_calc_2(:,1))
            title('Error (X1)')
            legend('Actual','Calculated','calc2')

            subplot(2,3,5)
            plot(E(:,2))
            hold on
            plot(E_calc(:,2))
            plot(E_calc_2(:,2))
            title('Error (X2)')
            legend('Actual','Calculated','calc2')

            subplot(2,3,3)
            plot(E_norm)
            hold on
            plot(E_norm_calc)
            plot(E_norm_bound)
            title('Error Norm')
            legend('Actual','Sum of Norms','Closed-Form Bound')
            % ylim([0, 2 * max(E_norm)])

            subplot(2,3,6)
            plot((E_norm - E_norm_calc)./E_norm)
            hold on
            % plot((E_norm - E_norm_bound)./E_norm)
            % ylim([5 * min((E_norm - E_norm_calc)./E_norm),0])
            title('Error Norm Calc & Bound % Error')
            legend('Sum of Norms Error','Closed-Form Bound Error')
            end

            % Additional figures
            if false
            figure
            subplot(2,1,1)
            plot(((E-E_calc)./E)*100)
            title('% Error of Actual vs Calculated Error')
            legend('X_1', 'X_2')

            subplot(2,1,2)
            plot(((E-E_calc_2)./E)*100)
            title('% Error of Actual vs Calculated Error 2')
            legend('X_1', 'X_2')
            end

            % Save Data
            R_data(:,:,idx_x_0,idx_real,idx_hat) = R;
        end
    end
end
end
%% All Residual Plot
if plotAllResiduals
    figure('Position',[0,0,2.5e3,1.25e3])
    sgtitle(['Residual with varied $\alpha$, $\hat{\alpha}$, and ',...
        '$\hat{x}_0 = x_0$'], 'Interpreter','latex')
    for idx_x_0 = 1:size(X_0,2)
        for idx_real = 1:size(ALPHA_real,2)
            subplot(size(ALPHA_real,2),size(X_0,2),size(ALPHA_real,2)*(idx_x_0-1) + idx_real)
            hold on
            for idx_hat = 1:size(ALPHA_hat,2)
                plot(R_data(:,:,idx_x_0,idx_real,idx_hat),'DisplayName',...
                    ['$\hat{\alpha} = [', regexprep(num2str(...
                    ALPHA_hat(:,idx_hat)',2),'\s+',','), ']^T$'])
            end
            title(['$\alpha = [', regexprep(num2str(ALPHA_real(:,idx_real)',2),...
                '\s+',','), ']^T$ and $\hat{x}_0 = x_0 = [',...
                regexprep(num2str(X_0(:,idx_x_0)',1),'\s+',','),']^T$'], 'Interpreter','latex')
            legend('Interpreter','latex')
        end
    end
end


%% Residual Response Statistic Plots
residual_norm_p = inf;
axes = [];
if plotResidualStatistics
    delta_A_data = zeros(size(X_0,2),size(ALPHA_real,2),size(ALPHA_hat,2));
    R_statistic_data = zeros(size(X_0,2),size(ALPHA_real,2),size(ALPHA_hat,2));
    figure('Position',[0,0,9e2,1.3e3])%2.5e3,1.25e3])
    sgtitle(['Residual response statistic vs norm of $\Delta A$'],...
        'Interpreter','latex')
    for idx_real = 1:size(ALPHA_real,2)
        Alpha_real = ALPHA_real(:,idx_real)';
        ax = subplot(size(ALPHA_real,2),1,idx_real);
        axes = [axes, ax];
        hold on
        for idx_x_0 = 1:size(X_0,2)
            x_0 = X_0(:,idx_x_0)';
            for idx_hat = 1:size(ALPHA_hat,2)
                Alpha_hat = ALPHA_hat(:,idx_hat)';
                Delta_alpha = Alpha_real - Alpha_hat;
                Delta_A = zeros(n);
                for i = 1:m
                    Delta_A = Delta_A + Delta_alpha(i) * A(:,:,i);
                end
                delta_A = norm(Delta_A, 'fro'); % F-norm
                R = R_data(:,:,idx_x_0,idx_real,idx_hat);
                R_statistic = norm(R,residual_norm_p); % p-norm
%                 R_statistic = norm(R,2); % 2-norm = RMS
%                 R_statistic = norm(R,1); % 1-norm = Absolute Sum
%                 R_statistic = norm(R,inf); % inf-norm = Max Residual
                delta_A_data(idx_x_0,idx_real,idx_hat) = delta_A;
                R_statistic_data(idx_x_0,idx_real,idx_hat) = R_statistic;
            end
            scatter(reshape(delta_A_data(idx_x_0,idx_real,:),1,[]),...
                reshape(R_statistic_data(idx_x_0,idx_real,:),1,[]),...
                15,'filled','LineWidth',0.5,...
                'DisplayName', ['$x_0 = [', regexprep(num2str(...
                    X_0(:,idx_x_0)',2),'\s+',','), ']^T$']);
        end
        title(['$\alpha = [', regexprep(num2str(ALPHA_real(:,idx_real)',2),...
                '\s+',','), ']^T$'],'Interpreter','latex');
        xlabel('$||\Delta A||_F$','Interpreter','latex');
        ylabel(['$||r_k||_',num2str(residual_norm_p),'$'],'Interpreter','latex');
        legend('Interpreter','latex','Location','northwest')
    end
    linkaxes(axes,'xy')
end


%% Independent Sim
if independentSim  
close all

% Setup
Alpha_real = rand(m,1);%ALPHA_real(:,1);
Alpha_real = Alpha_real / sum(Alpha_real);
Alpha_hat = ALPHA_hat(:,5);

A_real = zeros(n);
A_hat = zeros(n);
for i = 1:m
    A_real = A_real + Alpha_real(i) * A(:,:,i);
    A_hat = A_hat + Alpha_hat(i) * A(:,:,i);
end

% Independent Sim
X_all = zeros(N,n,m);
for i = 1:m
    X_temp = zeros(N,n);
    x = x_0;
    for k = 1:N
        u = zeros(p,1);
        x = A(:,:,i) * x + B * u;
        X_all(k,:,i) = x;
    end
end

% Regular Sim
X_real = zeros(N,n);
X_hat = zeros(N,n);
x = x_0;
x_hat = x_hat_0;
for k = 1:N
    y = C * x + D * u;
    x_hat = A_hat * x_hat + B * u + L * (y - C*x_hat);
    x = A_real * x + B * u;
    X_real(k,:) = x;
    X_hat(k,:) = x_hat;
end

% Alpha Calc
ALPHA_tilde = zeros(N,m);
ALPHA_tilde_hat = zeros(N,m);
P = zeros(m); % Transform Matrix (assuming n = m)
for k = 1:N
    for i = 1:m
        P(:,i) = X_all(k,:,i);
    end
    ALPHA_tilde(k,:) = P \ X_real(k,:)';
    ALPHA_tilde_hat(k,:) = P \ X_hat(k,:)';
end


% Reconstruction
X_tilde = zeros(size(X_real));
X_tilde_hat = zeros(size(X_hat));
for i = 1:m
    X_tilde = X_tilde + ALPHA_tilde(:,i) .* X_all(:,:,i);
    X_tilde_hat = X_tilde_hat + ALPHA_tilde_hat(:,i) .* X_all(:,:,i);
end

% Ploting
figure('Position',[0,0,9e2,1.3e3])
sgtitle('System States for different Subsystems')
for idx_n = 1:n
    subplot(n,1,idx_n)
    hold on
    for i = 1:m
        plot(X_all(:,idx_n,i)', 'DisplayName', ['m = ',num2str(i)])
    end
    plot(X_real(:,idx_n), 'DisplayName', 'Actual')
    title(['x_',num2str(idx_n)])
    legend
end

% This all assumes n = m
figure('Position',[0,0,2.2e3,1.2e3])
sgtitle('Subsystem Based w/ Alpha Parameter Estimates')
num_col = 3;
for idx_n = 1:n
    subplot(n,num_col,num_col*idx_n-2)
    hold on
    plot(ALPHA_tilde(:,i), 'DisplayName', '$\tilde{\alpha}$')
    plot(ALPHA_tilde_hat(:,i), 'DisplayName', '$\hat{\tilde{\alpha}}$')
    title('Alpha Parameters')
    legend('Interpreter','latex')
    
    subplot(n,num_col,num_col*idx_n-1)
    hold on
    plot(X_real(:,idx_n), 'DisplayName', 'Actual')
    plot(X_tilde(:,idx_n), 'DisplayName', 'Alpha Based')
    title('Real vs Reconstructed State')
    legend
    
    subplot(n,num_col,num_col*idx_n)
    hold on
    plot(X_hat(:,idx_n), 'DisplayName', 'Actual')
    plot(X_tilde_hat(:,idx_n), 'DisplayName', 'Alpha Based')
    title('Real vs Reconstructed Estimate')
    legend
    
end

end