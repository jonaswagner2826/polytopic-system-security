% Uncertain Observer Toy System Simulation
% Jonas Wagner
% 2021-08-11 @ 6:12PM

clear
close all

% Sim Parameters
N = 20;
x_0_mag = 1;
x_0_theta = pi/4;
x_0 = x_0_mag * [cos(x_0_theta); sin(x_0_theta)];
x_hat_0 = x_0;

% Optional Tests
testL = false;

% % Toy System Def
% A(:,:,1) = [-0.80, 0.25; 0.25,-0.30];
% A(:,:,2) = [ 0.30, 0.70; 0.70, 0.00];
% A(:,:,3) = [-0.30, 0.65; 0.55, 0.10];
% A(:,:,4) = [ 0.55,-0.20;-0.40,-0.30];
% 
% B = [1.5; -0.5];
% C = [1, 0];
% D = 0;

% Toy Sys Def 2
V = [0.2785, 0.9575;  0.5469, 0.9649]; % generated randomly
A(:,:,1) = V * diag([0.8,0.9]) * inv(V);
A(:,:,2) = V * diag([0.85,-0.95]) * inv(V);
A(:,:,3) = V * diag([-0.85,0.95]) * inv(V);
A(:,:,4) = V * diag([-0.9,-0.8]) * inv(V);

B = 0;
C = [0.2, 0.5];
D = 0;

% System Dimensions
n = size(A,1);
m = size(A,3);
p = size(B,2);
q = size(C,1);

% Parameter Definition
Alpha_real = [1,0,0,0];%[0.25, 0.25, 0.25, 0.25];
% Alpha_hat = [0,0,0,1];%[0.375, 0.125, 0.125, 0.125];
Alpha_hat = [0.999,0.001,0,0]


% Polytopic Calculations
A_real = zeros(n);
A_hat = zeros(n);
for i = 1:m
    A_real = A_real + Alpha_real(i) * A(:,:,i);
    A_hat = A_hat + Alpha_hat(i) * A(:,:,i);
end
Delta_A = A_real - A_hat
delta = norm(Delta_A)

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

% Sim Setup
X = zeros(N,n);
X_hat = zeros(N,n);
U = zeros(N,p);
Y = zeros(N,q);

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
    
    % Save Data
    X(k,:) = x;
    X_hat(k,:) = x_hat;
    U(k,:) = u;
    Y(k,:) = y;
    E(k,:) = e;
    
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


% Ploting
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

% Residual Plot
figure
R = zeros(N,1);
for i = 1:N
    R(i,:) = C*E(i,:)';
end
plot(R)
title('Residual')