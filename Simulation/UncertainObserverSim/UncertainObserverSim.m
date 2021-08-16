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
Alphas = eye(4);

% calc1 = false;
% calc2 = true;


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
A(:,:,3) = V * diag([0.85,-0.95]) * inv(V);
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
Alpha_real = [0.25, 0.25, 0.25, 0.25];
Alpha_hat = [0.375, 0.125, 0.125, 0.125];


% Polytopic Calculations
A_real = zeros(n);
A_hat = zeros(n);
for i = 1:m
    A_real = A_real + Alpha_real(i) * A(:,:,i);
    A_hat = A_hat + Alpha_hat(i) * A(:,:,i);
end
Delta_A = A_real - A_hat;

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
cvx_begin sdp
    variable X(n,q)
    variable Q(n,n) symmetric
    tol = 1e-6;
    maximize(trace(Q))
    subject to
        for i = 1:m
            Q >= tol*eye(n);
            [Q, (Q*A(:,:,i) - X*C)';
             Q*A(:,:,i) - X*C, Q] >= tol*eye(2*n);
        end
cvx_end

L = inv(Q) * X

% % Testing L
% for i = 1:m
%     A(:,:,i)
%     A_LC = A(:,:,i) - L * C
%     eig_A_LC = eig(A_LC)
% end
% A_real
% A_LC = A_real - L*C
% eig_A_LC = eig(A_LC)
% A_hat
% A_LC = A_hat - L*C
% eig_A_LC = eig(A_LC)

% Sim Setup
X = zeros(N,n);
X_hat = zeros(N,n);
% U = zeros(N,p);
Y = zeros(N,q);

E = zeros(N,n);
E_calc = zeros(N,n);
E_norm = zeros(N,1);
E_norm_calc = zeros(N,1);

% Inital setup (k=0)
x = x_0;
x_hat = x_hat_0;
e = x - x_hat;

% Simulation
for k = 1:N
    % Control
    u = 0; % No Input
    
    % Time Update
    x = A_real * x + B * u;
    y = C * x + D * u;
    x_hat = A_hat * x_hat + B * u + L * (y - C * x_hat);
    e = x - x_hat;
    
    
    % Save Data
    X(k,:) = x;
    X_hat(k,:) = x_hat;
%     U(k,:) = u;
    Y(k,:) = y;
    E(k,:) = e;
    
%     % Error Calc Estimate
%     e_calc = x_0 - x_hat_0;
%     for i = 0:(k-1)
%         e_calc = e_calc + (A_hat - L * C)^(i) * Delta_A * A_real^(k - i) * x_0;
%     end
% %     e_calc = e_calc + (A_hat - L * C)^k * x_0;
%     E_calc(k,:) = e_calc;
    
    % New error calc attempt
    x_hat_calc = x_hat_0;
    x_hat_calc =  L*C*A_real^k *x_0;
    if k > 1
        for i = 1:(k-1)
            x_hat_calc = (A_hat - L*C)^i * L*C*(A_real)^(k-1-i) * x_0;
        end
    end
    x_hat_calc = (A_hat - L*C)^k * x_0;
    
    e_calc = A_real^k * x_0 - x_hat_calc;
    E_calc(k,:) = e_calc;
        
    
%     % Error Norm
%     E_norm(k,1) = norm(x - x_hat);  
%     
%     %Error Norm Calc
%     e_norm_calc = 0;
%     for i = 0:k
%         if calc1
%             e_norm_calc = e_norm_calc + norm((A_hat - L * C)^(i) * Delta_A * A_real^(k - i)) * norm(x_0);
%         elseif calc2
%             e_norm_calc = e_norm_calc + (norm((A_hat - L * C))^(i) * norm(Delta_A) * norm(A_real)^(k - i)) * norm(x_0);
%         end
%     end
%     E_norm_calc(k,:) = e_norm_calc;
end


% Ploting
figure
sgtitle(num2str(x_0))

subplot(2,2,1)
plot(X(:,1));
hold on
plot(X_hat(:,1));
title('State and Estimate (X1)')

subplot(2,2,3)
plot(X(:,2));
hold on
plot(X_hat(:,2));
title('State and Estimate (X2)')

subplot(2,2,2)
plot(E(:,1))
hold on
plot(E_calc(:,1))
title('Error (X1)')
legend('actual','calc')
% plot(E_norm)
% hold on
% plot(E_norm_calc)
% if calc1
%     title('Error Norm Calc 1')
% elseif calc2
%     title('Error Norm Calc 2')
% end
% % title('Error Norm')
% legend('actual','calc')

subplot(2,2,4)
plot(E(:,2))
hold on
plot(E_calc(:,2))
title('Error (X2)')
legend('actual','calc')
% plot((E_norm - E_norm_calc)./E_norm)
% title('Error Norm Calc % Error')



% % Sim Setup
% X = zeros(1,n,N);
% X_hat = zeros(1,n,N);
% x = x_0;
% x_hat = x_hat_0;
% Y = zeros(1,q,N);
% Y_hat = zeros(1,q,N);
% R = zeros(1,q,N);
% r = 0;
% 
% U = zeros(1,p,N);
% u_0 = [0];
% f = 1/(N/size(Alphas,1));
% % U = reshape(u_0 * 0.5 * (sin(2*pi*f*(1:N))+1), size(u_0,1),p,[]);
% U_ref = reshape(u_0 * 0.5 * (square(2*pi*f*(1:N)) + 1), 1,p,[]);
% % u = K * x_hat;
% E = zeros(1,n,N);
% e = x - x_hat;
% u = 0;%U_ref(:,1);
% % Alpha_real = ALPHA(:,1);
% j = 0;
% E_star = zeros(1,N);
% E_star_rec = zeros(2,N);
% e_star_rec = e;
% 
% % Simulation
% for k = 1:N
%     % Measurment
%     y = C * x + D * u;
%     y_hat = C * x_hat + D * u;
%     r = y - y_hat;
%     
%     % Saving for t = k
%     X(:,:,k) = x;
%     X_hat(:,:,k) = x_hat;
%     Y(:,:,k) = y;
%     Y_hat(:,:,k) = y_hat;
%     R(:,:,k) = r;(
%     U(:,:,k) = u;
%     E(:,:,k) = e;
%     E_star(:,k) = norm(Delta_A) * norm(x_0) * ...
%         ((norm(A_hat - L * C))^k - (norm(A_hat + Delta_A))^k) /...
%         ((norm(A_hat - L * C)) - (norm(A_hat + Delta_A)));
%     E_star_rec(:,k) = e_star_rec;
% %     %switching alphas
% %     if ALPHA(:,k) ~= Alpha_real
% %         Alpha_real = ALPHA(:,k);
% %         A_real = zeros(n);
% %         for i = 1:m
% %             A_real = A_real + Alpha_real(i) * A(:,:,i);
% %         end
% % %         x = x_0;
% %     end
%     
%     % Time-Update
%     u = 0;
% %     u = U_ref(:,k) + K * x_hat;
% %     u = K * (x_hat + x_0 * U_ref(:,k));
%     x = A_real * x + B * u;
%     e = Delta_A * x + (A_hat + L * C) * e;
% %     x_hat = A_hat * x_hat + B * u + L * r;
% %     e = x - x_hat;
%     x_hat = x - e;
%     
%     %this asssumes that norm(x_0) = 1
%     e_star_rec = e_star_rec + ((A_hat - L*C)^(k-1+1) * Delta_A * (A_hat + Delta_A)^(N-k+1))*x_0;
%     
% end
% 
% 
% % Ploting
% X_real = reshape(X, n, [])';
% X_est = reshape(X_hat, n, [])';
% X1 = reshape([X(:,1,:); X_hat(:,1,:)],2,[])';
% X2 = reshape([X(:,2,:); X_hat(:,2,:)],2,[])';
% Y1 = reshape([Y(:,1,:); Y_hat(:,1,:)],2,[])';
% U = reshape(U, p, [])';
% E = reshape(E,n,[])';
% E_star_rec = reshape(E_star_rec,n,[])';
% R = reshape(R,q,[])';
% % Alpha = reshape(ALPHA, m, [])';
% 
% 
% % Unforced Plot Comparisson
% figure
% sgtitle(num2str(Delta_A))
% 
% subplot(2,2,1)
% plot(X1)
% title('X1 and X1_hat')
% 
% subplot(2,2,3)
% plot(X2)
% title('X2 and x2_hat')
% 
% subplot(2,2,2)
% plot(E(:,1,:))
% hold on
% % plot(E_star)
% plot(E_star_rec(:,1,:))
% title('Error Norm')
% % 
% % subplot(2,2,4)
% % plot(vecnorm(R,2,2))
% % hold on
% % % plot(norm(C) * E_star)
% % plot(norm(C) * E_star_rec)
% % title('Residual Norm')
% 
% 
% 
% % 4-plots
% % figure
% % subplot(2,2,1)
% % plot(X1)
% % hold on
% % plot(X2
% % title('States and Estimates')
% % 
% % subplot(2,2,2)
% % plot(U)
% % title('Inputs');
% % 
% % subplot(2,2,3)
% % plot(reshape(E,2,[])')
% % title('Error')
% % 
% % subplot(2,2,4)
% % plot(R)
% % title('Residual')
% 
% % % 6-plots
% % figure
% % subplot(3,2,1)
% % plot(X_real)
% % title('States')
% % 
% % subplot(3,2,3)
% % plot(X_est)
% % title('Estimates')
% % 
% % subplot(3,2,5)
% % plot(E)
% % title('Error')
% % 
% % subplot(3,2,2)
% % plot(Alpha)
% % title('Alpha')
% % legend('x1','x2','x3','x4')
% % ylim([-0.1, 1.1])
% % 
% % subplot(3,2,4)
% % plot(U)
% % title('Input')
% % 
% % subplot(3,2,6)
% % plot(R)
% % title('Residual')