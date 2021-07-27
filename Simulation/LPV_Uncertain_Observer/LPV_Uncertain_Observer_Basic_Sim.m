% Uncertain LPV Observer Toy System Simulation
% Jonas Wagner
% 2021-07-26 @ 6:10 PM

clear
close all

% Sim Parameters
N = 10;
Alpha_real = [0.125, 0.25, 0.375, 0.25];
Alpha_hat = [0.25, 0.375, 0.25, 0.125];
x_0 = [-0.25; 0.25];
x_hat_0 = x_0;

% Toy System Def
A(:,:,1) = [-0.80, 0.25; 0.25,-0.30];
A(:,:,2) = [ 0.30, 0.70; 0.70, 0.00];
A(:,:,3) = [-0.30, 0.65; 0.55, 0.10];
A(:,:,4) = [ 0.55,-0.20;-0.40,-0.30];

B = [1.5; -0.5];
C = [1, 0];
D = 0;

% System Dimensions
n = size(A,1);
m = size(A,3);
p = size(B,2);
q = size(C,1);

% Polytopic Calculations
A_real = zeros(n);
A_hat = zeros(n);
for i = 1:m
    A_real = A_real + Alpha_real(i) * A(:,:,i);
    A_hat = A_hat + Alpha_hat(i) * A(:,:,i);
end

% Control Design: u = K * x_hat
p_ctrl = [0.5, 0.2];
K = -place(A_hat, B, p_ctrl); % Arbritrarily designed...

% Observer Design: x_hat = A_hat * x_hat + B * u + L * (y - y_hat)
p_obsv = [-0.5, 0.5];
L = place(A_hat', C', p_obsv).';

% Sim Setup
X = zeros(1,n,N);
X_hat = zeros(1,n,N);
x = x_0;
x_hat = x_hat_0;
Y = zeros(1,q,N);
Y_hat = zeros(1,q,N);
R = zeros(1,q,N);
r = 0;
U = zeros(1,p,N);
u = K * x_hat;
E = zeros(1,n,N);
e = x - x_hat;

% Simulation
for k = 1:N
    % Measurment
    y = C * x + D * u;
    y_hat = C * x_hat + D * u;
    r = y - y_hat;
    
    % Saving for t = k
    X(:,:,k) = x;
    X_hat(:,:,k) = x_hat;
    Y(:,:,k) = y;
    Y_hat(:,:,k) = y_hat;
    R(:,:,k) = r;
    U(:,:,k) = u;
    E(:,:,k) = e;
    
    % Time-Update
    u = K * x_hat;
    x = A_real * x + B * u;
    x_hat = A_hat * x_hat + B * u + L * r;
    e = x - x_hat;
end


% Ploting
X1 = reshape([X(:,1,:); X_hat(:,1,:)],2,[])';
X2 = reshape([X(:,2,:); X_hat(:,2,:)],2,[])';
Y1 = reshape([Y(:,1,:); Y_hat(:,1,:)],2,[])';
U = reshape(U, p, []);
E = reshape(E,n,[]);
R = reshape(R,q,[]);

figure
subplot(2,2,1)
plot(X1)
hold on
plot(X2)
title('States and Estimates')

subplot(2,2,2)
plot(U)
title('Inputs');

subplot(2,2,3)
plot(reshape(E,2,[])')
title('Error')

subplot(2,2,4)
plot(R)
title('Residual')

