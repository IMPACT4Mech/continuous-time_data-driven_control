%% IMPACT4Mech - Continuous-Time Data-Driven Control
% Numerical example for Algorithm 1 of the paper:
% A. Bosso, M. Borghesi, A. Iannelli, G. Notarstefano, A. R. Teel
% "Derivative-Free Data-Driven Control of Continuous-Time Linear
% Time-Invariant Systems." 2025 European Control Conference (ECC).

% This file requires the installation of MOSEK and YALMIP
% MOSEK:  https://docs.mosek.com/10.2/toolbox/index.html
% YALMIP: https://yalmip.github.io

%% Startup functions

clear
clc

%% System definition
% Batch reactor example, see:
% G. C. Walsh and H. Ye, “Scheduling of networked control systems”

% state dimension
n = 4;
% input dimension
m = 2;

% system matrices
A = [   1.38, -0.2077,  6.715, -5.676;
     -0.5814,   -4.29,      0,  0.675;
       1.067,   4.273, -6.654,  5.893;
       0.048,   4.273,  1.343, -2.104];
B = [    0,      0;
     5.679,      0;
     1.136, -3.146;
     1.136,      0];

C = [1,  0,  1, -1;
     0,  1,  0,  0]; % not used

%% Algorithm parameters

% experiment duration
T  = 1.5;

% filter gains
lambda = 1;
gamma  = 1;

% sampling time
Ts     = 0.1;

%% Continuous-time dataset

% plant initial conditions
% x0 = 2*rand(n, 1) - 1; % random initial conditions
x0  = [0.3110; -0.6576; 0.4121; -0.9363]; % example for the paper

% applied input (sum of sinusoids)
omega1 = 5;
omega2 = 7;
t      = 0:T/1000000:T;
u      = [    sin(omega1*t) + 2*sin(2*omega1*t) +...
          3*sin(3*omega1*t) + 4*sin(4*omega1*t);
            4*sin(omega2*t) + 3*sin(2*omega2*t) +...
          2*sin(3*omega2*t) +   sin(4*omega2*t)];

% plant simulation
plant = ss(A, B, eye(n), []);
x     = lsim(plant, u, t, x0)';

%% Filtering

% filter initial conditions
zeta10  = zeros(n, 1);
zeta20  = zeros(m, 1);

% filter dynamics
filter1 = ss(-lambda*eye(n), gamma*eye(n), eye(n), []);
zeta1   = lsim(filter1, x, t, zeta10)';

filter2 = ss(-lambda*eye(m), gamma*eye(m), eye(m), []);
zeta2   = lsim(filter2, u, t, zeta20)';

% filter derivatives
dzeta1 = -lambda*zeta1 + gamma*x;
dzeta2 = -lambda*zeta2 + gamma*u;

% filter transient error
epsilon = lsim(filter1, zeros(n, length(t)), t, x0)';

%% Sampling

% sample points
samples = 0:Ts:T-Ts;
N       = size(samples, 2);

% input sampling
U    = interp1(t, u', samples, "nearest")';

% filter states sampling
Z1   = interp1(t, zeta1', samples, "nearest")';
Z2   = interp1(t, zeta2', samples, "nearest")';
Z    = [Z1; Z2];

% filter derivatives sampling
dZ1  = interp1(t, dzeta1', samples, "nearest")';
dZ2  = interp1(t, dzeta2', samples, "nearest")';
dZ   = [dZ1; dZ2];

% error sampling
E    = interp1(t, epsilon', samples, "nearest")';
D    = [gamma*eye(n); zeros(m, n)];
DE   = D*E;

%% Plotting data

subplot(3, 1, 1)
hold on
grid on
box on
plot(t, u, 'LineWidth', 1.5)
scatter(samples, U, 40, 'filled', 'o')
title('Input')

subplot(3, 1, 2)
hold on
grid on
box on
plot(t, zeta1, 'LineWidth', 1.5)
scatter(samples, Z1, 40, 'filled', 'o')
plot(t, zeta2, 'LineWidth', 1.5)
scatter(samples, Z2, 40, 'filled', 'o')
title('Filter states')

subplot(3, 1, 3)
hold on
grid on
box on
plot(t, dzeta1, 'LineWidth', 1.5)
scatter(samples, dZ1, 40, 'filled', 'o')
plot(t, dzeta2, 'LineWidth', 1.5)
scatter(samples, dZ2, 40, 'filled', 'o')
title('Filter derivatives')

%% Computing the control gain

% decision variables
Q = sdpvar(size(Z, 2), n + m);
P = sdpvar(n + m, n + m);

% LMI constraints
Lyap_LMI = (dZ - DE)*Q + Q'*(dZ - DE)' <= -eps;
P_LMI    = P >= eps;
symmetry = Z*Q == P;

% Solving the LMI
constr = Lyap_LMI + P_LMI + symmetry;
obj = 0;
ops = sdpsettings('solver', 'mosek');
optimize(constr, obj, ops);

% gain computation
Q = value(Q);
K = U*Q*pinv(Z*Q);

%% Stability check

F = [          A              B;
     zeros(m, n) -lambda*eye(m)];
G = [zeros(n, m);  gamma*eye(m)];

disp('Eigenvalues of F+GK:')
disp(eig(F + G*K))

A_augmented = [           A        zeros(n)    zeros(n, m);
               gamma*eye(n)  -lambda*eye(n)    zeros(n, m);
                zeros(m, n)     zeros(m, n) -lambda*eye(m)];
B_augmented = [B; zeros(n, m); gamma*eye(m)];

disp('Closed-loop eigenvalues:')
disp(eig(A_augmented + B_augmented*[zeros(m, n) K]))
