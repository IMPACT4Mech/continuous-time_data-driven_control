%% IMPACT4Mech - Continuous-Time Data-Driven Control
% Numerical example for Algorithm 2 of the paper:
% A. Bosso, M. Borghesi, A. Iannelli, G. Notarstefano, A. R. Teel
% "Derivative-Free Data-Driven Control for Continuous-Time Linear
% Time-Invariant Systems"

% This file requires the installation of MOSEK and YALMIP
% MOSEK:  https://docs.mosek.com/10.2/toolbox/index.html
% YALMIP: https://yalmip.github.io

%% Startup functions

clear
clc

%% System definition

% plant transfer function
s          = tf('s');
plant_tf   = (s - 1)/(s^2 + 4)/s; % replace with the desired plant
[num, den] = tfdata(plant_tf, 'v');

% state space realization (controllability canonical form)
n = size(den, 2) - 1; % order of the system
A = [zeros(n-1, 1) eye(n-1);
          -flip(den(2:end))];
b = [zeros(n-1, 1); 1];
c = [flip(num(2:end))]';

%% Algorithm parameters

% experiment duration
T  = 2;

% filter gains
tau    = 1;
lambda = -(1/tau)*(1:n)'; % descending order
Lambda = diag(lambda);
ell    = -lambda;

% sampling time
Ts     = 0.1;

%% Continuous-time dataset

% plant initial conditions
% x0 = 5*(2*rand(n, 1) - 1); % random initial conditions
x0  = [-3.9223; 4.0631; 3.7965]; % example for the paper 

% applied input (sum of sinusoids)
omega = 5;
t     = 0:T/1000000:T;
u     =    4*sin(omega*t) + 3*sin(2*omega*t) +...
         2*sin(3*omega*t) +   1*sin(omega*t);

% plant simulation
plant = ss(A, b, c', []);
y     = lsim(plant, u, t, x0)';

%% Non-minimal realization (needed for analysis)

% finding Pi and theta such that:
% 1) Pi*(blkdiag(Lambda, Lambda) + [ell; zeros(n, 1)]*theta') = A*Pi
% 2) Pi*[zeros(n, 1); ell] = b
% 3) theta' = c'*Pi
% see the proof of Lemma 3
psi   = place(A', c, flip(lambda))';
A1    = A - psi*c';
[H,J] = eig(A1);
[DJ, indices] = sort(diag(J), 'descend');
H    = H(:, indices); % permuting the columns of H
% Pi1
mu    = pinv(ctrb(Lambda, ell))*pinv(H)*psi;
X     = zeros(n);
for i = 1:n
    X = X  + mu(i)*Lambda^(i-1);
end
Pi1   = H*X;
% Pi2
mu    = pinv(ctrb(Lambda, ell))*pinv(H)*b;
X     = zeros(n);
for i = 1:n
    X = X  + mu(i)*Lambda^(i-1);
end
Pi2   = H*X;
Pi    = [Pi1 Pi2];

% theta
theta1 = Pi1' * c;
theta2 = Pi2' * c;
theta  = [theta1; theta2];

% nonminimal realization
F = [Lambda + ell*theta1' ell*theta2';
                 zeros(n)      Lambda];
g = [zeros(n, 1); ell];

%% Filtering

% filter initial conditions
zeta10 = zeros(n, 1);
zeta20 = zeros(n, 1);
chi0   = ones(n, 1);

% filter dynamics
filter = ss(Lambda, ell, eye(n), []);
zeta1  = lsim(filter, y, t, zeta10)';
zeta2  = lsim(filter, u, t, zeta20)';
chi    = lsim(filter, zeros(1, length(t)), t, chi0)';

% derivatives
dzeta1 = Lambda*zeta1 + ell*y;
dzeta2 = Lambda*zeta2 + ell*u;
dchi   = Lambda*chi;

%% Sampling

% sample points
samples = 0:Ts:T-Ts;
N       = size(samples, 2);

% input sampling
U    = interp1(t, u, samples, "nearest");

% filter states sampling
Z1   = interp1(t, zeta1', samples, "nearest")';
Z2   = interp1(t, zeta2', samples, "nearest")';
Chi  = interp1(t, chi', samples, "nearest")';
Za   = [Chi; Z1; Z2];

% filter derivatives sampling
dZ1  = interp1(t, dzeta1', samples, "nearest")';
dZ2  = interp1(t, dzeta2', samples, "nearest")';
dChi = interp1(t, dchi', samples, "nearest")';
dZa  = [dChi; dZ1; dZ2];

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
Q = sdpvar(size(Za, 2), 3*n);
P = sdpvar(3*n, 3*n);

% LMI constraints
Lyap_LMI = dZa*Q + Q'*dZa' <= -eps;
P_LMI    = P >= eps;
symmetry = Za*Q == P;

% Solving the LMI
constr = Lyap_LMI + P_LMI + symmetry;
obj = 0;
ops = sdpsettings('solver', 'mosek');
optimize(constr, obj, ops);

% gain computation
Q = value(Q);
K = U*Q*pinv(Za*Q)*[zeros(n, 2*n); eye(2*n)];

%% Stability check
disp('Eigenvalues of F+gK:')
disp(eig(F + g*K))

A_augmented = [       A  zeros(n)  zeros(n);
                 ell*c'    Lambda  zeros(n);
               zeros(n)  zeros(n)    Lambda];
B_augmented = [b; zeros(n, 1); ell];

disp('Closed-loop eigenvalues:')
disp(eig(A_augmented + B_augmented*[zeros(1, n) K]))