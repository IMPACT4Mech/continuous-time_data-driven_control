%% IMPACT4Mech - Continuous-Time Data-Driven Control
% Numerical example 1 of:
% A. Bosso, M. Borghesi, A. Iannelli, G. Notarstefano, A. R. Teel,
% "Data-Driven Control of Continuous-Time LTI Systems via Non-Minimal
% Realizations"

% This file requires the installation of MOSEK and YALMIP
% MOSEK:  https://docs.mosek.com/10.2/toolbox/index.html
% YALMIP: https://yalmip.github.io

%% Startup functions

clear
clc

% startup messages
disp('IMPACT4Mech - Continuous-Time Data-Driven Control')
disp('Batch reactor numerical example, version 1')

% plot settings
% fonts
font_size_s = 15;
font_size_l = 19;
% lines
line_width  = 2;
% plot position
plot_x      = 500;
plot_y      = 300;
% plot size
plot_width  = 600;
plot_height = 200;

%% Batch reactor model
% G. C. Walsh and H. Ye, “Scheduling of networked control systems”

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
     0,  1,  0,  0];

% state dimension
n    = width(A);

disp('Plant matrices:')
disp('A =')
disp(A)
disp('B =')
disp(B)
disp('C =')
disp(C)

%% Checking the assumptions

% Assumption 1
disp('Rank of the reachability matrix:')
disp(rank(ctrb(A, B)))
disp('Rank of the observability matrix:')
disp(rank(obsv(A, C)))
% Assumption 2
disp('Rank of [A B; C zeros(2)]:')
disp(rank([A B; C zeros(2)]))
% Assumption 3
disp('Observability indices:')
disp(indices(A', C')')

%% Continuous-time dataset
% dataset generated via simulation

% experiment duration
tau = 2;

% plant initial conditions
% x0  = 1*(2*rand(n, 1) - 1); % random initial conditions
x0  = [-0.1490; 0.2225; 0.7115; 0.3416]; % example for the paper

% exploration input (sum of sinusoids)
omega1 = 5;
omega2 = 7;
t      = 0:tau/1000000:tau;
u      = [1*sin(1*omega1*t) + 2*sin(2*omega1*t) +...
          3*sin(3*omega1*t) + 4*sin(4*omega1*t);
          4*sin(1*omega2*t) + 3*sin(2*omega2*t) +...
          2*sin(3*omega2*t) + 1*sin(4*omega2*t)];

% plant simulation
plant = ss(A, B, C, []);
y     = lsim(plant, u, t, x0)';

disp('Dataset initial conditions:')
disp(x0)

% plotting the dataset
figure(1)
hold on
grid on
box on
ax = gca;
ax.FontSize = font_size_s;
plot(t, u, 'LineWidth', line_width)
plot(t, y, 'LineWidth', line_width)
xlabel('$t$ [s]', FontSize=font_size_l, Interpreter='latex')
legend('$u_1$', '$u_2$', '$y_1$', '$y_2$',...
       Location='southeast', FontSize=font_size_s, Interpreter='latex')
set(gcf,'position', [plot_x, plot_y, plot_width, plot_height])

%% General tuning

% dimensions
% input dimension
m = height(u);
% output dimension
p = height(y);
% regulated output dimension
q = p;

% sampling parameters
% number of samples
N       = 50;
% sampling time
tau_s   = tau/N;
% sample points
samples = 0:tau_s:tau-tau_s;

%% Observability index estimation

disp('Computing the observability index from the dataset')

nu = find_nu(u, y, t, N); % Algorithm 3

disp('Observability index found:')
disp(nu)

%% Filter and internal model tuning

% filter tuning
Lambda  = -4*diag(1:nu);
ell     = (1:nu)';
% filter state dimension
mu      = nu*(p + m);
% filter matrices
F  = kron(eye(p + m), Lambda);
G  = [   zeros(p*nu, m);
      kron(eye(m), ell)];
L  = [kron(eye(p), ell);
         zeros(p*nu, p)];

% internal model tuning (integral action)
S_0     = 0;
Gamma_0 = 5;
% exosystem state dimension
d    = width(S_0);
% internal model matrices
Phi   = kron(eye(d*q), S_0);
Gamma = kron(eye(d*q), Gamma_0);

%% Non-minimal realization

[H, Pi] = realization(A, B, C, Lambda, ell);

disp('Parameters of the non-minimal realization')
disp('(not used for control design):')

disp('H =')
disp(H)

disp('Pi =')
disp(Pi)

disp('Dimension of the non-minimal realization:')
disp(mu)
disp('Rank of the reachability matrix:')
disp(rank(ctrb(F + L*H, G)))

%% Problem 1 - Stabilization

disp('Problem 1: data-driven stabilization')

%% 1.1 - Filtering

% initial conditions
chi0   = ell;
zeta0  = zeros(mu, 1);

% simulation
Dchi       = ss(Lambda, zeros(size(chi0)), eye(nu), []);
chi        = lsim(Dchi, zeros(1, length(t)), t, chi0)';
chi(:, 1)  = chi0;

Dzeta      = ss(F, [G L], eye(mu), []);
zeta       = lsim(Dzeta, [u; y], t, zeta0)';
zeta(:, 1) = zeta0;

% derivatives
dzeta = F*zeta + G*u + L*y;

%% 1.2 - Sampling

% input sampling
U       = interp1(t, u', samples, "nearest")';

% states sampling
X       = interp1(t, chi', samples, "nearest")';
Z       = interp1(t, zeta', samples, "nearest")';

% derivatives sampling
dZ      = interp1(t, dzeta', samples, "nearest")';

%% 1.3 - Controller matrices
% controller structure:
%
%   Dxi = A_c*xi + B_c*y
%     u = C_c*xi

% stabilizing gain
K1 = stabilize(U, X, Z, dZ);

disp('K =')
disp(K1)

% controller matrices
A_c1 = F + G*K1;
B_c1 = L;
C_c1 = K1;

%% 1.4 - Validation

% closed-loop matrix
A_cl1 = [     A   B*C_c1;
         B_c1*C     A_c1];

% stability check
disp('Closed-loop eigenvalues (stabilization):')
disp(eig(A_cl1))

%% Problem 2 - Control with integral action

disp('Problem 2: data-driven control with integral action')

%% 2.1 - Filtering

e = y; % regulated output (same data as in problem 1)

% initial conditions
chi0   = [Gamma_0; ell];
zeta0  = zeros(mu, 1);
eta0   = zeros(d*q, 1);

% simulation
Dchi      = ss(blkdiag(S_0, Lambda), zeros(size(chi0)), eye(nu + 1), []);
chi       = lsim(Dchi, zeros(1, length(t)), t, chi0)';
chi(:, 1) = chi0;

Dzeta      = ss(F, [G L], eye(mu), []);
zeta       = lsim(Dzeta, [u; e], t, zeta0)';
zeta(:, 1) = zeta0;

Deta       = ss(Phi, Gamma, eye(d*q), []);
eta        = lsim(Deta, e, t, eta0)';
eta(:, 1)  = eta0;

% derivatives
dzeta = F*zeta + G*u + L*e;
deta  = Phi*eta + Gamma*e;

%% 2.2 - Sampling

% sample points
samples = 0:tau_s:tau-tau_s;

% input sampling
U   = interp1(t, u', samples, "nearest")';

% states sampling
Z2  = interp1(t, [zeta' eta'], samples, "nearest")';
X2  = interp1(t, chi', samples, "nearest")';

% derivatives sampling
dZ2 = interp1(t, [dzeta' deta'], samples, "nearest")';

%% 2.3 - Controller matrices
% controller structure:
%
%   Dxi = A_c*xi + B_c*y
%     u = C_c*xi

% stabilizing gain
K2 = stabilize(U, X2, Z2, dZ2);

disp('K =')
disp(K2)

% controller matrices
A_c2 = blkdiag(F, Phi) + [G; zeros(d*q, m)]*K2;
B_c2 = [L; Gamma];
C_c2 = K2;

%% 2.4 - Validation

% closed-loop matrices
A_cl2 = [     A   B*C_c2;
         B_c2*C     A_c2];
B_cl2 = [zeros(n, m);
               -B_c2];
C_cl2 = [C zeros(p, mu+d*q)];

% 1 - stability check
disp('Closed-loop eigenvalues (control with integral action):')
disp(eig(A_cl2))

% 2 - regulation objective
% simulation time
tau_sim = 30;
t_sim   = 0:tau_sim/1000000:tau_sim;

% output reference
y_ref              = 3*ones(p, width(t_sim));
y_ref(2, 1:300000) = zeros(1, 300000);
y_ref(1, 1:700000) = zeros(1, 700000);

% random initial conditions
% x0_sim    = 2*(2*rand(n, 1) - 1);
% zeta0_sim = 0.1*(2*rand(mu, 1) - 1);
% eta0_sim  = 0.1*(2*rand(d*q, 1) - 1);

% example for the paper
x0_sim    = [ 1.4687; -1.1638;  0.9150;  0.4532];
zeta0_sim = [ 0.0026;  0.0919; -0.0631; -0.0599;
              0.0318;  0.0478; -0.0251;  0.0952];
eta0_sim  = [-0.0930; -0.0344];

% simulation
Cl_sys = ss(A_cl2, B_cl2, C_cl2, []);
y_cl   = lsim(Cl_sys, y_ref, t_sim, [x0_sim; zeta0_sim; eta0_sim])';

% plotting the results
figure(2)
hold on
grid on
box on
ax = gca;
ax.FontSize = font_size_s;
plot(t_sim, y_cl, 'LineWidth', line_width)
plot(t_sim, y_ref, '--', 'LineWidth', line_width)
xlabel('$t$ [s]', FontSize=font_size_l, Interpreter='latex')
legend('$y_{p1}$', '$y_{p2}$', '$y^\star_1$', '$y^\star_2$',...
       Location='southeast', FontSize=font_size_s, Interpreter='latex')
set(gcf,'position', [plot_x, plot_y, plot_width, plot_height])