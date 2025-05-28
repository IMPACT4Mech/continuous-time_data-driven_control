%% IMPACT4Mech - Continuous-Time Data-Driven Control
% Numerical example 2 of the paper:
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
disp('Surface vessel numerical example, version v1')

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

%% Surface vessel model
% A. Pyrkin and A. Isidori, "Adaptive output regulation of right-invertible
% MIMO LTI systems, with application to vessel motion control"

% system matrices
A = [-0.1   0.012  0.015  0  0   0.01;
     0.01 -0.0333  -0.05  0  0 -0.014;
     0.02    0.03  -0.18  0  0      0;
        1       0      0  0  0      0;
        0       1      0  0  0      0;
        0       0      1  0  0      0];
B = [  0 0.03 0.025;
       0 0.21  -0.2;
     0.1 0.03  0.02;
       0    0     0;
       0    0     0;
       0    0     0];
C = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

% exosystem
S_mod = [0       0        0;
         0       0 -2*pi/10;
         0 2*pi/10        0];
T     = ctrb(S_mod, [1; 0; 1]);
S     = (T)'*S_mod'*(T^-1)';

% disturbance model
P = [-0.001    0 0.002;
       0.02 0.01 -0.02;
          0    0     0;
          0    0     0;
        0.1    0     0;
        0.1  0.1  -0.1];
Q = [1 0  25/pi^2;
     0 0        0;
     0 0        0]*2;

% regulated output
C_e = [eye(2) zeros(2, 1)]*C;
Q_e = [eye(2) zeros(2, 1)]*Q;

% state dimension
n    = width(A);
% input dimension
m    = width(B);
% output dimension
p    = height(C);
% regulated output dimension
q    = height(C_e);
% observability indices
nu_i = indices(A', C');
nu   = max(nu_i);
% exosystem state dimension
d    = width(S);

disp('Plant matrices:')
disp('A =')
disp(A)
disp('B =')
disp(B)
disp('C =')
disp(C)
disp('S =')
disp(S)
disp('P =')
disp(P)
disp('Q =')
disp(Q)

%% Checking the assumptions
% assumption 1
disp('Rank of the reachability matrix:')
disp(rank(ctrb(A, B)))
disp('Rank of the observability matrix:')
disp(rank(obsv(A, C)))
% assumption 2
disp('Rank of [A B; C_e zeros(2, 3)]:')
disp(rank([A B; C_e zeros(2, 3)]))
disp('Rank of [A-2*pi*0.1i*eye(n) B; C_e zeros(2, 3)]:')
disp(rank([A-2*pi*0.1i*eye(n) B; C_e zeros(2, 3)]))
disp('Rank of [A+2*pi*0.1i*eye(n) B; C_e zeros(2, 3)]:')
disp(rank([A+2*pi*0.1i*eye(n) B; C_e zeros(2, 3)]))
% assumption 3
disp('Observability indices:')
disp(indices(A', C')')

%% Continuous-time dataset
% dataset generated via simulation

% experiment duration
tau = 35;

% plant initial conditions
% x0 = 1*(2*rand(n, 1) - 1); % random initial conditions
% example for the paper
x0 = [0.7297; -0.7195; 0.3143; 0.186; 0.0267; -0.5108];
w0 = [1; 1; 1];

% exploration input (sum of sinusoids)
omega1 = 1;
omega2 = 2;
omega3 = 3;
t      = 0:tau/1000000:tau;
u      = [-1*sin(1*omega1*t) + 2*sin(2*omega1*t) +...
           3*sin(3*omega1*t) - 4*sin(4*omega1*t);
           4*sin(1*omega2*t) - 3*sin(2*omega2*t) +...
          -2*sin(3*omega2*t) + 1*sin(4*omega2*t);
          -2*sin(1*omega3*t) + 2*sin(2*omega3*t) +...
           2*sin(3*omega3*t) - 2*sin(4*omega3*t)];

% exosystem simulation
exo = ss(S, zeros(d, 1), eye(d), []);
w   = lsim(exo, zeros(size(t)), t, w0)';

% plant simulation
plant = ss(A, [B P], C, [zeros(p, m) Q]);
y     = lsim(plant, [u; w], t, x0)';
e     = [eye(2) zeros(2, 1)]*y;

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
set(gcf,'position', [plot_x, plot_y, plot_width, plot_height])

%% Tuning

% sampling parameters
% number of samples
N       = 80;
% sampling time
tau_s   = tau/N;
% sample points
samples = 0:tau_s:tau-tau_s;

% filter gains
Lambda = [ 0    1;
          -2   -2];
ell    = [ 0; 0.5];

% auxiliary system gains
F_0    = Lambda;
G_0    = ell;

% filter state dimension
mu = nu*(p + m);

% filter matrices
F  = kron(eye(p + m), Lambda);
G  = [    zeros(p*nu, m);
      kron(eye(m), ell)];
L  = [kron(eye(p), ell);
          zeros(p*nu, p)];

% internal model matrices
S_0     = S;
Gamma_0 = [0; 0; 0.1];

Phi   = kron(eye(q), S);
Gamma = kron(eye(q), Gamma_0);

%% Non-minimal realization

[H, Pi] = realization(A, B, C, Lambda, ell);

disp('Parameters of the non-minimal realization:')
disp('(not used for control design):')

disp('H =')
disp(H)

disp('Pi =')
disp(Pi)

disp('Dimension of the non-minimal realization:')
disp(mu)
disp('Rank of the reachability matrix:')
disp(rank(ctrb(F + L*H, G)))

%% Filtering

% initial conditions
chi0   = [Gamma_0; G_0];
zeta0  = zeros(mu, 1);
eta0   = zeros(d*q, 1);

% simulation
Dchi       = ss(blkdiag(S_0, F_0), zeros(d + nu, 1), eye(d + nu), []);
chi        = lsim(Dchi, zeros(size(t)), t, chi0)';
chi(:, 1)  = chi0;

Dzeta      = ss(F, [G L], eye(mu), []);
zeta       = lsim(Dzeta, [u; y], t, zeta0)';
zeta(:, 1) = zeta0;

Deta       = ss(Phi, Gamma, eye(d*q), []);
eta        = lsim(Deta, e, t, eta0)';
eta(:, 1)  = eta0;

% derivatives
dzeta = F*zeta + G*u + L*y;
deta  = Phi*eta + Gamma*e;

%% Sampling

% input sampling
U  = interp1(t, u', samples, "nearest")';

% states sampling
Z  = interp1(t, [zeta' eta'], samples, "nearest")';
X  = interp1(t, chi', samples, "nearest")';

% derivatives sampling
dZ = interp1(t, [dzeta' deta'], samples, "nearest")';

%% Controller matrices
% controller structure:
%
%   Dxi = A_c*xi + B_c*y
%     u = C_c*xi

% stabilizing gain
K = stabilize(U, X, Z, dZ);

disp('K =')
disp(K)

% controller matrices
A_c = blkdiag(F, Phi) + [G; zeros(q*d, m)]*K;
B_c = [L; [Gamma zeros(d*q, p-q)]];
C_c = K;

%% Validation

% closed-loop matrix
A_cl = [    A   B*C_c;
        B_c*C     A_c];
B_cl = [P; B_c*Q];
C_cl = [C zeros(p, mu + d*q)];
D_cl = Q;

% 1 - stability check
disp('Closed-loop eigenvalues:')
disp(eig(A_cl))

% 2 - regulation objective
% simulation time
t_sim      = 0:100/1000000:100;

% initial conditions for testing
w0_sim  = [1; -3; 0];

% random initial conditions
% x0_sim    = 5*(2*rand(n, 1) - 1);
% zeta0_sim = 1*(2*rand(mu, 1) - 1);
% eta0_sim  = 1*(2*rand(d*q, 1) - 1);

% example for the paper
x0_sim    = [-0.3950; -3.4733; -2.4421; -1.0835;  4.6461; 0.1942];
zeta0_sim = [-0.7582;  0.8298; -0.0868;  0.8094; -0.4162; 0.5644;
             -0.8882; -0.5739;  0.5531;  0.3936;  0.8158; 0.2178];
eta0_sim  = [ 0.0891;  0.5062;  0.7545;  0.7875;  0.7117; 0.3665];

% simulation
w      = lsim(exo, zeros(size(t_sim)), t_sim, w0_sim)';

Cl_sys = ss(A_cl, B_cl, C_cl, D_cl);
y_cl   = lsim(Cl_sys, w, t_sim, [x0_sim; zeta0_sim; eta0_sim])';

% plotting the results
figure(2)
hold on
grid on
box on
ax = gca;
ax.FontSize = font_size_s;
plot(t_sim, y_cl, 'LineWidth', line_width)
xlabel('$t$ [s]', FontSize=font_size_l, Interpreter='latex')
legend('$e_1$', '$e_2$', '$y_r$',...
       Location='southeast', FontSize=font_size_s, Interpreter='latex')
set(gcf,'position', [plot_x, plot_y, plot_width, plot_height])