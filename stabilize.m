function [K, P, Q] = stabilize(U, X, Z, dZ)
%STABILIZE - Compute a stabilizing gain from a dataset
%   This function runs the following algorithm:
%
%   1) Find P and Q such that:
%               P > 0
%   dZ*Q + Q'*dZ' < 0
%             X*Q = 0
%             Z*Q = P
%
%   2) Compute the control gain:
%   K = U*Q*P^-1
%
%   Syntax
%           K = STABILIZE(U, X, Z, dZ)
%   [K, P, Q] = STABILIZE(U, X, Z, dZ)
%
%   Input Arguments
%   U - Input dataset
%     m x N matrix
%   X - Disturbance dataset
%     delta x N matrix
%   Z - State dataset
%     mu x N matrix
%   dZ - State derivative dataset
%     mu x N matrix

% decision variables
P = sdpvar(size(Z, 1), size(Z, 1));
Q = sdpvar(size(Z, 2), size(Z, 1));

% LMI constraints
P_LMI    =             P >= eps;
Lyap_LMI = dZ*Q + Q'*dZ' <= -eps;
Data_EQ1 =           X*Q == 0;
Data_EQ2 =           Z*Q == P;

% Solving the LMI
constraints = P_LMI + Lyap_LMI + Data_EQ1 + Data_EQ2;
objective   = 0;
ops         = sdpsettings('solver', 'mosek');
optimize(constraints, objective, ops);

% gain computation
P = value(P);
Q = value(Q);
K = U*Q*P^-1;

end