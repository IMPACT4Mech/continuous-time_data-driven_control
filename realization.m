function [H, Pi] = realization(A, B, C, Lambda, ell)
%REALIZATION - Canonical non-minimal realization
%   This function computes the canonical non-minimal realization of the
%   following plant (with dim(x) = n, dim(u) = m, dim(y) = p):
%
%   dx = A*x + B*u
%    y = C*x
%
%   In particular, assuming that all the observability indices of (C, A) 
%   are the same (= nu) and given a controllable pair (Lambda, ell), with
%   height(ell) = nu, define the matrices:
%
%   F =  kron(eye(2*nu), Lambda)
%   G = [         zeros(p*nu, m);
%             kron(eye(nu), ell)]
%   L = [     kron(eye(nu), ell);
%                 zeros(p*nu, p)]
% 
%   Then, REALIZATION returns the matrices Pi and H such that the system:
%
%   dzeta = (F + L*H)*zeta + G*u
%       y = H*zeta
%
%   is input-output equivalent to the original plant, where zeta and x are
%   related via x = Pi*zeta.
%
%   Syntax
%         H = REALIZATION(A, B, C, Lambda, ell)
%   [H, Pi] = REALIZATION(A, B, C, Lambda, ell)
%
%   Input Arguments
%        A - Matrix of dimension n x n
%          state matrix of the plant
%        B - Matrix of dimension n x m
%          input matrix of the plant
%        C - Matrix of dimension p x n
%          output matrix of the plant
%   Lambda - Matrix of dimension nu x nu
%          tuning matrix
%      ell - Matrix of dimension nu x 1
%          tuning matrix

% dimensions
n    = width(A);        % state dimension
m    = width(B);        % input dimension
p    = height(C);       % output dimension
nu_i = indices(A', C'); % observability indices
nu   = max(nu_i);       % observability index

% stops the function if the arguments do not satisfy the requirements
for i = 1:p
    if nu_i(i) ~= nu
        error('The observability indices must be all equal.')
    end
end
Q = ctrb(Lambda, ell);
if (height(ell) ~= nu) || (rank(Q) ~= nu)
    error('Incorrect tuning.')
end

% observer canonical form
[A_o, ~, T_o] = canonical(A', C');
A_o           = A_o';
T_o           = (T_o')^-1;
B_o           = T_o*B;
A_m           = A_o(:, cumsum(nu_i));

% construction of Psi
Lambda_o  = Q^-1*Lambda*Q;
th_lambda = -Lambda_o(:, end);
Psi       = A_m + kron(eye(p), th_lambda);

% computation of Y = T_o*Pi
Y_y = zeros(n, nu*p); % preallocation
for i = 1:p
    for j = 1:p
        Y_y(nu*(i-1)+1:nu*i, nu*(j-1)+1:nu*j) = ...
            (Q^-1)*X_solve(Lambda, ell, Q*Psi(nu*(i-1)+1:nu*i, j));
    end
end
Y_u = zeros(n, nu*m); % preallocation
for i = 1:p
    for j = 1:m
        Y_u(nu*(i-1)+1:nu*i, nu*(j-1)+1:nu*j) = ...
            (Q^-1)*X_solve(Lambda, ell, Q*B_o(nu*(i-1)+1:nu*i, j));
    end
end
Y  = [Y_y Y_u];

% outputs
Pi = (T_o^-1)*Y;
H  = C*Pi;

end