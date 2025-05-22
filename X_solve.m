function X = X_solve(Theta, beta, phi)
%X_SOLVE Solves the matrix equation X*Theta = Theta*X, X*beta = phi
%   This function returns the unique solution X when the pair (Theta, beta)
%   is controllable.
%
%   Syntax
%   X = X_SOLVE(Theta, beta, phi)
%
%   Input Arguments
%       Theta - matrix of dimension rxr
%        beta - matrix of dimension rx1
%         phi - matrix of dimension rx1

r = height(beta);
R = ctrb(Theta, beta);

if rank(R) == r
    rho = (R^-1)*phi;
    X   = zeros(r);
    for i = 1:r
        X = X  + rho(i)*Theta^(i-1);
    end
else
    error('The pair (Theta, beta) must be controllable.')
end