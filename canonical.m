function [A_c, B_c, P] = canonical(A, B)
%CANONICAL - Compute the multivariable controller canonical form
%   This function returns the controller canonical form of the pair (A, B)
%   in the general multivariable scenario. The function also returns the 
%   associated similarity transformation P such that:
%
%   A_c = P*A*P^-1
%   B_c = P*B
%
%   Syntax
%      [A_c, B_c] = CANONICAL(A, B)
%   [A_c, B_c, P] = CANONICAL(A, B)
%
%   To compute the observer canonical form of a pair (C, A), run:
%   [A_o, C_o, P] = CANONICAL(A', C');
%   A_o           = A_o';
%   C_o           = C_o';
%   P             = (P')^-1;
%
%   Input Arguments
%   A - Matrix of dimension n x n
%     state matrix
%   B - Matrix of dimension n x m
%     input matrix

n       = width(A);        % state dimension
m       = width(B);        % input dimension
[mu, T] = indices(A, B);   % controllability indices

sigma   = cumsum(mu);
T_inv   = T^-1;
Q       = T_inv(sigma, :); % extracting m rows of T_inv

P       = zeros(n);        % preallocation
for i = 1:m                % i: row of Q
    for j = 1:mu(i)        % j-1: power of A
        P(sigma(i) - mu(i) + j, :) = Q(i, :)*A^(j - 1);
    end
end

A_c = P*A*P^-1;
B_c = P*B;

end