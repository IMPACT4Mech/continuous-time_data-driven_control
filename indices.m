function [mu, T] = indices(A, B)
%INDICES - Controllability indices of (A, B)
%   This function returns the controllability indices of the pair (A, B)
%   and a non-singular square matrix that collects the first n linearly
%   independent vectors of the reachability matrix, appropriately
%   reordered according to the columns of B.
%
%   Syntax
%        mu = INDICES(A, B)
%   [mu, T] = INDICES(A, B)
%
%   Input Arguments
%   A - Matrix of dimension n x n
%     state matrix
%   B - Matrix of dimension n x m
%     input matrix

R  = ctrb(A, B);  % reachability matrix
n  = width(A);    % state dimension
m  = width(B);    % input dimension

if rank(R) ~= n
    error('Pair (A, B) is uncontrollable.')
end

if rank(B) ~= m
    error('rank(B) is not equal to m.')
end

% computing the controllability indices
mu = zeros(m, 1); % preallocation
T  = zeros(n);    % preallocation
k  = 1;           % column # of T
for i = 0:n-1              % i: power of A
    for j = 1:m            % j: column of B
        v = R(:, i*m + j); % v = A^i B(:, j)
        if rank([T v]) == rank(T) + 1
            % if v is linearly independent from
            % the previous columns of R:
            T(:, k) = v;         % update T
            k       = k + 1;     % next column
            mu(j)   = mu(j) + 1; % update the index
        end
    end
end

% reordering the linearly independent columns
id = cumsum(mu) - (mu - 1); % first column of each index
for i = 1:m                 % i: column of B
    for j = 0:mu(i)-1       % j: power of A
        T(:, id(i) + j) = A^j*B(:, i);
    end
end

end
