function nu = find_nu(u, y, t, N)
%FIND_NU Observability index estimation from a dataset
%   This function uses the input-output trajectory of a MIMO system to
%   compute its observability index. The algorithm assumes that the index
%   is uniform across all outputs and that the data are informative.
%
%   Syntax
%       nu = FIND_NU(u, y, t, N)
%
%   Input Arguments
%   u - Input vector
%   y - Output vector
%   t - Time vector
%   N - Number of samples

% dimensions
m = height(u); % input dimension
p = height(y); % output dimension

% sampling parameters
tau     = t(end);
tau_s   = tau/N;
samples = 0:tau_s:tau-tau_s;

% algorithm tuning
nu_max = 10;
lambda = 1:nu_max;
gamma  = 1:nu_max;

% algorithm initialization
nu_found = false;
nu_hat   = 1;
tol      = 1e-10;

% running the algorithm
% first filter:
% auxiliary dynamics
aux        = ss(-lambda(1), 0, 1, []);
chi        = lsim(aux, zeros(1, length(t)), t, gamma(1))';
X          = interp1(t, chi', samples, "nearest");
% filter dynamics
filter     = ss(-lambda(1)*eye(p+m), gamma(1)*eye(p+m), eye(p+m), []);
zeta       = lsim(filter, [y; u], t, zeros(p+m, 1))';
Z          = interp1(t, zeta', samples, "nearest")';

while (nu_found == false) && (nu_hat <= nu_max)
    nu_hat = nu_hat + 1;
    
    % filter # nu_hat:
    % auxiliary dynamics
    aux        = ss(-lambda(nu_hat), 0, 1, []);
    chi        = lsim(aux, zeros(1, length(t)), t, gamma(nu_hat))';
    chi(:, 1)  = gamma(nu_hat);
    X_new      = interp1(t, chi', samples, "nearest");
    % filter dynamics
    filter     = ss(-lambda(nu_hat)*eye(p+m), gamma(nu_hat)*eye(p+m), ...
                    eye(p+m), []);
    zeta       = lsim(filter, [y; u], t, zeros(p+m, 1))';
    zeta(:, 1) = zeros(p+m, 1);
    Z_new      = interp1(t, zeta', samples, "nearest")';

    % adding the new row to X
    X_hat      = [X; X_new];
    % adding the new row to Z and reordering the rows
    Z_hat      = kron(eye(p+m), [     eye(nu_hat-1); ...
                                 zeros(1, nu_hat-1)]) * Z;
    Z_new      = kron(eye(p+m), [zeros(nu_hat-1, 1); ...
                                                  1]) * Z_new;
    Z_hat      = Z_hat + Z_new;

    % dataset
    XZ     = [X_hat; Z_hat];

    disp('nu_hat:')
    disp(nu_hat)
    disp('Rank of [X; Z]:')
    disp(rank(XZ, tol))
    disp('Number of rows of [X; Z]')
    disp(height(XZ))

    if rank(XZ, tol) < height(XZ)
        % if the data lose rank:
        nu_found = true;
    else
        % otherwise:
        % saving the data for the next iteration
        X = X_hat;
        Z = Z_hat;
    end
end

% observability index
nu = nu_hat - 1;

end