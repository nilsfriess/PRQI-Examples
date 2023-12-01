function [e, varargout] = classic_rqi_general(A, M, x0, mu, tol)
% CLASSIC_RQI_GENERAL  Performs the classic Rayleigh Quotient
%                      Iteration (RQI) for the eigenvalue problem
%                      Av = \lambda Mv
%
%   e = CLASSIC_RQI_GENERAL(A, M, x0) applies RQI with the initial
%                                     vector x0 until the residual
%                                     norm is less than or equal 10^-8
%
%   e = CLASSIC_RQI_GENERAL(A, M, x0, mu, tol) applies RQI with the
%                                              initial vector x0  and
%                                              initial shift mu until
%                                              the residual norm is less
%                                              than or equal to tol
%
x0 = x0(:);
v = x0 / sqrt(x0'*M*x0);

if nargin < 4
    e = v'*A*v; % Rayleigh quotient
else
    e = mu;
end

res = norm((A - e*M)*v); % Residual

if nargin < 5
    tol = 1e-8;
end

its = 0;
while res > tol
    v = (A - e*M) \ M*v;
    v = v / sqrt(v'*M*v);

    e = v'*A*v;
    res = norm((A - e*M)*v);

    its = its + 1;
end

varargout{1} = v;
varargout{2} = its;
end
