function [e, varargout] = prqi(A, x0, tol)
% PRQI  Performs the Projected Rayleigh Quotient
%       Iteration (RQI) for the eigenvalue problem
%       Av = \lambda v
x0 = x0(:);
v = x0 / norm(x0);

I = speye(size(A));

e = v'*A*v;              % Rayleigh quotient
res = norm((A - e*I)*v); % Residual

its = 0;
while res > tol
    v = (A - (e - 1i*res)*I) \ v;
    v = v / norm(v);

    e = v'*A*v;
    res = norm((A - e*I)*v);

    its = its + 1;

    if its > 15
        break;
    end
end

e = real(e);

varargout{1} = v;
varargout{2} = its;
end
