function [A,B,M,x,varargout] =...
    assemble_schroedinger(X, N, comp_spectrum, comp_pert_spectrum)

if nargin < 3
    comp_spectrum = false;
end

if nargin < 4
    comp_pert_spectrum = false;
end

x = linspace(0, X, N + 2);
x = x(:);

G = @(xs) sin(xs) - 40./(1+xs.^2);
F = @(x) ones(size(x));

BC = [0,1,0; 1,0,0]; % boundary conditions as expected by sturm.m

[A, B, M] = sturm(x, BC, F, G, F);


opts = struct;
opts.p = 600;

if comp_spectrum
    E = eigs(A+B, M, 250, -50, opts);
end

if comp_pert_spectrum
    %% Construct initial vector as described in the paper and
    %% compute spectrum of A + i*gamma*(I - x*x')
    n_osc = 6;
    R = 35;

    x0 = initial_vec(x, n_osc, R);
    x0 = x0 / norm(x0);

    Ep = eigs(A+B + 1i*M*(speye(size(A,1)) - x0*x0'), M, 250, -50, opts);
end

if comp_spectrum
    varargout{1} = E;
end

if comp_pert_spectrum
    if comp_spectrum
        varargout{2} = Ep;
    else
        varargout{1} = Ep;
    end

end

end
