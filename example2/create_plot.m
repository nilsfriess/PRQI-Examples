function create_plot(test, size, weight_target_max, target_index, n_tests, filename)

[A,n] = test_matrix(test, size);
[V,D] = eigs(A, n, "smallestreal");
D     = diag(D);

add_weights = linspace(1, weight_target_max, n_tests);

classic_rqi_res = zeros(n_tests, 2);
prqi_res        = zeros(n_tests, 2);

% Select target eigenvector and eigenvalue
target_v = V(:, target_index);
target_e = D(target_index);

for i = 1 : n_tests
    weights = rand(n, 1);
    weights(target_index) = add_weights(i);

    x0 = V * weights;
    x0 = x0 / norm(x0);

    e_rqi  = classic_rqi(A, x0, 1e-15);
    e_prqi = prqi(A, x0, 1e-15);

    classic_rqi_res(i,1) = acos(abs(x0'*target_v));
    prqi_res(i,1)        = acos(abs(x0'*target_v));

    classic_rqi_res(i,2) = e_rqi;
    prqi_res(i,2)        = e_prqi;
end

rqi_corr_ind  = find(abs(classic_rqi_res(:,2) - target_e) < 1e-14);
rqi_wrong_ind = find(abs(classic_rqi_res(:,2) - target_e) >= 1e-14);

prqi_corr_ind  = find(abs(prqi_res(:,2) - target_e) < 1e-14);
prqi_wrong_ind = find(abs(prqi_res(:,2) - target_e) >= 1e-14);

success_colour = [0.9020, 0.6235, 0];
fail_colour    = [0, 0.4471, 0.6980];

fig = figure();
set(fig, 'Visible', 'off');

hold on;

scatter(classic_rqi_res(rqi_corr_ind,1),...
        classic_rqi_res(rqi_corr_ind,2),...
        30,...
        "marker", '^', ...
        "markeredgecolor", success_colour,...
        "LineWidth", 1.5);

scatter(classic_rqi_res(rqi_wrong_ind,1),...
        classic_rqi_res(rqi_wrong_ind,2),...
        30,...
        "marker", '^', ...
        "markeredgecolor", fail_colour,...
        "LineWidth", 1.5);

scatter(prqi_res(prqi_corr_ind,1),...
        prqi_res(prqi_corr_ind,2),...
        80,...
        "marker", 'o', ...
        "markeredgecolor", success_colour,...
        "LineWidth", 1.5);

scatter(prqi_res(prqi_wrong_ind,1),...
        prqi_res(prqi_wrong_ind,2),...
        80,...
        "marker", 'o', ...
        "markeredgecolor", fail_colour,...
        "LineWidth", 1.5);

plot([min(prqi_res(:,1)) - 1, pi/2 + 1], [target_e, target_e], 'k');

ylabel('Computed eigenvalue', 'interpreter', 'latex')
pbaspect([ 1 1 1 ])
xlabel('Angle between initial vector and target', 'interpreter', 'latex')

legend({'Classic RQI (correct eigenvalue)', ...
        'Classic RQI (wrong eigenvalue)' ...
        'PRQI (correct eigenvalue)', ...
        'PRQI (wrong eigenvalue)' ...
        'Target eigenvalue'...
        },...
        'interpreter','latex', ...
        'location', 'northeast');

xlim([pi/12 - 0.1, pi/2 + 0.1]);
set(gca,'TickLabelInterpreter','latex');

set(gca,...
    "Xtick", [pi/12, pi/6, pi/4, pi/3, 5*pi/12, pi/2]);
set(gca,...
    "XtickLabel", ["$\pi/12$", "$\pi/6$", "$\pi/4$", "$\pi/3$", "$5\pi/12$", "$\pi/2$"]);

curr_ylim = ylim;
ylim([curr_ylim(1) - 0.1, curr_ylim(2)]);

box on;
exportgraphics(fig, filename, 'ContentType', 'vector');
end
