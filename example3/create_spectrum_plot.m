function create_spectrum_plot(E_unpert, E_pert, filename, sp_idx, tgt_idx, tgt_idx_pert)
fig = figure();
fig.Units               = 'centimeters';
fig.Position(3)         = 18;
fig.Position(4)         = 6;

orange = [0.9020, 0.6235, 0];
blue   = [0, 0.4471, 0.6980];
red    = [238/255, 75/255, 43/255];

if nargin > 3
    ylim([-0.1,1.1]);
    ax = gca;
    ax.XTick = -1 : 0.1 : 2;

    obj.X = [0.59480, 0.91806, 0.91806, 0.59480];
    obj.Y = [-1, -1, 1.11, 1.11];
    obj.FaceColor	= 'k';
    obj.EdgeColor = 'k';
    obj.FaceAlpha = 0;

    p1 = patch(obj);
    hatchfill2(p1, 'single',...
               'HatchAngle', 75,...
               'HatchDensity', 30,...
               'HatchColor', 'black');
    hold on;

    xlim([0.3, 1]);
    E_unpert2 = E_unpert;
    E_unpert2(sp_idx) = -100;
    E_unpert2(tgt_idx) = -100;

    E_pert2 = E_pert;
    E_pert2(tgt_idx_pert) = -100;

    p2 = plot(real(E_unpert2), imag(E_unpert2), 'o', 'markersize', 7, 'markeredgecolor', orange, 'markerfacecolor', orange);
    p3 = plot(E_unpert(tgt_idx), 0, 'x', 'markerfacecolor', red, 'markeredgecolor', red, 'markersize', 7, 'LineWidth', 1.5);
    p4 = plot(E_unpert(sp_idx), 0, 'o', 'markerfacecolor', orange, 'markeredgecolor', 'black', 'markersize', 7);

    p5 = plot(real(E_pert2), imag(E_pert2), '^', 'markersize', 8, 'markeredgecolor', blue, 'markerfacecolor', blue);
    p6 = plot(real(E_pert(tgt_idx_pert)), imag(E_pert(tgt_idx_pert)), '+', 'markerfacecolor', blue, 'markeredgecolor', blue, 'markersize', 5, 'LineWidth', 1.5);

    plots = [p1,p2,p3,p4,p5,p6];
    legendData = {'Spectral band',...
                  'Spectrum of the original problem',...
                  'Target eigenvalue',...
                  'Spurious eigenvalue',...
                  'Spectrum of the perturbed problem',...
                  'Perturbed target'};
    [legend_h, object_h, plot_h, text_str] = legendflex(plots, legendData,...
        'anchor', [4 4], 'buffer', [-12 8]);
    hatchfill2(object_h(7), 'single',...
               'HatchAngle', 70,...
               'HatchDensity', 30,...
               'HatchColor', 'black');
    object_h(7).FaceColor = 'w';
    for i = 1:6
        set(object_h(i), 'interpreter', 'latex');
    end
else
    fig.Position(4) = 5;

    plot(real(E_unpert), imag(E_unpert), 'o', 'markersize', 6, 'markeredgecolor', orange, 'markerfacecolor', orange);
    hold on;
    plot(real(E_pert), imag(E_pert), '^', 'markersize', 6, 'markeredgecolor', blue, 'markerfacecolor', blue);
    legend("Spectrum of the original problem", "Spectrum of the perturbed problem", "location", "east", 'Interpreter','latex');
    xlim([-0.56, 2]);
    ylim([-0.09,1.09])

    ax = gca;
    ax.XTick = -1 : 0.5 : 2;
end

xlabel("Real",'interpreter','latex');
ylabel("Imaginary",'interpreter','latex');

set(gca,'TickLabelInterpreter','latex')

box on;
exportgraphics(fig, filename, 'ContentType', 'vector');
end

