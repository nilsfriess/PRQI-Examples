function create_table_results(A, M, E, x)

R_max = 80;
norm_tol = 0.4;

fprintf(['\\begin{tabular}{llcccccc}\\toprule\n',...
         '  && \\multicolumn{3}{c}{PRQI} & \\multicolumn{3}{c}{Classic RQI}\n',...
         '  \\\\\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\n',...
         '  $n_{\\text{osc}}$ & $R$ & Eigenvalue & Index & \\# Its. & Eigenvalue & Index & \\# Its. \\\\\\midrule\n']);

for cutoff = [35, 55]
    for k = 3 : 8
        fprintf("  %3.2g & %i &", k / 2, cutoff);
        skipped = 0;

        x0 = initial_vec(x, k, cutoff);
        v = x0 / sqrt(x0'*M*x0);

        %% Perform PRQI
        its = 0;
        mu = v'*A*v;
        res = norm((A - mu*M)*v);
        while res > 1e-8
            its = its + 1;
            gam = res;

            v = (A - (mu - 1i*gam)*M) \ M*v;
            v = v / sqrt(v'*M*v);

            % Check if current vector is not localized. If yes, break and
            % generate new initial vector
            w = v .* (x(:) > R_max);
            if norm(w) / norm(v) > norm_tol
                skipped = 1;
                break;
            end

            mu = v'*A*v;
            res = norm((A - (mu)*M)*v);
        end

        %% Print results  
        if skipped == 0
            % Compute classic RQI result for comparison
            [mu_rqi, ~, cits] = classic_rqi_general(A, M, x0);

            % Compute the index of the eigenvalue for PRQI
            dists = abs(E - mu);
            [distmin, idxmin] = min(dists);
            if distmin > 1e-7
                idx_prqi = -1;
            else
                idx_prqi = idxmin;
            end

            % Compute the index of the eigenvalue for classic RQI
            dists = abs(E - mu_rqi);
            [distmin, idxmin] = min(dists);
            if distmin > 1e-7
                idx_rqi = -1;
            else
                idx_rqi = idxmin;
            end
    
            fprintf(' % .5f & %3d & %2d & %.5f & %3d & %2d \\\\\n',...
                    mu, idx_prqi, its, mu_rqi, idx_rqi, cits);
        else
            fprintf(" -- & -- & -- & -- & -- & --\\\\\n");
        end
    end
end

fprintf('  \\bottomrule \n\\end{tabular}');

end
