%% ANALYZE ALGORITHMS: BARS, RADAR, CONVERGENCE, AND STATS
% Metrics: Best Objective (min), Mean Objective, Std Objective, Mean Time (s), Std Time (s)
% Plus: Convergence curves (median + IQR), nonparametric stats with FDR-adjusted p-values,
%       and Cliff's delta effect sizes.
clear; clc;

%% ------------------ File discovery ------------------
prettyName = @(s) strrep(strrep(s,'_Experimental.mat',''),'_','-');
files = dir('*_Experimental.mat');
if isempty(files)
    error('No *_Experimental.mat files found in the current folder.');
end

%% ------------------ Load & extract metrics ------------------
K = numel(files);
names   = strings(K,1);
bestFO  = nan(K,1);   % Best objective (min across runs)
meanFO  = nan(K,1);   % Mean objective across runs
stdFO   = nan(K,1);   % Std objective across runs
meanT   = nan(K,1);   % Mean time (s)
stdT    = nan(K,1);   % Std time (s)

% Raw arrays (for statistical tests)
Ploss_all  = cell(K,1);   % each cell: vector of best per run (Ploss_min)
Times_all  = cell(K,1);   % each cell: vector of time per run
Conv_all   = cell(K,1);   % each cell: convergence matrix [iters x runs]
Conv_names = strings(0,1);

for k = 1:K
    S = load(files(k).name);
    names(k) = prettyName(files(k).name);

    hasSolution = isfield(S,'Solution')   && ~isempty(S.Solution);
    hasPloss    = isfield(S,'Ploss_min')  && ~isempty(S.Ploss_min);
    hasTimes    = isfield(S,'Times')      && ~isempty(S.Times);

    % --- Final metrics from Solution / raw arrays ---
    if hasSolution
        Nv = 6; vec = S.Solution(:).';
        if numel(vec) >= Nv+4
            bestFO(k) = vec(Nv+1);
            meanFO(k) = vec(Nv+2);
            if hasPloss
                stdFO(k) = std(S.Ploss_min(:));
                Ploss_all{k} = S.Ploss_min(:);
            else
                if numel(vec) >= Nv+3 && ~isnan(vec(Nv+3)) && vec(Nv+3) ~= 0
                    stdFO(k) = (vec(Nv+3)/100) * meanFO(k);
                else
                    stdFO(k) = NaN;
                end
            end
            if hasTimes
                meanT(k) = mean(S.Times(:));
                stdT(k)  = std(S.Times(:));
                Times_all{k} = S.Times(:);
            else
                meanT(k) = vec(Nv+4);
                stdT(k)  = NaN;
            end
        else
            if hasPloss
                bestFO(k) = min(S.Ploss_min(:));
                meanFO(k) = mean(S.Ploss_min(:));
                stdFO(k)  = std(S.Ploss_min(:));
                Ploss_all{k} = S.Ploss_min(:);
            end
            if hasTimes
                meanT(k) = mean(S.Times(:));
                stdT(k)  = std(S.Times(:));
                Times_all{k} = S.Times(:);
            end
        end
    else
        if hasPloss
            bestFO(k) = min(S.Ploss_min(:));
            meanFO(k) = mean(S.Ploss_min(:));
            stdFO(k)  = std(S.Ploss_min(:));
            Ploss_all{k} = S.Ploss_min(:);
        else
            error('File %s lacks both Solution and Ploss_min.', files(k).name);
        end
        if hasTimes
            meanT(k) = mean(S.Times(:));
            stdT(k)  = std(S.Times(:));
            Times_all{k} = S.Times(:);
        end
    end

    % --- Convergence extraction (optional) ---
    C = extract_convergence_struct(S);
    if ~isempty(C)
        Conv_all{k} = C;  % [iters x runs]
        Conv_names(end+1,1) = names(k); %#ok<AGROW>
    end
end

%% ------------------ Bar charts ------------------
figure('Name','Algorithm comparison: bar charts','Color','w','Position',[100 100 1300 720]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; 
bar(bestFO); 
title('Best Objective (min)'); ylabel('Objective , lower is better');
set(gca,'XTick',1:K,'XTickLabel',names,'XTickLabelRotation',30); grid on;

nexttile; 
bar(meanFO); 
title('Mean Objective'); ylabel('Objective , lower is better');
set(gca,'XTick',1:K,'XTickLabel',names,'XTickLabelRotation',30); grid on;

nexttile; 
bar(stdFO); 
title('Std Objective'); ylabel('Std , lower is better');
set(gca,'XTick',1:K,'XTickLabel',names,'XTickLabelRotation',30); grid on;

nexttile; 
bar(meanT); 
title('Mean Time (s)'); ylabel('Time (s) , lower is better');
set(gca,'XTick',1:K,'XTickLabel',names,'XTickLabelRotation',30); grid on;

nexttile; 
bar(stdT); 
title('Std Time (s)'); ylabel('Std (s) , lower is better');
set(gca,'XTick',1:K,'XTickLabel',names,'XTickLabelRotation',30); grid on;

nexttile; axis off;
text(0,0.7,'Metrics shown in English:','FontWeight','bold');
text(0,0.5,'- Best Objective (min)','Interpreter','none');
text(0,0.4,'- Mean Objective','Interpreter','none');
text(0,0.3,'- Std Objective','Interpreter','none');
text(0,0.2,'- Mean Time (s)','Interpreter','none');
text(0,0.1,'- Std Time (s)','Interpreter','none');

%% ------------------ Radar chart ------------------
M = [bestFO, meanFO, stdFO, meanT, stdT];
radar_labels = {'Best Objective (min)','Mean Objective','Std Objective','Mean Time (s)','Std Time (s)'};

% Handle NaNs , replace by worst observed to keep fair normalization
for j = 1:size(M,2)
    col = M(:,j);
    if all(isnan(col))
        col = zeros(size(col));
    else
        col(isnan(col)) = max(col(~isnan(col)));
    end
    M(:,j) = col;
end

% Min–max normalize and invert (all metrics are “lower is better”)
Mn = zeros(size(M));
for j = 1:size(M,2)
    cmin = min(M(:,j)); cmax = max(M(:,j));
    if cmax > cmin
        Mn(:,j) = 1 - (M(:,j) - cmin) / (cmax - cmin); % invert so higher is better
    else
        Mn(:,j) = ones(K,1); % equal -> 1 after inversion
    end
end

figure('Name','Algorithm comparison: radar','Color','w','Position',[100 100 780 780]);
spider_plot_basic(Mn, radar_labels, names);

%% ------------------ Convergence plots ------------------
% Combined convergence (median + IQR) if any algorithm provides per-iteration data
if ~isempty(Conv_names)
    figure('Name','Convergence , median with IQR band','Color','w','Position',[100 100 1200 700]);
    hold on; grid on; box on;
    cmap = lines(numel(Conv_names));
    leg  = strings(0,1);
    for i = 1:K
        C = Conv_all{i};
        if isempty(C), continue; end
        % C expected shape: [iters x runs]
        [q25, med, q75] = quantiles_over_runs(C);
        x = (1:numel(med))';
        shaded_iqr(x, q25, q75, cmap(i,:));
        plot(x, med, 'LineWidth', 1.8, 'Color', cmap(i,:));
        leg(end+1,1) = names(i); %#ok<AGROW>
    end
    xlabel('Iteration'); ylabel('Best-so-far objective'); 
    title('Convergence , median with interquartile band (lower is better)');
    if ~isempty(leg), legend(leg,'Location','northeastoutside'); end
    hold off;

    % Small multiples , one tile per algorithm
    figure('Name','Convergence per algorithm','Color','w','Position',[100 100 1300 720]);
    valid_idx = find(~cellfun(@isempty, Conv_all));
    nV = numel(valid_idx); 
    nrows = ceil(nV/3); ncols = min(3,nV);
    tiledlayout(nrows,ncols,'Padding','compact','TileSpacing','compact');
    for ii = 1:nV
        i = valid_idx(ii);
        nexttile; grid on; box on; hold on;
        C = Conv_all{i};
        [q25, med, q75] = quantiles_over_runs(C);
        x = (1:numel(med))';
        shaded_iqr(x, q25, q75, [0.3 0.5 0.9]);
        plot(x, med, 'LineWidth', 1.8, 'Color', [0 0.2 0.6]);
        xlabel('Iteration'); ylabel('Best-so-far objective');
        title(names(i), 'Interpreter','none');
        hold off;
    end
end

%% ------------------ Statistical analysis ------------------
% Build descriptive summary table
Summary = table(names, bestFO, meanFO, stdFO, meanT, stdT, ...
    'VariableNames', {'Algorithm','BestObjective','MeanObjective','StdObjective','MeanTime_s','StdTime_s'});

% Kruskal–Wallis on final best objective and times
% (More robust than ANOVA when distributions are non-normal)
[KW_FO_p, ~, KW_FO_stats] = kruskalwallis_grouped(Ploss_all, names);
[KW_T_p,  ~, KW_T_stats ] = kruskalwallis_grouped(Times_all, names);

% Pairwise ranksum with FDR adjustment and Cliff's delta
[PairsFO, P_FO_adj, DeltaFO] = pairwise_ranksum_cliff(Ploss_all, names);
[PairsT,  P_T_adj,  DeltaT ] = pairwise_ranksum_cliff(Times_all, names);

% Display quick console summary
disp('----- Descriptive Summary -----'); disp(Summary);
fprintf('\nKruskal–Wallis (Final Objective): p = %.3g (chi2 = %.3f, df = %d)\n', KW_FO_p, KW_FO_stats.chi2, KW_FO_stats.df);
fprintf('Kruskal–Wallis (Time):            p = %.3g (chi2 = %.3f, df = %d)\n', KW_T_p,  KW_T_stats.chi2,  KW_T_stats.df);

% Export CSVs
writetable(Summary, 'summary_algorithms.csv');
writetable(PairsFO, 'pairwise_objective_pvalues_FDR.csv');
writetable(DeltaFO, 'pairwise_objective_cliffs_delta.csv');
writetable(PairsT,  'pairwise_time_pvalues_FDR.csv');
writetable(DeltaT,  'pairwise_time_cliffs_delta.csv');

fprintf('\nSaved CSVs:\n  - summary_algorithms.csv\n  - pairwise_objective_pvalues_FDR.csv\n  - pairwise_objective_cliffs_delta.csv\n  - pairwise_time_pvalues_FDR.csv\n  - pairwise_time_cliffs_delta.csv\n');

%% ------------------ Radar function (Cartesian, supports patch) ------------------
function spider_plot_basic(data, labels, series)
    if isstring(labels), labels = cellstr(labels); end
    if isstring(series), labels = cellstr(labels); end %#ok<NASGU>
    if isstring(series), series = cellstr(series); end

    [N, D] = size(data);
    if length(labels) ~= D
        error('Labels length must match data columns.');
    end

    theta = linspace(0, 2*pi, D + 1);
    theta(end) = theta(1);

    ax = axes('NextPlot','add'); %#ok<LAXES>
    axis(ax, 'equal'); axis(ax, 1.25*[-1 1 -1 1]); ax.Visible = 'off'; hold(ax, 'on');
    title('Normalized performance radar , higher is better');

    nRings = 5;
    rticks = linspace(0,1,nRings+1);
    thFine = linspace(0,2*pi,360);
    for r = rticks
        [xg, yg] = pol2cart(thFine, r);
        plot(ax, xg, yg, ':', 'LineWidth', 0.8);
    end

    for d = 1:D
        plot(ax, [0 cos(theta(d))], [0 sin(theta(d))], ':', 'LineWidth', 0.8);
        [xl, yl] = pol2cart(theta(d), 1.12);
        text(ax, xl, yl, labels{d}, 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle','FontWeight','bold');
    end

    cmap = lines(max(N,7));
    for i = 1:N
        r = [data(i,:), data(i,1)];
        [x, y] = pol2cart(theta, r);
        patch('XData', x, 'YData', y, 'FaceColor', cmap(i,:), ...
              'FaceAlpha', 0.12, 'EdgeColor', 'none', 'Parent', ax);
        plot(ax, x, y, 'LineWidth', 1.5, 'Color', cmap(i,:));
        plot(ax, x, y, 'o', 'MarkerSize', 5, 'Color', cmap(i,:));
    end

    for k = 2:length(rticks)
        txt = sprintf('%.1f', rticks(k));
        text(ax, rticks(k)*cos(pi/2), rticks(k)*sin(pi/2), txt, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
    end

    legend(series, 'Location','southoutside','NumColumns', max(2,ceil(N/2)));
    hold(ax, 'off');
end

%% ------------------ Convergence extraction helper ------------------
% Tries various common field names; returns [] if not found.
% Expected: either a vector per run (best-so-far per iteration) or a matrix.
% Returns a matrix [iters x runs].
function C = extract_convergence_struct(S)
    C = [];
    % Common field candidates
    cands = {'report_conver','ReportConver','convergence','Convergence', ...
             'BestHistory','History','BestPerIter','Best_curve','curve'};
    fn = fieldnames(S);
    hit = '';
    for i = 1:numel(cands)
        if any(strcmpi(fn, cands{i}))
            hit = cands{i}; break;
        end
    end
    if isempty(hit), return; end

    raw = S.(hit);
    if isempty(raw), return; end

    % Accept cell array of runs or numeric matrix
    if iscell(raw)
        % Each cell can be a row/column vector; pad to the same length
        L = cellfun(@numel, raw);
        T = max(L);
        R = numel(raw);
        C = nan(T, R);
        for r = 1:R
            v = raw{r}(:);
            C(1:numel(v), r) = v;
        end
    elseif isnumeric(raw)
        A = raw;
        if isvector(A)
            C = A(:);  % [iters x 1]
        else
            % Heuristic: if one axis equals number of runs saved elsewhere, keep as is
            % We prefer [iters x runs]; if detected transposed, transpose.
            if size(A,1) >= size(A,2)
                C = A;
            else
                C = A.'; % transpose to [iters x runs]
            end
        end
    else
        C = [];
    end
end

%% ------------------ Quantiles helper ------------------
% Input: C [iters x runs]; Output: q25, med, q75 (column vectors)
function [q25, med, q75] = quantiles_over_runs(C)
    if isvector(C), C = C(:); end
    q25 = nan(size(C,1),1);
    med = q25; q75 = q25;
    for t = 1:size(C,1)
        v = C(t, :); v = v(~isnan(v));
        if isempty(v), q25(t)=NaN; med(t)=NaN; q75(t)=NaN; continue; end
        q = quantile(v, [0.25 0.5 0.75]);
        q25(t) = q(1); med(t) = q(2); q75(t) = q(3);
    end
end

%% ------------------ Shaded IQR helper ------------------
function shaded_iqr(x, y1, y2, col)
    X = [x(:); flipud(x(:))];
    Y = [y1(:); flipud(y2(:))];
    patch('XData', X, 'YData', Y, 'FaceColor', col, 'EdgeColor','none', 'FaceAlpha', 0.12);
end

%% ------------------ Kruskal–Wallis wrapper ------------------
% Input: cell array of vectors and labels. Missing/empty groups are ignored.
function [p, tbl, stat] = kruskalwallis_grouped(groups, labels)
    vals = []; grp = [];
    for i = 1:numel(groups)
        g = groups{i};
        if isempty(g), continue; end
        vals = [vals; g(:)];
        grp  = [grp; repmat(i, numel(g), 1)];
    end
    if isempty(vals)
        p = NaN; tbl = []; stat = struct('chi2',NaN,'df',NaN);
        return
    end
    [p, ~, stats] = kruskalwallis(vals, grp, 'off');
    % Recover chi2 and df from stats (via ANOVA table if needed)
    % Quick approximation:
    try
        t = stats; %#ok<NASGU>
        % MATLAB doesn't return chi2 directly; we leave NaN if not trivial to recover
        stat = struct('chi2', NaN, 'df', numel(unique(grp))-1);
    catch
        stat = struct('chi2', NaN, 'df', numel(unique(grp))-1);
    end
    tbl = table(labels, 'VariableNames', {'Algorithms'});
end

%% ------------------ Pairwise ranksum + FDR + Cliff's delta ------------------
function [Ptable, P_adj_table, DeltaTable] = pairwise_ranksum_cliff(groups, names)
    K = numel(groups);
    pairs = {}; pvals = []; A = []; B = []; D = [];
    for i = 1:K
        xi = groups{i}; xi = xi(:);
        if isempty(xi), continue; end
        for j = i+1:K
            xj = groups{j}; xj = xj(:);
            if isempty(xj), continue; end
            % Ranksum p-value
            p = ranksum(xi, xj);
            % Cliff's delta
            d = cliffs_delta(xi, xj);
            pairs{end+1,1} = sprintf('%s vs %s', names(i), names(j)); %#ok<AGROW>
            pvals(end+1,1) = p; %#ok<AGROW>
            A(end+1,1) = i; B(end+1,1) = j; %#ok<AGROW>
            D(end+1,1) = d; %#ok<AGROW>
        end
    end
    % FDR adjustment
    p_adj = fdr_bh(pvals);
    Ptable = table(pairs, pvals, 'VariableNames', {'Pair','p_value'});
    P_adj_table = table(pairs, p_adj, 'VariableNames', {'Pair','p_value_FDR'});
    DeltaTable  = table(pairs, D,     'VariableNames', {'Pair','CliffsDelta'});
end

%% ------------------ Benjamini–Hochberg FDR ------------------
function p_adj = fdr_bh(p)
    p = p(:); n = numel(p);
    [sp, ix] = sort(p);
    adj = sp .* n ./ (1:n)';
    % ensure monotonicity
    for i = n-1:-1:1
        adj(i) = min(adj(i), adj(i+1));
    end
    p_adj = zeros(n,1);
    p_adj(ix) = min(adj, 1);
end

%% ------------------ Cliff's delta ------------------
% d in [-1,1]; 0 ~ no effect, +/-1 maximal effect
function d = cliffs_delta(x, y)
    x = x(:); y = y(:);
    nx = numel(x); ny = numel(y);
    % Efficient pairwise compare via sorting ranks
    % Simpler O(n*m) version (adequate for moderate sizes):
    greater = 0; less = 0;
    for i = 1:nx
        greater = greater + sum(x(i) > y);
        less    = less    + sum(x(i) < y);
    end
    d = (greater - less) / (nx*ny);
end
