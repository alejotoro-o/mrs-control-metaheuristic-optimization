clc
clear; close all;

for T = 1:1
    tic;

    %% Multi-Verse Optimizer (MVO) – Generic/Commented Template
    % Ni   : population size (universes)
    % Nv   : number of decision variables (dimension)
    % tmax : max iterations
    % iterno : patience window (no-improvement early stop)
    % Bounds: xmin, xmax  (row vectors of length Nv)
    % Objective: replace "funcionCD_UR3_Jac_orientacion" by your function if needed

    Ni    = 50;        % population size
    Nv    = 6;         % number of variables
    tmax  = 1000;      % max iterations
    iterno = 1000;     % no-improvement iterations (early stop)
    tx    = 0;         % no-improvement counter

    xmin = deg2rad([-180 -180 -180 -180 -180 -180]);
    xmax = deg2rad([ 180  180  180  180  180  180]);

    % --- MVO hyperparameters (common choices) ---
    % WEP grows from small to large (exploration -> exploitation)
    WEP_min = 0.2;     % you had min > max ; set a standard schedule
    WEP_max = 1.0;
    p       = 3;       % parameter for TDR (higher p -> slower decrease)

    % --- Initial population (universes), with last column for fitness ---
    X  = xmin + (xmax - xmin).*rand(Ni, Nv);
    fo = zeros(Ni,1);

    % Evaluate initial population
    parfor rr = 1:Ni
        fo(rr,1) = funcionCD_UR3_Jac_orientacion(X(rr,:));
    end
    X = [X, fo];  % append fitness

    % Best universe (elite)
    [~, gbest] = min(X(:, end));
    U_best     = X(gbest, :);   % [1 x (Nv+1)]
    X_best1    = U_best;        % best-so-far for early stopping

    % --- Main loop ---
    for t = 1:tmax
        % Wormhole existence probability (WEP) increases with t
        WEP = WEP_min + (WEP_max - WEP_min) * (t / tmax);
        % Travelling distance rate (TDR) decreases with t
        TDR = 1 - ( (t)^(1/p) / (tmax)^(1/p) );

        % Sort population by fitness (ascending)
        SU = sortrows(X, Nv+1);

        % --- Build selection probabilities (roulette wheel) ---
        fvals = SU(:, end);
        fmin  = fvals(1);
        fmax  = fvals(end);
        if fmax > fmin
            w = (fmax - fvals) ./ (fmax - fmin + eps); % best->1, worst->0
        else
            w = ones(Ni,1) / Ni;                       % all equal fitness
        end
        prob = w / sum(w);
        cdf  = cumsum(prob);                           % valid CDF in (0,1]

        % Keep the best in position 1 (elitism)
        X_new = SU;            % copy (will overwrite decision vars)
        X_new(1, :) = SU(1, :);

        % --- Update universes 2..Ni ---
        % Two mechanisms:
        %  (1) White/black hole (roulette-based) replacement
        %  (2) Wormhole (exploit around the elite using WEP and TDR)
        for i = 2:Ni
            for j = 1:Nv
                % 1) White/black hole interaction
                if rand < prob(i)   % probability biased by universe quality
                    r = rand;
                    idx = find(cdf >= r, 1, 'first');
                    X_new(i, j) = SU(idx, j);
                end

                % 2) Wormhole flying around the elite (with probability WEP)
                if rand < WEP
                    if rand < 0.5
                        X_new(i, j) = U_best(1, j) + TDR * ((xmax(j) - xmin(j))*rand + xmin(j));
                    else
                        X_new(i, j) = U_best(1, j) - TDR * ((xmax(j) - xmin(j))*rand + xmin(j));
                    end
                end
            end
        end

        % --- Enforce variable bounds (saturation) ---
        X_new(:, 1:Nv) = min(max(X_new(:, 1:Nv), xmin), xmax);

        % --- Re-evaluate fitness ---
        fo = zeros(Ni,1);
        parfor rr = 1:Ni
            fo(rr,1) = funcionCD_UR3_Jac_orientacion(X_new(rr,1:Nv));
        end
        X_new(:, Nv+1) = fo;

        % --- Update global best ---
        [curr_best_val, idxb] = min(X_new(:, end));
        if curr_best_val < U_best(1, end)
            U_best(1, :) = X_new(idxb, :);
        end

        % --- Early-stopping by no-improvement ---
        if U_best(1, end) < X_best1(1, end)
            tx      = 0;
            X_best1 = U_best;
        else
            tx = tx + 1;
        end
        if tx >= iterno
            break;
        end

        % Prepare for next iteration
        X = X_new;
        % Optional: fprintf('Iter %4d, Best %.6e\n', t, X_best1(1,end));
    end

    % --- Bookkeeping like your other scripts ---
    time = toc;
    Mejor_Sol(T, :) = X_best1(1, 1:Nv);  % best decision vector
    Ploss_min(T, 1) = X_best1(1, end);   % best objective
    Times(T, 1)     = time;
end

% --- Summary stats (same fields you use elsewhere) ---
MediaTime          = mean(Times);
DesviacionestTime  = std(Times);
PorcentDesTime     = (DesviacionestTime/MediaTime)*100;

Mediamejorsol      = mean(Ploss_min);
DesviacionestPloss = std(Ploss_min);
PorcentDesSol      = (DesviacionestPloss/Mediamejorsol)*100;

bestsolution1 = min(Ploss_min(:,1));
[x, ~]        = find(Ploss_min(:,1) == bestsolution1, 1, 'first');
bestsolution  = Mejor_Sol(x,:);

Solution = [bestsolution bestsolution1 Mediamejorsol PorcentDesSol MediaTime];

save results/MVO_Experimental
