clc
clear all

for T = 1:1
    tic;

    %% JADE (Adaptive Differential Evolution) with p-best and archive
    Ni     = 50;      % population size
    Nv     = 6;       % number of variables
    tmax   = 1000;    % max iterations
    iterno = 10;      % no-improvement window
    tx     = 0;

    xmin = deg2rad([-180 -180 -180 -180 -180 -180]);
    xmax = deg2rad([ 180  180  180  180  180  180]);

    % JADE hyperparameters
    p      = 0.1;     % percentile for p-best
    c_adapt = 0.1;    % learning rate for mu_F and mu_CR
    mu_F   = 0.5;     % initial mean of F (Cauchy)
    mu_CR  = 0.5;     % initial mean of CR (Gaussian)

    % Initialize population
    X = xmin.*ones(Ni, Nv) + (xmax - xmin).*rand(Ni, Nv);
    f = zeros(Ni,1);
    parfor i = 1:Ni
        f(i) = funcionCD_UR3_Jac_orientacion(X(i,:));
    end

    % Archive for JADE
    A = zeros(0, Nv);

    % Track best
    [fbest, ibest] = min(f);
    X_best  = [X(ibest,:), fbest];
    X_best1 = X_best;

    for t = 0:tmax
        if t > 0
            % Determine p-best set
            pnum = max(2, round(p * Ni)); % at least 2
            [~, order] = sort(f, 'ascend');
            pbest_idx = order(1:pnum);

            % Successful parameter memories
            SF = [];  % successful F
            SCR = []; % successful CR
            dF = [];  % for Lehmer mean numerator

            Xnew = X; fnew = f;
            for i = 1:Ni
                % --- Parameter adaptation ---
                % F ~ Cauchy(mu_F, 0.1), truncated to (0,1]
                Fi = -1;
                while Fi <= 0
                    Fi = mu_F + 0.1 * tan(pi*(rand - 0.5));
                end
                Fi = min(Fi, 1);

                % CR ~ N(mu_CR, 0.1), clipped to [0,1]
                CRi = min(max(mu_CR + 0.1*randn, 0), 1);

                % --- Mutation: "current-to-pbest/1" ---
                % pick p-best
                kp = pbest_idx(randi(pnum));
                % pick r1 != i
                r1 = randi(Ni);
                while r1 == i
                    r1 = randi(Ni);
                end
                % pick r2 from union of pop and archive
                if ~isempty(A)
                    % choose from [X; A]
                    XA = [X; A];
                    r2 = randi(size(XA,1));
                    V2 = XA(r2,:);
                else
                    % fallback to population
                    r2i = randi(Ni);
                    while r2i == i || r2i == r1
                        r2i = randi(Ni);
                    end
                    V2 = X(r2i,:);
                end

                Vi = X(i,:) + Fi*(X(kp,:) - X(i,:)) + Fi*(X(r1,:) - V2);

                % --- Crossover: binomial ---
                Ui = X(i,:);
                jrand = randi(Nv);
                for j = 1:Nv
                    if rand <= CRi || j == jrand
                        Ui(j) = Vi(j);
                    end
                end

                % --- Bounds handling: projection (simple and robust) ---
                Ui = min(max(Ui, xmin), xmax);

                % --- Selection ---
                fi = funcionCD_UR3_Jac_orientacion(Ui);
                if fi <= f(i)
                    % Success
                    Xnew(i,:) = Ui;
                    fnew(i)   = fi;

                    % Append the replaced parent to archive
                    A = [A; X(i,:)];

                    % Record success params
                    SF(end+1)  = Fi; %#ok<*AGROW>
                    SCR(end+1) = CRi;
                    dF(end+1)  = Fi^2;
                end
            end

            % Limit archive size
            if size(A,1) > Ni
                idx_remove = randperm(size(A,1), size(A,1) - Ni);
                A(idx_remove,:) = [];
            end

            % Update population
            X = Xnew;
            f = fnew;

            % Update parameter means
            if ~isempty(SF)
                % Lehmer mean for F
                mu_F  = (1 - c_adapt)*mu_F  + c_adapt * (sum(dF) / sum(SF));
                % Arithmetic mean for CR
                mu_CR = (1 - c_adapt)*mu_CR + c_adapt * mean(SCR);
            end

            % Global best update
            [fbest, ibest] = min(f);
            if fbest < X_best(1,end)
                X_best = [X(ibest,:), fbest];
            end

            % No-improvement counter
            if X_best(1,end) < X_best1(1,end)
                tx = 0;
                X_best1 = X_best;
            else
                tx = tx + 1;
            end
            if tx == iterno
                tx = 0;
                break
            end
        end
        % Optional:
        % fprintf('Iter %4d  f=%.6e  muF=%.3f muCR=%.3f\n', t, X_best1(1,end), mu_F, mu_CR);
    end

    time = toc;
    Mejor_Sol(T,:) = X_best1(1,1:Nv);
    Ploss_min(T,1) = X_best1(1,end);
    Times(T,1)     = time;
end

% Summary statistics
MediaTime          = mean(Times);
DesviacionestTime  = std(Times);
PorcentDesTime     = (DesviacionestTime/MediaTime)*100;

Mediamejorsol      = mean(Ploss_min);
DesviacionestPloss = std(Ploss_min);
PorcentDesSol      = (DesviacionestPloss/Mediamejorsol)*100;

bestsolution1 = min(Ploss_min(:,1));
[x, ~] = find(Ploss_min(:,1) == bestsolution1, 1, 'first');
bestsolution = Mejor_Sol(x,:);

Solution = [bestsolution bestsolution1 Mediamejorsol PorcentDesSol MediaTime];

save results/JADE_Experimental
