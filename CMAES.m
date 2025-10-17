clc
clear all

for T = 1:1
    tic;

    %% CMA-ES settings (adapted to your scaffold)
    Nv     = 6;           % number of variables
    tmax   = 1000;        % max iterations
    iterno = 10;          % no-improvement window
    tx     = 0;

    xmin = deg2rad([-180 -180 -180 -180 -180 -180]);
    xmax = deg2rad([ 180  180  180  180  180  180]);

    % Population size consistent with your "Ni" role
    Ni = 10;                      % λ (population size). You can use 4+floor(3*log(Nv)).
    mu = floor(Ni/2);             % number of parents
    % Recombination weights
    w  = log(mu + 0.5) - log(1:mu);
    w  = w / sum(w);
    mu_eff = 1/sum(w.^2);

    % Strategy parameter settings: adaptation
    cc = (4 + mu_eff/Nv) / (Nv + 4 + 2*mu_eff/Nv);     % time constant for cumulation for C
    cs = (mu_eff + 2) / (Nv + mu_eff + 5);             % time constant for cumulation for sigma control
    c1 = 2 / ((Nv + 1.3)^2 + mu_eff);                  % learning rate for rank-one update
    cmu = min(1 - c1, 2 * (mu_eff - 2 + 1/mu_eff) / ((Nv + 2)^2 + mu_eff)); % rank-mu update
    damps = 1 + 2*max(0, sqrt((mu_eff - 1)/(Nv + 1)) - 1) + cs;  % damping for sigma

    % Initial mean and step-size
    m  = (xmin + xmax)/2;           % mean at box mid-point
    sigma = 0.3;                    % initial global step size (rad)

    % Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(1,Nv);
    ps = zeros(1,Nv);
    C  = eye(Nv);
    B  = eye(Nv);
    D  = ones(Nv,1);
    BD = B .* D';                   % for sampling
    inv_sqrtC = eye(Nv);

    % Expectation of |N(0,I)| along a random direction
    chiN = sqrt(Nv) * (1 - 1/(4*Nv) + 1/(21*Nv^2));

    % Track best
    X_best  = [m, inf];
    X_best1 = X_best;

    for t = 0:tmax
        if t == 0
            % Evaluate initial mean (optional, not required by CMA-ES, but aligns with your pattern)
            f_m = funcionCD_UR3_Jac_orientacion(m);
            X_best = [m, f_m];
            X_best1 = X_best;
        else
            % Sample λ offspring
            arz = randn(Ni, Nv);                 % N(0, I)
            ary = arz * (B*diag(D))';            % transform to N(0, C)
            Xcand = m + sigma * ary;             % offspring in R^Nv

            % Resample out-of-bounds (simple resampling for feasibility)
            for i = 1:Ni
                tries = 0;
                while any(Xcand(i,:) < xmin) || any(Xcand(i,:) > xmax)
                    % resample that candidate
                    arz(i,:) = randn(1, Nv);
                    ary(i,:) = arz(i,:) * (B*diag(D))';
                    Xcand(i,:) = m + sigma * ary(i,:);
                    tries = tries + 1;
                    if tries > 20
                        % fallback to projection to keep going
                        Xcand(i,:) = min(max(Xcand(i,:), xmin), xmax);
                        break
                    end
                end
            end

            % Evaluate
            fvals = zeros(Ni,1);
            parfor i = 1:Ni
                fvals(i) = funcionCD_UR3_Jac_orientacion(Xcand(i,:));
            end

            % Sort by fitness
            [fvals_sorted, idx] = sort(fvals, 'ascend');
            Xcand_sorted = Xcand(idx,:);
            arz_sorted   = arz(idx,:);   % needed for update of ps

            % Update best
            if fvals_sorted(1) < X_best(1,end)
                X_best = [Xcand_sorted(1,:), fvals_sorted(1)];
            end

            % Recombination: new mean
            m_old = m;
            m = w * Xcand_sorted(1:mu,:);

            % Update evolution paths
            % zmean in isotropic coordinates
            zmean = w * arz_sorted(1:mu,:);
            ps = (1 - cs)*ps + sqrt(cs*(2 - cs)*mu_eff) * (zmean * B);
            hsig = norm(ps) / sqrt(1 - (1 - cs)^(2*(t+1))) / chiN < (1.4 + 2/(Nv + 1));
            pc = (1 - cc)*pc + hsig * sqrt(cc*(2 - cc)*mu_eff) * ((m - m_old) / sigma);

            % Adapt covariance
            artmp = (Xcand_sorted(1:mu,:) - m_old) / sigma;
            C = (1 - c1 - cmu)*C + ...
                c1*(pc'*pc + (1 - hsig)*cc*(2 - cc)*C) + ...
                cmu * (artmp' * diag(w) * artmp);

            % Adapt step size
            sigma = sigma * exp((cs/damps) * (norm(ps)/chiN - 1));

            % Update B and D from C every few iterations or each time for simplicity
            [B, Ddiag] = eig((C + C')/2);
            D = sqrt(max(diag(Ddiag), 1e-30));
            inv_sqrtC = B * diag(1./D) * B';

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
        % Optional progress print
        % fprintf('Iter %4d  f=%.6e  sigma=%.3g\n', t, X_best1(1,end), sigma);
    end

    time = toc;
    Mejor_Sol(T,:) = X_best1(1,1:Nv);
    Ploss_min(T,1) = X_best1(1,end);
    Times(T,1)     = time;
end

% Summary
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

save CMAES_Experimental
