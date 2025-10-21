clc
clear all

for T = 1:1
    tic;

    %% Whale Optimization Algorithm (WOA)
    Ni     = 50;     % Number of whales
    Nv     = 6;      % Number of variables
    tmax   = 1000;   % Maximum iterations
    iterno = 10;     % No-improvement stopping window
    tx     = 0;

    xmin = deg2rad([-180 -180 -180 -180 -180 -180]); % Lower bounds
    xmax = deg2rad([ 180  180  180  180  180  180]); % Upper bounds

    rep = 1;
    Guardar = zeros(rep, Nv+1);
    Tiempo  = zeros(rep, 1);

    for kx = 1:rep
        for t = 0:tmax
            if t == 0
                % Initialize population
                X = xmin.*ones(Ni, Nv) + (xmax - xmin).*rand(Ni, Nv);
                % Append objective column
                X = [X, zeros(Ni,1)];

                % Evaluate objective
                fo = zeros(Ni,1);
                parfor rr = 1:Ni
                    RMSE = funcionCD_UR3_Jac_orientacion(X(rr,1:Nv));
                    fo(rr,1) = RMSE;
                end
                X(:,Nv+1) = fo;

                % Best solution trackers
                Xsorted = sortrows(X, Nv+1);
                X_best  = Xsorted(1,:);   % current global best [1 x (Nv+1)]
                X_best1 = X_best;         % best-so-far for early stop

            else
                % Control parameter a decreases from 2 to 0
                a = 2*(1 - t/tmax);
                b = 1; % spiral shape parameter

                Xnew = X(:,1:Nv); % preallocate new positions

                for i = 1:Ni
                    Xi = X(i,1:Nv);

                    r1 = rand(1,Nv);
                    r2 = rand(1,Nv);
                    A  = 2*a.*r1 - a;     % coefficient vector
                    C  = 2.*r2;           % coefficient vector
                    p  = rand;            % probability for phase selection

                    if p < 0.5
                        if all(abs(A) < 1)
                            % Encircling prey around the best
                            D  = abs(C.*X_best(1,1:Nv) - Xi);
                            X1 = X_best(1,1:Nv) - A.*D;
                            Xnew(i,:) = X1;
                        else
                            % Exploration using a random whale
                            r_idx = randi(Ni);
                            Xrand = X(r_idx,1:Nv);
                            D  = abs(C.*Xrand - Xi);
                            X2 = Xrand - A.*D;
                            Xnew(i,:) = X2;
                        end
                    else
                        % Spiral updating around the best
                        Dprime = abs(X_best(1,1:Nv) - Xi);
                        l = -1 + 2*rand(1,Nv); % in [-1,1]
                        X3 = Dprime.*exp(b.*l).*cos(2*pi.*l) + X_best(1,1:Nv);
                        Xnew(i,:) = X3;
                    end
                end

                % Enforce bounds
                X(:,1:Nv) = max(Xnew, xmin);
                X(:,1:Nv) = min(X(:,1:Nv), xmax);

                % Re evaluate objective
                fo = zeros(Ni,1);
                parfor rr = 1:Ni
                    RMSE = funcionCD_UR3_Jac_orientacion(X(rr,1:Nv));
                    fo(rr,1) = RMSE;
                end
                X(:,Nv+1) = fo;

                % Update best
                [curr_best_val, idx] = min(X(:,Nv+1));
                if curr_best_val < X_best(1,end)
                    X_best(1,:) = X(idx,:);
                end

                % No-improvement counter
                if X_best(1,end) < X_best1(1,end)
                    tx = 0;
                    X_best1 = X_best;
                else
                    tx = tx + 1;
                end

                % Early stop by stagnation
                if tx == iterno
                    tx = 0;
                    break
                end
            end

            % Optional display
            % fprintf('Iteracion= %d, Error= %2.6f \n', t, X_best1(1,end));
            % report_conver(t+1,1) = X_best1(1,end);
        end
    end

    time = toc;
    Mejor_Sol(T,:) = X_best1(1,1:Nv); % Best decision vector
    Ploss_min(T,1) = X_best1(1,end);  % Best objective value
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

save results/WOA_Experimental
