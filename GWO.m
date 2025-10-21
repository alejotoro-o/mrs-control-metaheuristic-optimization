clc
clear all

for T = 1:1
    tic;

    %% Grey Wolf Optimizer
    Ni    = 50;     % Number of wolves
    Nv    = 6;      % Number of variables
    tmax  = 1000;   % Maximum iterations
    iterno = 10;    % No-improvement stopping window
    tx    = 0;

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
                % Append FO column
                X = [X, zeros(Ni, 1)];

                % Evaluate objective
                fo = zeros(Ni,1);
                parfor rr = 1:Ni
                    RMSE = funcionCD_UR3_Jac_orientacion(X(rr,1:Nv));
                    fo(rr,1) = RMSE;
                end
                X(:, Nv+1) = fo;

                % Identify alpha, beta, delta
                Xsorted = sortrows(X, Nv+1);
                Alpha = Xsorted(1,:);   % best
                Beta  = Xsorted(2,:);
                Delta = Xsorted(3,:);

                % Best solution trackers
                X_best  = Alpha;
                X_best1 = X_best;

            else
                % Decreasing coefficient a
                a = 2*(1 - t/tmax);

                % Positions update
                Xnew = X(:,1:Nv); % preallocate

                for i = 1:Ni
                    Xi = X(i,1:Nv);

                    % With respect to Alpha
                    A1 = 2*a*rand(1,Nv) - a;
                    C1 = 2*rand(1,Nv);
                    Dalpha = abs(C1.*Alpha(1,1:Nv) - Xi);
                    X1 = Alpha(1,1:Nv) - A1.*Dalpha;

                    % With respect to Beta
                    A2 = 2*a*rand(1,Nv) - a;
                    C2 = 2*rand(1,Nv);
                    Dbeta = abs(C2.*Beta(1,1:Nv) - Xi);
                    X2 = Beta(1,1:Nv) - A2.*Dbeta;

                    % With respect to Delta
                    A3 = 2*a*rand(1,Nv) - a;
                    C3 = 2*rand(1,Nv);
                    Ddelta = abs(C3.*Delta(1,1:Nv) - Xi);
                    X3 = Delta(1,1:Nv) - A3.*Ddelta;

                    % New position
                    Xnew(i,:) = (X1 + X2 + X3)/3;
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
                X(:, Nv+1) = fo;

                % Update alpha, beta, delta
                Xsorted = sortrows(X, Nv+1);
                Alpha = Xsorted(1,:);
                Beta  = Xsorted(2,:);
                Delta = Xsorted(3,:);

                % Global best
                if Alpha(1,end) < X_best(1,end)
                    X_best(1,:) = Alpha;
                end

                % No improvement counter
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
            % Optional reporting
            % fprintf('Iteracion= %d, Error= %2.6f \n', t, X_best1(1,end));
            % report_conver(t+1,1) = X_best1(1,end);
        end
    end

    time = toc;
    Mejor_Sol(T,:)  = X_best1(1,1:Nv);  % Best decision vector
    Ploss_min(T,1)  = X_best1(1,end);   % Best objective value
    Times(T,1)      = time;

end

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

save results/GWO_Experimental
