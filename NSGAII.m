clc
clear all

for T = 1:1
    tic;

    %% NSGA-II (estructura compatible con sus scripts previos)
    Ni     = 50;      % Tamaño de población (debe ser par)
    Nv     = 6;       % Número de variables
    tmax   = 1000;    % Máximo de generaciones
    iterno = 10;      % Ventana de no-mejora para parada temprana
    tx     = 0;

    xmin = deg2rad([-180 -180 -180 -180 -180 -180]);
    xmax = deg2rad([ 180  180  180  180  180  180]);

    % Parámetros genéticos (SBX + mutación polinómica)
    pc    = 0.9;      % probabilidad de cruce
    pm    = 1/Nv;     % probabilidad de mutación por variable
    etaC  = 15;       % índice de distribución para SBX
    etaM  = 20;       % índice de distribución para mutación polinómica

    % ---------- Inicialización ----------
    % Población: Ni x Nv, más columna para FO
    X = xmin.*ones(Ni,Nv) + (xmax - xmin).*rand(Ni,Nv);

    % Evaluación inicial
    f = zeros(Ni,1);
    parfor i = 1:Ni
        f(i) = funcionCD_UR3_Jac_orientacion(X(i,:));
    end

    % Para mantener estructura de “mejor”
    [bestf, ibest] = min(f);
    X_best  = [X(ibest,:), bestf];   % [x*, f*]
    X_best1 = X_best;

    % ---------- Ciclo evolutivo ----------
    for t = 0:tmax
        if t > 0
            % 1) Selección por torneo binario con Rank y Crowding
            % Con 1 objetivo, Rank = 1 para todos; se desempata por crowding = diversidad
            rank = ones(Ni,1);
            crowd = computeCrowding(f); % mayor crowd preferible

            % Padres (índices) por torneos
            parents = zeros(Ni,1);
            for k = 1:Ni
                a = randi(Ni); b = randi(Ni);
                if rank(a) < rank(b)
                    parents(k) = a;
                elseif rank(b) < rank(a)
                    parents(k) = b;
                else
                    % mismo rank, usar crowding (mayor es mejor)
                    if crowd(a) > crowd(b)
                        parents(k) = a;
                    else
                        parents(k) = b;
                    end
                end
            end

            % 2) Cruce SBX + mutación polinómica para generar descendencia
            Y = zeros(Ni, Nv);
            for k = 1:2:Ni
                p1 = X(parents(k),   :);
                p2 = X(parents(k+1), :);

                c1 = p1; c2 = p2;
                if rand < pc
                    [c1, c2] = sbx_crossover(p1, p2, xmin, xmax, etaC);
                end

                % Mutación polinómica
                c1 = poly_mutation(c1, xmin, xmax, pm, etaM);
                c2 = poly_mutation(c2, xmin, xmax, pm, etaM);

                % Respetar cotas
                c1 = min(max(c1, xmin), xmax);
                c2 = min(max(c2, xmin), xmax);

                Y(k,:)   = c1;
                Y(k+1,:) = c2;
            end

            % 3) Evaluar descendencia
            fy = zeros(Ni,1);
            parfor i = 1:Ni
                fy(i) = funcionCD_UR3_Jac_orientacion(Y(i,:));
            end

            % 4) Unión, ordenamiento no-dominado y poda por crowding
            % Con una sola FO, basta ordenar por f ascendente, pero se implementa
            % el esquema estándar para que pueda pasar a multiobjetivo sin reescribir.
            % Unir
            U  = [X; Y];
            fu = [f; fy];

            % Non-dominated sorting “trivial” para 1 objetivo:
            % crear frentes simplemente por orden ascendente
            [fu_sorted, idx] = sort(fu, 'ascend');
            U_sorted = U(idx,:);

            % Quedarse con los mejores Ni aplicando crowding en el último corte
            % Como es 1 objetivo, crowding se computa en fu_sorted
            % Si hay empate en el borde del corte, usar crowding
            X_next = zeros(Ni, Nv);
            f_next = zeros(Ni, 1);

            % Si múltiples con igual valor en el borde, aplicar crowding en ese grupo
            % Implementación simple: tomar los primeros Ni, si queremos refinar:
            if length(fu_sorted) > Ni
                % Buscar todos con el mismo valor que el último incluido
                cutoff_val = fu_sorted(Ni);
                cut_idx = fu_sorted <= cutoff_val;
                K = find(cut_idx);
                if numel(K) > Ni
                    % demasiados con el mismo valor, usar crowding local
                    subU = U_sorted(K,:);
                    subf = fu_sorted(K);
                    subCrowd = computeCrowding(subf);
                    % ordenar subgrupo por crowding descendente
                    [~, ordC] = sort(subCrowd, 'descend');
                    pick = ordC(1:Ni);
                    X_next = subU(pick,:);
                    f_next = subf(pick);
                else
                    % tomar todos los previos a cutoff y completar hasta Ni
                    X_next(1:numel(K),:) = U_sorted(1:numel(K),:);
                    f_next(1:numel(K))   = fu_sorted(1:numel(K));

                    % si faltan, completar con los siguientes en orden
                    if numel(K) < Ni
                        rest = (numel(K)+1):Ni;
                        X_next(rest,:) = U_sorted(numel(K)+1:Ni,:);
                        f_next(rest)   = fu_sorted(numel(K)+1:Ni);
                    end
                end
            else
                X_next = U_sorted(1:Ni,:);
                f_next = fu_sorted(1:Ni);
            end

            % Actualizar población
            X = X_next;
            f = f_next;

            % 5) Actualizar mejor global
            [curr_best, ic] = min(f);
            if curr_best < X_best(1,end)
                X_best = [X(ic,:), curr_best];
            end

            % 6) Parada por no-mejora
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

        % % Opcional:
        % fprintf('Gen %4d  f*=%.6e\n', t, X_best1(1,end));
    end

    time = toc;

    % --- Salidas con su misma convención ---
    Mejor_Sol(T,:) = X_best1(1,1:Nv); % mejor vector decisión
    Ploss_min(T,1) = X_best1(1,end);  % mejor FO
    Times(T,1)     = time;
end

% Estadísticos finales (mismo formato que sus scripts)
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

save results/NSGAII_Experimental


%% ===================== FUNCIONES AUXILIARES =====================

% Crowding distance 1D para una sola FO
function cd = computeCrowding(f)
    n  = numel(f);
    cd = zeros(n,1);
    [fs, idx] = sort(f, 'ascend');
    cd(idx(1))   = inf;
    cd(idx(end)) = inf;
    fmin = fs(1); fmax = fs(end);
    if fmax > fmin
        for k = 2:n-1
            cd(idx(k)) = (fs(k+1) - fs(k-1)) / (fmax - fmin);
        end
    else
        % todos iguales, dar pequeña diversidad aleatoria
        cd = rand(n,1)*1e-6;
    end
end

% SBX crossover (Simulated Binary Crossover) con cotas
function [c1, c2] = sbx_crossover(p1, p2, xmin, xmax, etaC)
    Nv = numel(p1);
    c1 = zeros(1,Nv);
    c2 = zeros(1,Nv);
    for j = 1:Nv
        u = rand;
        if abs(p1(j) - p2(j)) < eps
            c1(j) = p1(j);
            c2(j) = p2(j);
        else
            x1 = min(p1(j), p2(j));
            x2 = max(p1(j), p2(j));
            lb = xmin(j); ub = xmax(j);

            beta  = 1 + (2*(x1 - lb)/(x2 - x1));
            alpha = 2 - beta^(-(etaC+1));
            if u <= 1/alpha
                betaq = (u*alpha)^(1/(etaC+1));
            else
                betaq = (1/(2 - u*alpha))^(1/(etaC+1));
            end
            c1j = 0.5*((x1 + x2) - betaq*(x2 - x1));

            beta  = 1 + (2*(ub - x2)/(x2 - x1));
            alpha = 2 - beta^(-(etaC+1));
            if u <= 1/alpha
                betaq = (u*alpha)^(1/(etaC+1));
            else
                betaq = (1/(2 - u*alpha))^(1/(etaC+1));
            end
            c2j = 0.5*((x1 + x2) + betaq*(x2 - x1));

            % En cota
            c1(j) = min(max(c1j, lb), ub);
            c2(j) = min(max(c2j, lb), ub);
        end
    end
end

% Mutación polinómica con cotas
function y = poly_mutation(x, xmin, xmax, pm, etaM)
    Nv = numel(x);
    y  = x;
    for j = 1:Nv
        if rand < pm
            lb = xmin(j); ub = xmax(j);
            if ub - lb < eps
                y(j) = min(max(x(j), lb), ub);
                continue;
            end
            delta1 = (x(j) - lb)/(ub - lb);
            delta2 = (ub - x(j))/(ub - lb);
            r = rand;
            mut_pow = 1/(etaM + 1);
            if r < 0.5
                xy = 1 - delta1;
                val = 2*r + (1 - 2*r)*(xy^(etaM + 1));
                deltaq = val^mut_pow - 1;
            else
                xy = 1 - delta2;
                val = 2*(1 - r) + 2*(r - 0.5)*(xy^(etaM + 1));
                deltaq = 1 - val^mut_pow;
            end
            y(j) = x(j) + deltaq*(ub - lb);
            y(j) = min(max(y(j), lb), ub);
        end
    end
end
