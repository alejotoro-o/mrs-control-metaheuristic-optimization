function PSO
clc; clear; close all; format longG;
cfg = common_load_config();

% Hiperparámetros
Ni   = 50;
wMax = 0.95; wMin = 0.35; C1 = 1.4; C2 = 1.4;
vMax = 0.2*(cfg.xmax - cfg.xmin); vMin = -vMax;

Mejor_Sol = nan(cfg.T, cfg.Nv);
BestObj   = nan(cfg.T, 1);
Times     = nan(cfg.T, 1);
report_conver_all = cell(cfg.T,1);

for run = 1:cfg.T
    tStart = tic;

    X = cfg.xmin + (cfg.xmax - cfg.xmin).*rand(Ni, cfg.Nv);
    V = zeros(Ni, cfg.Nv);
    f = zeros(Ni,1);
    parfor i=1:Ni, f(i)=cfg.objfun(X(i,:)); end
    Pbest_pos = X; Pbest_val = f;
    [gbest_val, gidx] = min(f);
    Gbest_pos = X(gidx,:);

    rc = nan(cfg.tmax+1,1); rc(1)=gbest_val; tx=0;
    for t=1:cfg.tmax
        w = wMax - (wMax - wMin)*(t/cfg.tmax);
        rp = rand(Ni,cfg.Nv); rg = rand(Ni,cfg.Nv);
        V = w.*V + C1*rp.*(Pbest_pos - X) + C2*rg.*(Gbest_pos - X);
        V = max(min(V, vMax), vMin);
        X = min(max(X + V, cfg.xmin), cfg.xmax);

        parfor i=1:Ni, f(i)=cfg.objfun(X(i,:)); end
        improved = f < Pbest_val;
        Pbest_pos(improved,:) = X(improved,:);
        Pbest_val(improved)   = f(improved);

        [curr, idx] = min(Pbest_val);
        if curr < gbest_val, gbest_val=curr; Gbest_pos=Pbest_pos(idx,:); tx=0; else, tx=tx+1; end
        rc(t+1)=gbest_val; if tx>=cfg.iterno, break; end
    end

    li=find(~isnan(rc),1,'last'); if isempty(li), li=1; end
    report_conver_all{run}=rc(1:li);

    Times(run,1)   = toc(tStart);
    Mejor_Sol(run,:)=Gbest_pos;
    BestObj(run,1) = gbest_val;
end

Ploss_min = BestObj; % alias
[MediaTime,DesviacionestTime] = deal(mean(Times), std(Times));
[MeanBest,StdBest]            = deal(mean(BestObj), std(BestObj));
PorcentDesTime = 100*DesviacionestTime/MediaTime;
PorcentDesSol  = 100*StdBest/MeanBest;
[bestsolution1,ix] = min(BestObj); bestsolution = Mejor_Sol(ix,:);
Solution = [bestsolution bestsolution1 MeanBest PorcentDesSol MediaTime];

report_conver = report_conver_all;
save results/PSO_Experimental Mejor_Sol BestObj Ploss_min Times Solution report_conver
end
