function GA
clc; clear; close all; format longG;
cfg = common_load_config();

Ni = 50; pc = 0.9; eta_c = 15; pm = 1/cfg.Nv; eta_m = 20;

Mejor_Sol=nan(cfg.T,cfg.Nv); BestObj=nan(cfg.T,1); Times=nan(cfg.T,1);
report_conver_all=cell(cfg.T,1);

for run=1:cfg.T
    tStart=tic;

    Pop = cfg.xmin + (cfg.xmax - cfg.xmin).*rand(Ni,cfg.Nv);
    f = zeros(Ni,1); parfor i=1:Ni, f(i)=cfg.objfun(Pop(i,:)); end
    [best_val, ib]=min(f); best_pos=Pop(ib,:);
    rc=nan(cfg.tmax+1,1); rc(1)=best_val; tx=0;

    for t=1:cfg.tmax
        % torneo
        M=zeros(Ni,cfg.Nv);
        for k=1:Ni
            i1=randi(Ni); i2=randi(Ni);
            M(k,:)=Pop( ternary(f(i1)<f(i2), i1, i2), :);
        end
        % crossover
        Off=M;
        for k=1:2:Ni-1
            if rand<pc
                [c1,c2]=sbx(M(k,:),M(k+1,:),cfg.xmin,cfg.xmax,eta_c);
                Off(k,:)=c1; Off(k+1,:)=c2;
            end
        end
        % mutación
        for k=1:Ni
            Off(k,:)=poly_mut(Off(k,:),cfg.xmin,cfg.xmax,pm,eta_m);
        end
        % evaluar
        fOff=zeros(Ni,1); parfor i=1:Ni, fOff(i)=cfg.objfun(Off(i,:)); end
        % elitismo simple
        U=[Pop;Off]; fU=[f;fOff];
        [f_sorted, idx]=sort(fU,'ascend'); Pop=U(idx(1:Ni),:); f=f_sorted(1:Ni);

        if f(1)<best_val, best_val=f(1); best_pos=Pop(1,:); tx=0; else, tx=tx+1; end
        rc(t+1)=best_val; if tx>=cfg.iterno, break; end
    end

    li=find(~isnan(rc),1,'last'); if isempty(li), li=1; end
    report_conver_all{run}=rc(1:li);

    Times(run,1)=toc(tStart);
    Mejor_Sol(run,:)=best_pos;
    BestObj(run,1)=best_val;
end

Ploss_min=BestObj;
[MediaTime,DesviacionestTime]=deal(mean(Times),std(Times));
[MeanBest,StdBest]=deal(mean(BestObj),std(BestObj));
PorcentDesTime=100*DesviacionestTime/MediaTime;
PorcentDesSol =100*StdBest/MeanBest;
[bestsolution1,ix]=min(BestObj); bestsolution=Mejor_Sol(ix,:);
Solution=[bestsolution bestsolution1 MeanBest PorcentDesSol MediaTime];

report_conver=report_conver_all;
save results/GA_Experimental Mejor_Sol BestObj Ploss_min Times Solution report_conver

%% helpers
function [c1,c2]=sbx(p1,p2,lb,ub,eta)
u=rand(size(p1)); beta=zeros(size(p1));
beta(u<=0.5)=(2*u(u<=0.5)).^(1/(eta+1));
beta(u>0.5) =(2*(1-u(u>0.5))).^(-1/(eta+1));
c1=0.5*((1+beta).*p1 + (1-beta).*p2);
c2=0.5*((1-beta).*p1 + (1+beta).*p2);
c1=min(max(c1,lb),ub); c2=min(max(c2,lb),ub);
end
function y=poly_mut(x,lb,ub,pm,eta)
y=x;
for j=1:numel(x)
    if rand<pm
        delta1=(x(j)-lb(j))/(ub(j)-lb(j));
        delta2=(ub(j)-x(j))/(ub(j)-lb(j));
        u=rand; mut_pow=1/(eta+1);
        if u<=0.5
            xy=1-delta1; val=2*u + (1-2*u)*(xy^(eta+1)); deltaq=val^mut_pow - 1;
        else
            xy=1-delta2; val=2*(1-u) + 2*(u-0.5)*(xy^(eta+1)); deltaq=1 - val^mut_pow;
        end
        y(j)=x(j) + deltaq*(ub(j)-lb(j));
        y(j)=min(max(y(j),lb(j)),ub(j));
    end
end
end
function r=ternary(c,a,b), if c, r=a; else, r=b; end; end
end
