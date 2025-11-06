%% RUNNER ROBUSTO: ejecuta todos los algoritmos metaheurísticos
clear; clc;

%% --- CONFIGURACIÓN GLOBAL ---
cfg = struct();
cfg.objfun_name = 'function_lqr'; % nombre de la función objetivo
cfg.Nv          = 6;                                % número de variables
cfg.xmin        = [0 0 0 0 0 0] + 1;
cfg.xmax        = [10 10 10 10 10 10];

% Corridas e iteraciones
cfg.T           = 10;     % número de corridas por algoritmo
cfg.tmax        = 10;   % iteraciones máximas
cfg.iterno      = 1;    % paciencia sin mejora 

% Poblaciones por algoritmo (puedes ajustar por algoritmo dentro del .m)
cfg.Ni_default  = 50;

% Guarda la configuración en un MAT que leerán los algoritmos
save('global_run_config.mat','cfg');

%% --- Lista de algoritmos (ajusta nombres si difieren) ---
%algos = {'PSO','GA','MVO','GWO','CMAES','JADE','NSGAII','WOA','SSA','OCPHSOLANA'};

algos = {'PSO','GA'};
workdir = pwd;

%% --- Comprobación de archivos ---
req = [{'common_load_config.m','global_run_config.mat', cfg.objfun_name} strcat(algos, '.m')];
missing = req(~cellfun(@(f) exist(f,'file')>0, req));
if ~isempty(missing)
    fprintf(2,'Missing files:\n'); fprintf(2,'  - %s\n', missing{:});
    error('Add the missing files or fix names in "algos".');
end

%% --- Carpeta de logs ---
logdir = fullfile(workdir,'logs_run_all');
if ~exist(logdir,'dir'), mkdir(logdir); end

%% --- Carpeta de resultados ---
resdir = fullfile(workdir,'results');
if ~exist(resdir,'dir'), mkdir(resdir); end

%% --- Detectar soporte de -batch ---
hasBatch = ~verLessThan('matlab','9.6'); % 9.6 = R2019a

%% --- Detectar ejecutable de MATLAB ---
if ispc, matlabExe = fullfile(matlabroot,'bin','matlab.exe');
else,     matlabExe = fullfile(matlabroot,'bin','matlab'); end

%% --- Ejecutar ---
fprintf('=== Launching algorithms (%s mode) ===\n', ternary(hasBatch,'-batch','-r'));

results = struct('algo',algos,'status',[],'elapsed',[],'outfile',[],'logfile',[]);
for i = 1:numel(algos)
    algo = algos{i};
    logf = fullfile(logdir, sprintf('log_%s_%s.txt', algo, datestr(now,'yyyymmdd_HHMMSS')));
    fprintf('\n[%2d/%2d] %s ...\n', i, numel(algos), algo);

    if hasBatch
        batchCmd = sprintf('cd(''%s''); %s', workdir, algo);
        cmd = sprintf('"%s" -batch "%s" > "%s" 2>&1', matlabExe, batchCmd, logf);
    else
        wdEsc = strrep(workdir, '''', '''''');
        rCmd = sprintf(['try, cd(''%s''); %s;', ...
                        'catch ME, disp(getReport(ME,''extended'')); exit(1); end; exit'], wdEsc, algo);
        cmd = sprintf('"%s" -nosplash -nodesktop -r "%s" > "%s" 2>&1', matlabExe, rCmd, logf);
    end

    t0 = tic;
    code = system(cmd);
    elapsed = toc(t0);

    % Archivo de salida esperado
    expectedOut = ['results/' algo '_Experimental.mat'];
    if ~exist(expectedOut,'file')
        D = dir('results/*_Experimental.mat');
        if ~isempty(D)
            [~,idxNew] = max([D.datenum]);
            expectedOut = D(idxNew).name;
        else
            expectedOut = '';
        end
    end

    results(i).status  = code;
    results(i).elapsed = elapsed;
    results(i).outfile = expectedOut;
    results(i).logfile = logf;

    if code == 0
        fprintf('  Done in %.1f s. Output: %s\n  Log: %s\n', elapsed, expectedOut, logf);
    else
        fprintf(2,'  ERROR (code %d). See log: %s\n', code, logf);
    end
end

%% --- Resumen ---
fprintf('\n=== Summary ===\n');
for i = 1:numel(results)
    tag = ternary(results(i).status==0,'OK',sprintf('ERR(%d)',results(i).status));
    fprintf('%-8s | %-8s | %6.1f s | %s\n', results(i).algo, tag, results(i).elapsed, results(i).outfile);
end

save('results/run_all_summary.mat','results');
fprintf('\nSaved summary to run_all_summary.mat\n');

%% --- Funciones auxiliares ---
function out = ternary(cond,a,b)
    if cond, out=a; else, out=b; end
end

