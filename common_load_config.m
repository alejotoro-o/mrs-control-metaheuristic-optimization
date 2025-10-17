function cfg = common_load_config()
% Carga la configuración global y devuelve un struct listo para usar.

S = load('global_run_config.mat','cfg');
if ~isfield(S,'cfg')
    error('global_run_config.mat no contiene el struct cfg.');
end
cfg = S.cfg;

% Construye el handle de la FO a partir del nombre
if ~isfield(cfg,'objfun_name') || isempty(cfg.objfun_name)
    error('cfg.objfun_name no definido.');
end
if exist(cfg.objfun_name,'file')==0
    error('No se encuentra la función objetivo "%s" en el path.', cfg.objfun_name);
end
cfg.objfun = str2func(cfg.objfun_name);

% Comprobaciones básicas
fields = {'Nv','xmin','xmax','T','tmax','iterno'};
for k=1:numel(fields)
    if ~isfield(cfg,fields{k}), error('cfg.%s no definido.',fields{k}); end
end
if numel(cfg.xmin)~=cfg.Nv || numel(cfg.xmax)~=cfg.Nv
    error('Longitud de xmin/xmax debe ser Nv.');
end
end
