function [error_f] = funcion_pid(var)
% Funcion objetivo para sintonizar un PID en un sistema de 2º orden
% var = [Kp, Ki, Kd]

    % -----------------------------------------------------------
    %           PARÁMETROS DEL CONTROLADOR PID
    % -----------------------------------------------------------
    Kp = var(1);
    Ki = var(2);
    Kd = var(3);

    % -----------------------------------------------------------
    %           PARÁMETROS DEL SISTEMA (masa-resorte-amortiguador)
    % -----------------------------------------------------------
    m = 1.0;      % masa [kg]
    c = 0.5;      % coef. amortiguamiento [Ns/m]
    k = 2.0;      % rigidez [N/m]

    % -----------------------------------------------------------
    %           CONFIGURACIÓN DE SIMULACIÓN
    % -----------------------------------------------------------
    tf = 5;                  % tiempo final [s]
    dt = 0.01;               % paso de integración
    tspan = 0:dt:tf;

    % referencia (ejemplo: escalón unitario desde t=0.5s)
    y_ref = @(t) (t >= 0.5);

    % condiciones iniciales [pos, vel, int_error]
    y0 = [0; 0; 0];

    % -----------------------------------------------------------
    %           SIMULACIÓN DEL SISTEMA
    % -----------------------------------------------------------
    [T, Y] = ode45(@(t,y) dynamics(t,y,m,c,k,Kp,Ki,Kd,y_ref), tspan, y0);

    pos = Y(:,1);
    e = arrayfun(y_ref,T) - pos;

    % -----------------------------------------------------------
    %           FUNCIÓN DE COSTO
    % -----------------------------------------------------------
    % Integral del error cuadrático (ISE)
    error_f = trapz(T, e.^2);

end

% ======================================================
% SUBFUNCIÓN: DINÁMICA DEL SISTEMA CON PID
% ======================================================
function dydt = dynamics(t,y,m,c,k,Kp,Ki,Kd,y_ref)
    % Estados:
    % y(1) = posición
    % y(2) = velocidad
    % y(3) = integral del error

    pos = y(1);
    vel = y(2);
    int_err = y(3);

    % Error de seguimiento
    ref = y_ref(t);
    err = ref - pos;

    % Control PID
    u = Kp*err + Ki*int_err - Kd*vel; % derivada sobre medición

    % Dinámica del sistema
    acc = (u - c*vel - k*pos)/m;

    % Derivadas de los estados
    dydt = [vel;
            acc;
            err];   % d/dt del error integral
end
