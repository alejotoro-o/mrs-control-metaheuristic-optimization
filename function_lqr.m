function [error_f] = function_lqr(var)

    % Controller weights (Q: states, R: inputs)
    % var = [Q1 Q2 Q3 R1 R2 R3]
    Q = diag(var(1:3));
    R = diag(var(4:6));

    % continuous LQR
    K = lqr(zeros(3),eye(3),Q,R);

    % Simulation parameters
    tf = 10;
    dt = 0.1;
    tspan = 0:dt:tf;

    % initial state
    q0 = [0;0;0];
    y0 = q0;

    % simulate
    [T,Y] = ode45(@(t,y) dynamics(t,y,K), tspan, y0);

    % Parse state history
    q = Y; % [x y theta]

    % compute reference at each time
    [x_ref, y_ref, theta_ref] = circular_reference(T);

    % position error
    pos_error = sqrt((x_ref - q(:,1)).^2 + (y_ref - q(:,2)).^2);

    % orientation error
    theta_error = wrapToPi(theta_ref - q(:,3));

    % Cost
    w_pos = 1.0;
    w_theta = 0.5;

    e_total = w_pos*pos_error.^2 + w_theta*theta_error.^2;
    error_f = trapz(T, e_total);

end

% System dynamics
function dydt = dynamics(t, q, K)

    % state: [x; y; theta]
    x = q(1); y = q(2); theta = q(3);

    % reference
    [xr, yr, thetar] = circular_reference(t);

    % state error
    e = [x - xr; y - yr; wrapToPi(theta - thetar)];

    % control law
    u = -K*e;

    % robot dynamics
    vx = u(1);
    vy = u(2);
    w = u(3);
    % Dynamics w.r.t. global frame
    dydt = [vx;
            vy;
            w];
end

% Circular trajectory reference
function [x_ref, y_ref, theta_ref] = circular_reference(t)

    R = 1.0;        % radius [m]
    w = 0.1;        % angular velocity [rad/s]

    x_ref = R*cos(w*t);
    y_ref = R*sin(w*t);

    % Orientation tangent to circle (robot heads along path)
    theta_ref = wrapToPi(w*t + pi/2);

    % Fixed orientation
    % theta_ref = 0; % Faces +x always
end
