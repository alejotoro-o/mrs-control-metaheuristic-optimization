function [results] = sim_leader_follower(params, goal, initial_conditions, sim_time, leader_input)

    arguments
        params (1,3) double = [1.0 1.0 1.0] % K = [Kd, Kalpha, Ktheta];
        goal (1,2) double = [1.0 0.0] % [d_goal, alpha_goal]
        initial_conditions (3,2) double = [0.0 0.5; % [q1_0, q2_0]
                                           0.0 0.5;
                                           0.0 0.0]
        sim_time double = 20
        leader_input (3,1) double = [0.0;0.0;0.0] % u_l = [v_x; v_y, w]
    end

    % Controller params
    Kd = params(1);
    Kalpha = params(2);
    Ktheta = params(3);
    K = [Kd, Kalpha, Ktheta];

    % System params
    d_goal = goal(1);          % Desired distance between robots [m]
    alpha_goal = goal(2);      % Desired angle between robots [rad]
    tf = sim_time;               % Sim time [s]
    dt = 0.01;             % dt
    tspan = 0:dt:tf;

    % Initial conditions
    % q1 = [x1, y1, theta1] (leader)
    % q2 = [x2, y2, theta2] (follower)
    q1_0 = initial_conditions(:,1);
    q2_0 = initial_conditions(:,2);
    y0 = [q1_0; q2_0]; 

    % System sim
    [T, Y] = ode45(@(t,y) dynamics(t,y,K,d_goal,alpha_goal, leader_input), tspan, y0);

    % Leader and follower pose
    q1 = Y(:,1:3);
    q2 = Y(:,4:6);

    % Tracking error of distance between robots
    d = sqrt((q1(:,1) - q2(:,1)).^2 + (q1(:,2) - q2(:,2)).^2);
    d_error = d - d_goal;

    % Orientation error
    theta_F_goal = q1(:,3) - pi;
    theta_error = wrapToPi(q2(:,3) - theta_F_goal);

    % Cost function
    % error weights
    w_d = 1.0;
    w_theta = 0.5;

    e_total = w_d*(d_error.^2) + w_theta*(theta_error.^2);
    results.error_f = trapz(T, e_total);
    results.tspan = tspan;
    results.q1 = q1;
    results.q2 = q2;   

end

% System dynamics
function dydt = dynamics(t,y,K,d_goal,alpha_goal,leader_input)
    % Get pose
    q1 = y(1:3);
    q2 = y(4:6);

    % Leader parameters (circular trajectory)
    u_l = leader_input;
    u_l_G = [u_l(1)*cos(q1(3));
           u_l(2)*sin(q1(3));
           u_l(3)];

    % Leader dynamics
    dq1 = [u_l_G(1);
           u_l_G(2);
           u_l_G(3)];

    % Follower input using controller
    u_f = follower_control(q1,q2,u_l,K,d_goal,alpha_goal);

    % Follower dynamics w.r.t. global frame
    dq2 = [u_f(1)*cos(q2(3)) - u_f(2)*sin(q2(3));
           u_f(1)*sin(q2(3)) + u_f(2)*cos(q2(3));
           u_f(3)];

    dydt = [dq1; dq2];
end

% Controller
function u_f = follower_control(q1, q2, u_l, K, d_goal, alpha_goal)
    x1 = q1(1); y1 = q1(2); th1 = q1(3);
    x2 = q2(1); y2 = q2(2); th2 = q2(3);

    % Distance and relative angle
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
    alpha = th2 - atan2((y1 - y2),(x1 - x2));
    alpha = wrapToPi(alpha);
    gamma = th1 - th2 + alpha;
    gamma = wrapToPi(gamma);

    % Matrices A and B
    A = [cos(gamma), -sin(gamma), 0;
        -(1/d)*sin(gamma), -(1/d)*cos(gamma), 0;
        0, 0, 0];

    B = [cos(alpha), -sin(alpha), 0;
        -(1/d)*sin(alpha), -(1/d)*cos(alpha), -1;
        0, 0, -1];

    % Errors
    alpha_error = wrapToPi(alpha - alpha_goal);
    theta_error = wrapToPi(th2 - (th1 - pi));
    d_error = d - d_goal;

    % Vector p_d
    p_d = [-K(1)*d_error;
           -K(2)*alpha_error;
            u_l(3) - K(3)*theta_error];

    % Control law
    u_f = B \ (A*u_l - p_d);
end