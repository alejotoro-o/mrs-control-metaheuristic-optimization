close all
clear all
clc

params = [0.15 0.08 0.17]; % K = [Kd, Kalpha, Ktheta];
goal = [1.0 0.0]; % [d_goal, alpha_goal]
initial_conditions = [0.0 1.5; % [q1_0, q2_0]
                      0.0 0.5;
                      0.0 0.0];
sim_time = 20;
leader_input = [0.5;0.5;0.0]; % u_l = [v_x; v_y, w]

results = sim_leader_follower(params, goal, initial_conditions, sim_time, leader_input);

%% Plot
figure; hold on; grid on;
axis equal; % ✅ prevents arrow distortion
xlabel('X [m]'); ylabel('Y [m]');
title('Leader–Follower Simulation Results');

% Parameters
frame_scale = 0.2;
lw = 1;

% Plot trajectories
leader = plot(results.q1(:,1), results.q1(:,2), 'b-');
follower = plot(results.q2(:,1), results.q2(:,2), 'r-');

% Helper to draw heading frame
draw_frame = @(x,y,th,color) ...
    [ quiver(x, y, frame_scale*cos(th), frame_scale*sin(th), 0, ...
             'Color',color, 'LineWidth',1.2, 'MaxHeadSize',1.2), ...
      quiver(x, y, -frame_scale*sin(th), frame_scale*cos(th), 0, ...
             'Color',color, 'LineWidth',1.2, 'MaxHeadSize',1.2) ];

% Draw frames only at beginning & end
% Leader
draw_frame(results.q1(1,1), results.q1(1,2), results.q1(1,3), 'b');
draw_frame(results.q1(end,1), results.q1(end,2), results.q1(end,3), 'b');

% Follower
draw_frame(results.q2(1,1), results.q2(1,2), results.q2(1,3), 'r');
draw_frame(results.q2(end,1), results.q2(end,2), results.q2(end,3), 'r');

legend([leader, follower], {'Leader', 'Follower'}, 'Location', 'best');
hold off;