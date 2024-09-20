% parameters
num_polygons = 4;
sigma = 3;
sides = 3;
particles_per_side = 1;
total_particles = sides * particles_per_side;
m = ones(1, total_particles);
k_int = 10000;
b_trans = 2;
b_ang = 2;

dt = 0.001;
max_t = 100;

% BINARY SEARCH PARAMS
L = 15;
r_scale = 0.99;
P_t = 10;
P_l = 0.995 * P_t;
P_h = 1.005 * P_t;
L_l = -1;
L_h = -1;
L_tol = 0.00000000001;

% Define xmin, xmax, ymin, ymax
xmin = -L/2;
xmax = L/2;
ymin = -L/2;
ymax = L/2;

polygons = cell(1, num_polygons);

% generate initial conditions without overlap
disp('Generating initial conditions...');
pressure = 1;
while pressure > 0
    for i = 1:num_polygons
        q = [rand()*(xmax-xmin)+xmin, rand()*(ymax-ymin)+ymin];
        v = rand(1, 2)*5 - 2.5;
        w = 0;
        polygons{i} = regular_polygon(sigma, sides, m, q, v, w, particles_per_side);
    end
    box = periodic_box_polygon(polygons, k_int, b_trans, b_ang, -L/2, L/2, -L/2, L/2);
    pressure = box.iterate_time(dt);
end
disp('Found Starting Config');

old_state = box;

% Define E_thresh, num_steps, and colors
E_thresh = 1e-6;
num_steps = round(max_t / dt);
colors = lines(num_polygons);

figure;
subplot(2, 1, 1);
title('Simulation');
xlabel('X');
ylabel('Y');
axis equal;

subplot(2, 1, 2);
line_kinetic_trans = animatedline('Color', 'r');
line_kinetic_rot = animatedline('Color', 'g');
line_potential = animatedline('Color', 'b');
line_total = animatedline('Color', 'k');
title('Energy');
xlabel('Time');
ylabel('Energy');
legend('Kinetic (Trans)', 'Kinetic (Rot)', 'Potential', 'Total');

progress_bar = waitbar(0, 'Running binary search...');
while true
    [t_values, q_series, v_series, theta_series, w_series,...
    pressure_series, vertices_series, state, polygons, particles_per_side,...
    xmin_series, xmax_series, ymin_series, ymax_series] =...
    energy_minimize(box, L, L_l, L_h, E_thresh, dt, max_t);

    L_prior = L;
    P = pressure_series(end);
    if P < P_l
        old_state = state;
        L_old = L;
        L_h = L;
        if L_l > 0
            L = (L_h + L_l)/2;
            L_l = -1;
        else
            L = L * r_scale;
        end
    else
        if P>P_h
            L_l = L;
            state = old_state;
            L_prior = L_old;
            L = (L_h + L_l)/2;
        else
            for step = 1:num_steps
                if mod(step,1000) == 0
                    % Update the simulation plot
                    subplot(2, 1, 1);
                    cla;
                    rectangle('Position', [xmin_series(step), ymin_series(step),...
                        xmax_series(step)-xmin_series(step),...
                        ymax_series(step)-ymin_series(step)], ...
                        'EdgeColor', 'k', 'LineWidth', 2);
                    
                    for polygon = 1:num_polygons
                        sigma = polygons{polygon}.sigma;
                        vertices = vertices_series{polygon, step};
                        plot([vertices(:,1); vertices(1,1)], [vertices(:,2); vertices(1,2)], ...
                            '-', 'Color', colors(polygon,:), 'LineWidth', 2);
                        
                        for j = 1:size(vertices, 1)
                            rectangle('Position', [vertices(j,1)-sigma/2, vertices(j,2)-sigma/2, sigma, sigma], ...
                                'Curvature', [1, 1], 'FaceColor', colors(polygon,:));
                        end
                    end
                    
                    title(sprintf('Time: %.2f, Pressure: %.2f', t_values(step), pressure_series(step)));
                    
                    % Update the energy plot
                    subplot(2, 1, 2);
                    set(line_kinetic_trans, 'XData', t_values(1:step), 'YData', E_series(1:step, 1));
                    set(line_kinetic_rot, 'XData', t_values(1:step), 'YData', E_series(1:step, 2));
                    set(line_potential, 'XData', t_values(1:step), 'YData', E_series(1:step, 3));
                    set(line_total, 'XData', t_values(1:step), 'YData', sum(E_series(1:step,:), 2));
                end
    
                drawnow;
            end
            break;
        end
    end
    waitbar((L_h - L) / (L_h - L_l), progress_bar);
end
close(progress_bar);

function [t_values, q_series, v_series, theta_series, w_series,...
    pressure_series, vertices_series, box, polygons, particles_per_side,...
    xmin_series, xmax_series, ymin_series, ymax_series] =...
    energy_minimize(system, L, L_l, L_h, E_thresh, dt, max_t)
    
    box = system.copy();
    num_polygons = numel(box.polygons);

    max_t = round(max_t / dt) * dt;
    t_values = 0:dt:max_t;
    num_steps = length(t_values);
    q_series = cell(1, num_polygons);
    v_series = cell(1, num_polygons);
    theta_series = cell(1, num_polygons);
    w_series = cell(1, num_polygons);
    pressure_series = zeros(1, num_steps);
    vertices_series = cell(num_polygons, num_steps);
    xmin_series = ones(1, num_steps) * (-L/2);
    xmax_series = ones(1, num_steps) * (L/2);
    ymin_series = ones(1, num_steps) * (-L/2);
    ymax_series = ones(1, num_steps) * (L/2);
    E_series = zeros(num_steps, 3);
    
    progress_bar = waitbar(0, 'Running energy minimization...');
    for step = 1:num_steps
        pressure_series(step) = box.iterate_time(dt);
    
        xmin_series(step) = box.xmin;
        xmax_series(step) = box.xmax;
        ymin_series(step) = box.ymin;
        ymax_series(step) = box.ymax;
    
        [E_series(step,1),E_series(step,2),E_series(step,3)] = box.get_energy();
        for polygon = 1:num_polygons
            poly = box.polygons(polygon);
            q_series{polygon}(step, :) = box.apply_pbc2d(poly.q);
            v_series{polygon}(step, :) = poly.v;
            theta_series{polygon}(step) = poly.theta;
            w_series{polygon}(step) = poly.w;
            vertices_series{polygon, step} = poly.vertices_relative + q_series{polygon}(step, :);
        end

        if E_series(step,3) <= E_thresh
            polygons = box.polygons;
            particles_per_side = polygons(1).particles_per_side;
            break;
        end
        waitbar(step / num_steps, progress_bar);
    end
    if L_l > 0 && abs((L_h / L_l) - 1) < L_tol
        error('No jammed packing found.')
    end
    close(progress_bar);
    disp('Energy minimization complete');

end