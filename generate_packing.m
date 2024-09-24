% parameters
num_polygons = 10;
sigma = 1;
sides = 3;
particles_per_side = 2;
total_particles = sides * particles_per_side;
m = ones(1, total_particles);
k_int = 10000;
b_trans = 10;
b_ang = 10;
dt = 0.001;
max_t = 20;

% BINARY SEARCH PARAMS
L = 10;
L_old = L;
r_scale = 0.99;
P_t = 10; 
P_l = 0.995 * P_t;
P_h = 1.005 * P_t;
L_l = -1;
L_h = -1;
L_tol = .01;

polygons = cell(1, num_polygons);

% Generate initial conditions without overlap
pressure = 101; 
while pressure > 100
    % Calculate grid size based on polygon size and number of polygons
    grid_size = ceil(sqrt(num_polygons));
    polygon_size = sigma * particles_per_side;
    spacing = (L - grid_size * polygon_size) / (grid_size + 1);
    
    for i = 1:num_polygons
        row = floor((i - 1) / grid_size);
        col = mod(i - 1, grid_size);
        q = [-L/2 + (col + 1) * spacing + col * polygon_size, ...
             -L/2 + (row + 1) * spacing + row * polygon_size];
        v = rand(1,2)*1 - 0.5;
        w = 0;
        polygons{i} = regular_polygon(sigma, sides, m, q, v, w, particles_per_side);
    end
    box = periodic_box_polygon(polygons, k_int, b_trans, b_ang, -L/2, L/2, -L/2, L/2);
    pressure = box.iterate_time(dt);
    
    if pressure > 100
        % increase box size if pressure is too high
        L_old = L;
        L = L * 1.1;
        
        % update grid spacing based on new box size
        spacing = (L - grid_size * polygon_size) / (grid_size + 1);
        
        % update polygon positions based on new spacing
        for i = 1:num_polygons
            row = floor((i - 1) / grid_size);
            col = mod(i - 1, grid_size);
            q_old = polygons{i}.q;
            q_new = [-L/2 + (col + 1) * spacing + col * polygon_size, ...
                     -L/2 + (row + 1) * spacing + row * polygon_size];
            polygons{i}.q = q_new;
            polygons{i}.vertices_relative = polygons{i}.vertices_relative + (q_new - q_old);
        end
        
        box.xmin = -L/2;
        box.xmax = L/2;
        box.ymin = -L/2;
        box.ymax = L/2;
    end
    
    disp(pressure)
end

old_state = box;
E_thresh = 1e-6;
num_steps = round(max_t / dt);

colors = lines(num_polygons);

fprintf('Starting binary search...\n');
iteration = 0;
while true
    iteration = iteration + 1;    
    [t_values, q_series, v_series, theta_series, w_series,...
    pressure_series, vertices_series, state, polygons, particles_per_side,...
    xmin_series, xmax_series, ymin_series, ymax_series] =...
    energy_minimize(box, L, E_thresh, dt, max_t);
    
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
            box = old_state.copy();
            L = (L_h + L_l)/2;
        else
            break;
        end
    end
    % RESIZE SYSTEM HERE
    box.xmin = -L/2;
    box.xmax = L/2;
    box.ymin = -L/2;
    box.ymax = L/2;

    if L_l > 0 && abs((L_h / L_l) - 1) < L_tol
        error('No jammed packing found.')
    end 
    fprintf('Binary search iteration: %d, Current L: %.6f, Current P: %.6f\n', iteration, L, P);
end
fprintf('Binary search complete.\n');

function [t_values, q_series, v_series, theta_series, w_series,...
    pressure_series, vertices_series, box, polygons, particles_per_side,...
    xmin_series, xmax_series, ymin_series, ymax_series] =...
    energy_minimize(system, L, E_thresh, dt, max_t)
    
    box = system.copy();
    num_polygons = numel(box.polygons);
    polygons = box.polygons;
    particles_per_side = polygons(1).particles_per_side;

    colors = lines(num_polygons);
    figure;
    subplot(2,1,1);
    hold on;
    axis equal;
    buffer = particles_per_side * polygons(1).sigma;
    xlim([-L/2-buffer, L/2+buffer]);
    ylim([-L/2-buffer, L/2+buffer]);
    
    subplot(2,1,2);
    line_kinetic_trans = animatedline('Color', 'r');  
    line_kinetic_rot = animatedline('Color', 'g');
    line_potential = animatedline('Color', 'b');
    line_total = animatedline('Color', 'k');
    title('Energy'); 
    xlabel('Time');
    ylabel('Energy');
    legend('Kinetic (Trans)', 'Kinetic (Rot)', 'Potential', 'Total');

    max_t = round(max_t / dt) * dt;
    t_values = 0:dt:max_t;    
    num_steps = length(t_values);
    q_series = cell(1, num_polygons);
    v_series = cell(1, num_polygons);
    theta_series = cell(1, num_polygons);
    w_series = cell(1, num_polygons);
    pressure_series = zeros(1, num_steps);
    vertices_series = cell(num_polygons, num_steps);
    xmin_series = -L/2 * ones(1, num_steps);
    xmax_series = L/2 * ones(1, num_steps); 
    ymin_series = -L/2 * ones(1, num_steps);
    ymax_series = L/2 * ones(1, num_steps);
    E_series = zeros(num_steps, 3);
    
    fprintf('Starting energy minimization...\n');
    for step = 1:num_steps
        pressure_series(step) = box.iterate_time(dt);
    
        [E_series(step,1),E_series(step,2)] = box.get_kinetic();
        E_series(step,3) = box.get_potential();
        for polygon = 1:num_polygons
            poly = box.polygons(polygon);
            q_series{polygon}(step, :) = box.apply_pbc2d(poly.q);
            v_series{polygon}(step, :) = poly.v;
            theta_series{polygon}(step) = poly.theta;
            w_series{polygon}(step) = poly.w;
            vertices_series{polygon, step} = poly.vertices_relative + q_series{polygon}(step, :);
        end
  
        if mod(step, 10) == 0
            if mod(step,100) == 0
                fprintf('Energy minimization step: %d / %d, Current energy: %.6f\n', step, num_steps, E_series(step,3));
            end
            % Update the simulation plot
            subplot(2,1,1);
            cla;
            rectangle('Position', [-L/2, -L/2, L, L],...
                'EdgeColor', 'k', 'LineWidth', 2);
            
            for polygon = 1:num_polygons
                sigma = polygons(polygon).sigma;
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
            subplot(2,1,2);
            addpoints(line_kinetic_trans, t_values(step), E_series(step, 1));
            addpoints(line_kinetic_rot, t_values(step), E_series(step, 2));
            addpoints(line_potential, t_values(step), E_series(step, 3));  
            addpoints(line_total, t_values(step), sum(E_series(step,:)));
            set(gca, 'YScale','log');
        end
        
        drawnow; % Flush graphics pipeline at every step
        
        if E_series(step,3) <= E_thresh
            polygons = box.polygons;
            particles_per_side = polygons(1).particles_per_side;
            break; 
        end
    end
    fprintf('Energy minimization complete.\n');
    close all;
end