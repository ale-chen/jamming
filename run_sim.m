% % parameters
% num_polygons = 15;
% sigma = 1;
% sides = 3;
% particles_per_side = 2;
% total_particles = sides * particles_per_side;
% m = ones(1, total_particles);
% k_int = 10000;
% b_trans = 0;
% b_ang = 0;
% xmin = -10;
% xmax = 10;
% ymin = -10;
% ymax = 10;
% dt = 0.001;
% max_t = 20;
% box_shrink_rate = 0.00; % per second
% 
% polygons = cell(1, num_polygons);
% 
% % generate initial conditions without overlap
% pressure = 1;
% while pressure > 0
%     for i = 1:num_polygons
%         q = [rand()*(xmax-xmin)+xmin, rand()*(ymax-ymin)+ymin];
%         v = rand(1, 2)*7 - 3.5;
%         w = 0;
%         polygons{i} = regular_polygon(sigma, sides, m, q, v, w, particles_per_side);
%     end
%     box = periodic_box_polygon(polygons, k_int, b_trans, b_ang, xmin, xmax, ymin, ymax);
%     pressure = box.iterate_time(dt); % calculate initial pressure, is it
%     % ok to just do this? I think you're adding a time step in the final
%     % run.
% end
% 
% disp('Found Starting Config')

% SIMPLE SETUP FOR DEBUGGING
tiny = 0.0000000001;

num_polygons = 2;
sigma = 3;
sides = 3;
particles_per_side = 1;
total_particles = sides * particles_per_side;
m = [1,1,1]; % some particles have no mass, cofm is in first particle
k_int = 100;
b_trans = 0;
b_ang = 0;
xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
dt = 0.001;
max_t = 20;

box_shrink_rate = 0; % per second

polygons = cell(1, num_polygons);

polygons{1} = regular_polygon(sigma, sides, m, [-6,0], [2,0], 0, particles_per_side);
polygons{2} = regular_polygon(sigma, sides, m, [0,0], [0,0], 0, particles_per_side);

box = periodic_box_polygon(polygons, k_int, b_trans, b_ang, xmin, xmax, ymin, ymax);
box.polygons(2).rotate_vertices(pi);

% store the data time series
t_values = 0:dt:max_t;
num_steps = length(t_values);
q_series = cell(1, num_polygons);
v_series = cell(1, num_polygons);
theta_series = cell(1, num_polygons);
w_series = cell(1, num_polygons);
pressure_series = zeros(1, num_steps);
vertices_series = cell(num_polygons, num_steps);
xmin_series = zeros(1, num_steps);
xmax_series = zeros(1, num_steps);
ymin_series = zeros(1, num_steps);
ymax_series = zeros(1, num_steps);
E_series = zeros(num_steps, 3);

max_t = round(max_t / dt) * dt;
progress_step = floor(num_steps / 100); % Update progress every 1% of total steps
box_shrink_rate_dt = box_shrink_rate * max_t / num_steps;

for step = 1:num_steps
    [pressure_series(step),...
    E_series(step,1),E_series(step,2),E_series(step,3)]...
        = box.iterate_time(dt);

    xmin_series(step) = box.xmin;
    xmax_series(step) = box.xmax;
    ymin_series(step) = box.ymin;
    ymax_series(step) = box.ymax;

    for polygon = 1:num_polygons
        [poly] = deal(box.polygons(polygon));
        q_series{polygon}(step, :) = box.apply_pbc2d(poly.q);
        v_series{polygon}(step, :) = poly.v;
        theta_series{polygon}(step) = poly.theta;
        w_series{polygon}(step) = poly.w;
        vertices_series{polygon, step} = poly.vertices_relative + q_series{polygon}(step, :);
    end

    % shrink box
    box.xmin = box.xmin * (2 - box_shrink_rate_dt) / 2;
    box.xmax = box.xmax * (2 - box_shrink_rate_dt) / 2;
    box.ymin = box.ymin * (2 - box_shrink_rate_dt) / 2;
    box.ymax = box.ymax * (2 - box_shrink_rate_dt) / 2;

    if mod(step, progress_step) == 0
        progress = step / num_steps * 100;
        fprintf('Sim Progress: [%-50s] %.1f%%\r', repmat('#', 1, floor(progress/2)), progress);
    end
end
% complete the progress bar
fprintf('Progress: [%-50s] 100.0%%\n', repmat('#', 1, 50));

disp('Finished Simulation')
save('periodic_box_polygon_simulation.mat','t_values','E_series', 'q_series', 'v_series',...
    'theta_series', 'w_series', 'pressure_series', 'vertices_series',...
    'box', 'polygons','particles_per_side', 'xmin_series', 'xmax_series', 'ymin_series', 'ymax_series');

figure;
hold on;
axis equal;
buffer = particles_per_side * sigma;
xlim([xmin-buffer, xmax+buffer]);
ylim([ymin-buffer, ymax+buffer]);

colors = lines(num_polygons);

for step = 1:num_steps
    if mod(step,100) == 0
        cla;
        rectangle('Position', [xmin_series(step), ymin_series(step),...
            xmax_series(step)-xmin_series(step),...
            ymax_series(step)-ymin_series(step)],...
            'EdgeColor', 'k', 'LineWidth', 2);

        for polygon = 1:num_polygons
            vertices = vertices_series{polygon, step};
            
            % polygon edges
            plot([vertices(:,1); vertices(1,1)], [vertices(:,2); vertices(1,2)], '-', 'Color', colors(polygon,:), 'LineWidth', 2);
    
            for j = 1:size(vertices, 1)
                rectangle('Position', [vertices(j,1)-sigma/2, vertices(j,2)-sigma/2, sigma, sigma], ...
                    'Curvature', [1, 1], 'FaceColor', colors(polygon,:));
            end
        end
    end
    title(sprintf('Time: %.2f, Pressure: %.2f', t_values(step), pressure_series(step)));

    drawnow()
end

% % params
% num_polygons = 5;
% sigma = 1;
% sides = 4;
% particles_per_side = 1;
% k = 1;
% total_particles = sides * particles_per_side;
% m = ones(1,total_particles);
% v = [0,0];
% w = 0;
% k_int = 1;
% b_trans = 0.1;
% b_ang = 0.1;
% xmin_init = -10;
% xmax_init = 10;
% ymin_init = -10;
% ymax_init = 10;
% dt = 0.01;
% max_t = 10;
% shrink_factor = 0.9;
% 
% max_t = round(max_t / dt) * dt;
% 
% polygons = cell(1, num_polygons);
% 
% while true
%     for i = 1:num_polygons
%         q = [rand()*(xmax_init-xmin_init)+xmin_init, rand()*(ymax_init-ymin_init)+ymin_init];
%         polygons{i} = regular_polygon(sigma, sides, k, m, q, v, w, particles_per_side);
%     end
% 
%     box = periodic_box_polygon(polygons, k_int, b_trans, b_ang, xmin_init, xmax_init, ymin_init, ymax_init);
%     pressure = box.iterate_time(dt); % Calculate initial pressure
% 
%     if pressure == 0
%         break; % no overlap, exit the loop
%     end
% end
% disp('Starting Value Found')
% 
% t_values = 0:dt:max_t;
% num_steps = length(t_values);
% q_series = cell(1, num_polygons);
% v_series = cell(1, num_polygons);
% theta_series = cell(1, num_polygons);
% w_series = cell(1, num_polygons);
% pressure_series = zeros(1, num_steps);
% 
% % Run the simulation
% for step = 1:num_steps
%     pressure_series(step) = box.iterate_time(dt);
% 
%     for i = 1:num_polygons
%         poly = box.polygons(i);
%         q_series{i}(step, :) = poly.q;
%         v_series{i}(step, :) = poly.v;
%         theta_series{i}(step) = poly.theta;
%         w_series{i}(step) = poly.w;
%     end
% 
%     % shrink
%     box.xmin = box.xmin * shrink_factor;
%     box.xmax = box.xmax * shrink_factor;
%     box.ymin = box.ymin * shrink_factor;
%     box.ymax = box.ymax * shrink_factor;
% end
% 
% % Save the data time series and polygon parameters
% save('periodic_box_polygon_simulation.mat', 'q_series', 'v_series', 'theta_series', 'w_series', 'pressure_series', 'box', 'polygons');