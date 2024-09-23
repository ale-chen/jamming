%% Display Sim, Energy
% Load the saved simulation data
load('periodic_box_polygon_simulation.mat');
% load('previous_runs/two_triangles_periodic.mat')

% Create a figure with two subplots
fig = figure;
subplot(2, 1, 1);
hold on;
axis equal;
sigma = polygons{1}.sigma;
buffer = particles_per_side * sigma;
xlim([xmin_series(1)-buffer, xmax_series(1)+buffer]);
ylim([ymin_series(1)-buffer, ymax_series(1)+buffer]);
num_polygons = size(polygons,2);
colors = lines(num_polygons);
num_steps = size(pressure_series,2);

subplot(2, 1, 2);
hold on;
xlabel('Time');
ylabel('Energy');

% Initialize energy plot lines
line_kinetic_trans = plot(nan, nan, 'b-', 'LineWidth', 1.5);
line_kinetic_rot = plot(nan, nan, 'r-', 'LineWidth', 1.5);
line_potential = plot(nan, nan, 'g-', 'LineWidth', 1.5);
line_total = plot(nan, nan, 'k-', 'LineWidth', 1.5);
legend('Kinetic (Trans)', 'Kinetic (Rot)', 'Potential', 'Total');

for step = 1:num_steps
    if mod(step,100) == 0
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
%% Display Sim
% % Load the saved simulation data
% load('periodic_box_polygon_simulation.mat');
% 
% fig_sim = figure;
% hold on;
% axis equal;
% buffer = particles_per_side * sigma;
% xlim([xmin-buffer, xmax+buffer]);
% ylim([ymin-buffer, ymax+buffer]);
% colors = lines(num_polygons);
% 
% fig_energy = figure;
% hold on;
% xlabel('Time');
% ylabel('Kinetic Energy');
% 
% kinetic_energy_total = zeros(1, num_steps);
% kinetic_energy_translational = zeros(1, num_steps);
% kinetic_energy_rotational = zeros(1, num_steps);
% line_total = plot(nan, nan, 'b-');
% line_translational = plot(nan, nan, 'r-');
% line_rotational = plot(nan, nan, 'g-');
% 
% legend([line_total, line_translational, line_rotational], ...
%        'Total', 'Translational', 'Rotational');
% 
% for step = 1:num_steps
%     for i = 1:num_polygons
%         mass = sum(box.polygons(i).m);
%         moment_of_inertia = box.polygons(i).mofi;
%         velocity = v_series{i}(step, :);
%         angular_velocity = w_series{i}(step);
% 
%         kinetic_energy_translational(step) = kinetic_energy_translational(step) + ...
%             0.5 * mass * dot(velocity, velocity);
%         kinetic_energy_rotational(step) = kinetic_energy_rotational(step) + ...
%             0.5 * moment_of_inertia * angular_velocity^2;
%     end
%     kinetic_energy_total(step) =...
%         kinetic_energy_translational(step) + kinetic_energy_rotational(step);
% 
%     figure(fig_energy);
%     set(line_total, 'XData', t_values(1:step), 'YData', kinetic_energy_total(1:step));
%     set(line_translational, 'XData', t_values(1:step), 'YData', kinetic_energy_translational(1:step));
%     set(line_rotational, 'XData', t_values(1:step), 'YData', kinetic_energy_rotational(1:step));
% 
%     if mod(step, 10) == 0
%         figure(fig_sim);
%         cla;
%         rectangle('Position', [xmin, ymin, xmax-xmin, ymax-ymin], ...
%                   'EdgeColor', 'k', 'LineWidth', 2);
% 
%         for i = 1:num_polygons
%             vertices = vertices_series{i, step};
% 
%             plot([vertices(:,1); vertices(1,1)], [vertices(:,2); vertices(1,2)], ...
%                  '-', 'Color', colors(i,:), 'LineWidth', 2);
% 
%             for j = 1:size(vertices, 1)
%                 rectangle('Position', [vertices(j,1)-sigma/2, vertices(j,2)-sigma/2, sigma, sigma], ...
%                           'Curvature', [1, 1], 'FaceColor', colors(i,:));
%             end
%         end
% 
%         title(sprintf('Time: %.2f, Pressure: %.2f', t_values(step), pressure_series(step)));
%         drawnow;
%     end
% end
%% Drawing Polygon Test Code
%{
% Set the polygon parameters
sigma = .01;
sides = 3;
k = 1;
particles_per_side = 5;
total_particles = sides * particles_per_side;
m = ones(1,total_particles);  % Mass for each particle
q = [1, 1];  % center of mass position
v = [0, 0];  % initial velocity (2D vector)
w = 0;  % angular velocity (scalar)

% generate the polygon
poly = regular_polygon(sigma, sides, k, m, q, v, w, particles_per_side);

% get the vertices for plotting
vertices = poly.get_vertices();

% plot the polygon edges
figure;
plot([vertices(:,1); vertices(1,1)], [vertices(:,2); vertices(1,2)], 'b-', 'LineWidth', 2);
hold on;

% plot the particles as circles

for i = 1:size(vertices, 1)
    rectangle('Position', [vertices(i,1)-sigma/2, vertices(i,2)-sigma/2, sigma, sigma], ...
        'Curvature', [1, 1], 'FaceColor', 'r');
end

% Set the axis limits and aspect ratio
axis equal;
xlim([-2, 2]);
ylim([-2, 2]);

% Add labels and title
xlabel('X');
ylabel('Y');
title('Regular Polygon with Particles');
%}
%% Old


% % Set up initial conditions
% sigma = 1;
% sides = 3;
% particles_per_side = 1;
% num_polygons = 20;
% polygons = regular_polygon.empty();
% for i = 1:num_polygons
%  q = 10 * rand(1,2); % random position in box
%  v = 0.1 * (rand(1,2) - 0.5); % random initial velocity
%  polygons(i) = regular_polygon(sigma, sides, 500, 0.5, ones(1, sides*particles_per_side), q, v, 0, particles_per_side);
% end
% initial_box_bounds = [[-10 -10]; [10 10]];
% k_int = 100; % interaction spring constant
% box = periodic_box_polygon(polygons, initial_box_bounds, k_int);
% dt = 0.001;
% t_end = 10;
% pressure_threshold = 1e3;
% pressure_var_threshold = 1e1;
% pressures = [];
% % Run simulation
% figure;
% while box.t < t_end
%  box
%  box.iterate_verlet(dt);
%  pressure = box.get_pressure();
%  pressures = [pressures pressure];
% % contract box size
% if mod(round(box.t / dt), 100) == 0
%  box.box_bounds = box.box_bounds * 0.99;
% end
% % check jamming condition
% if length(pressures) > 100
% if mean(pressures(end-99:end)) > pressure_threshold && var(pressures(end-99:end)) < pressure_var_threshold
%  disp('Jammed packing achieved');
% break;
% end
% end
% % Plot polygons
%  cla;
%  hold on;
% for i = 1:length(box.polygons)
%  vertices = box.polygons(i).vertices;
%  plot(vertices(:,1), vertices(:,2), 'b-', 'LineWidth', 2);
% % Plot polygon centers
%  plot(box.polygons(i).q(1), box.polygons(i).q(2), 'ro', 'MarkerSize', 10);
% % Plot particles
% for j = 1:size(vertices, 1)
%  rectangle('Position', [vertices(j,1)-sigma/2, vertices(j,2)-sigma/2, sigma, sigma], ...
% 'Curvature', [1, 1], 'FaceColor', 'r');
% end
% end
%  xlim(box.box_bounds(:,1)');
%  ylim(box.box_bounds(:,2)');
%  xlabel('X');
%  ylabel('Y');
%  title(sprintf('Polygons in Periodic Box, t = %.3f', box.t));
%  drawnow;
% % Print the current time step
%  fprintf('Time step: %.3f\n', box.t);
% end