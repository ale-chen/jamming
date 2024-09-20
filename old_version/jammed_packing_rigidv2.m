% Main simulation code
% Simulation parameters
num_bodies = 10;
total_time = 1;
dt = 0.0001;
k_inter = 1000;
sigma_micro = 1;
k_damping = 0.9;  % Damping coefficient
pressure_tolerance = 0.01;  % Tolerance for pressure convergence

% Create a single rigid body shape to be repeated
fprintf('Drawing the Rigid Body shape to be repeated:\n');
template_rigid_body = create_rigid_body();
close;  % Close the figure after the rigid body is drawn
fprintf('Template Rigid Body created successfully.\n');

% Initialize rigid bodies with randomized initial positions
rigid_bodies = cell(num_bodies, 1);
for i = 1:num_bodies
    rigid_bodies{i} = copy_rigid_body(template_rigid_body);
    fprintf('Rigid Body %d created successfully.\n', i);
end

% Ensure non-overlapping initial positions
rigid_bodies = initialize_positions(rigid_bodies, sigma_micro);

% Initialize velocities
for i = 1:num_bodies
    rigid_bodies{i}.velocity_history = zeros(1, 2);
    rigid_bodies{i}.angular_velocity_history = 0;
end

% Binary search for jammed packing
box_size_low = 0;
box_size_high = 100;
target_pressure = 1e-3;  % Desired pressure for jammed packing
max_iterations = 100;  % Maximum number of iterations for binary search

% Initialize figure for in-time visualization
fig = figure;
hold on;
axis equal;
title('Jammed Packing Simulation');

for iteration = 1:max_iterations
    box_size = (box_size_low + box_size_high) / 2;

    % Reset time and collision count
    t = 0;
    pressure_history = [];

    % Simulation loop
    while t < total_time
        % Update positions and orientations of rigid bodies
        for i = 1:num_bodies
            rigid_bodies{i} = update_rigid_body(rigid_bodies{i}, dt);
        end

        % Calculate forces and torques on each rigid body
        for i = 1:num_bodies
            [force_total, torque_total] = calculate_total_force_torque(rigid_bodies, i, k_inter, sigma_micro, box_size);
            rigid_bodies{i} = apply_force_torque(rigid_bodies{i}, force_total, torque_total, k_damping, dt);
        end

        % Calculate instantaneous pressure
        pressure = calculate_pressure(rigid_bodies, box_size, k_inter, sigma_micro);
        pressure_history = [pressure_history; pressure];

        % Visualize current state
        if mod(round(t/dt), 100) == 0  % Update visualization every 100 steps
            visualize_system(fig, rigid_bodies, box_size);
            drawnow;
        end

        % Increment time
        t = t + dt;
    end

    % Calculate average pressure over the simulation
    avg_pressure = mean(pressure_history);

    % Check if pressure is within tolerance
    if abs(avg_pressure - target_pressure) < pressure_tolerance
        fprintf('Jammed packing found with box size: %.4f\n', box_size);
        break;
    elseif avg_pressure < target_pressure
        box_size_high = box_size;  % Decrease box size
    else
        box_size_low = box_size;  % Increase box size
    end

    fprintf('Iteration %d: box_size = %.4f, avg_pressure = %.6f\n', iteration, box_size, avg_pressure);
end

if iteration == max_iterations
    fprintf('Warning: Maximum iterations reached. Jammed packing may not be found.\n');
end

% Final visualization
visualize_system(fig, rigid_bodies, box_size);
title('Final Jammed Packing');

% Helper functions

function rigid_body = create_rigid_body()
    disp('Please draw a polygon.');
    disp('Left-click to add points, right-click or press Enter to finish drawing.');
    figure;
    hold on;
    axis equal;
    xlim([-10 10]);
    ylim([-10 10]);
    sigma = 2;

    x = [];
    y = [];
    while true
        [xi, yi, button] = ginput(1);
        if button == 1  % Left-click
            x = [x; xi];
            y = [y; yi];
            plot(x, y, 'b-', 'LineWidth', 2);
        else  % Right-click, Enter key, or any other button
            break;
        end
    end

    [x_opt, y_opt] = optimize_vertices(x, y, sigma);
    particle_centers = generate_particle_centers(x_opt, y_opt, sigma);

    mass_vector = 0.5 + rand(size(particle_centers, 1), 1);
    rigid_body = RigidBody(particle_centers, mass_vector);
end

function new_rigid_body = copy_rigid_body(template_rigid_body)
    new_rigid_body = RigidBody(template_rigid_body.particle_centers, template_rigid_body.mass_vector);
end

function rigid_body = update_rigid_body(rigid_body, dt)
    % Update position and orientation
    rigid_body.position_history(end+1, :) = rigid_body.position_history(end, :) + rigid_body.velocity_history(end, :) * dt;
    rigid_body.orientation_history(end+1) = rigid_body.orientation_history(end) + rigid_body.angular_velocity_history(end) * dt;
end

function coords = get_particle_coordinates(rigid_body, frame_index)
    rotation_matrix = [cos(rigid_body.orientation_history(frame_index)), -sin(rigid_body.orientation_history(frame_index)); 
                       sin(rigid_body.orientation_history(frame_index)), cos(rigid_body.orientation_history(frame_index))];
    relative_positions = rigid_body.particle_centers - repmat(rigid_body.center_of_mass, size(rigid_body.particle_centers, 1), 1);
    rotated_positions = (rotation_matrix * relative_positions')';
    coords = rotated_positions + repmat(rigid_body.position_history(frame_index, :), size(rigid_body.particle_centers, 1), 1);
end

function moment_of_inertia = calculate_moment_of_inertia(centers, mass_vector)
    moment_of_inertia = sum(mass_vector .* sum(centers.^2, 2));
end

function [force_total, torque_total, collision_count] = calculate_total_force_torque(rigid_bodies, body_index, k_inter, sigma_micro, box_size)
    num_bodies = length(rigid_bodies);
    force_total = zeros(2, 1);
    torque_total = 0;
    collision_count = 0;

    % Interaction forces between particles of different rigid bodies
    for j = 1:num_bodies
        if body_index ~= j
            coords_i = get_particle_coordinates(rigid_bodies{body_index}, size(rigid_bodies{body_index}.position_history, 1));
            coords_j = get_particle_coordinates(rigid_bodies{j}, size(rigid_bodies{j}.position_history, 1));

            for k = 1:size(coords_i, 1)
                for l = 1:size(coords_j, 1)
                    diff = coords_i(k, :)' - coords_j(l, :)';
                    distance = norm(diff);

                    if distance <= sigma_micro
                        direction = diff / distance;
                        force = k_inter * direction;

                        force_total = force_total + force;
                        r = coords_i(k, :) - rigid_bodies{body_index}.center_of_mass;
                        torque = cross([r, 0], [force', 0]);  % Ensure 3D vectors for cross product
                        torque_total = torque_total + torque(3);
                    end
                end
            end
        end
    end
    
    % Collision forces with walls
    coords_i = get_particle_coordinates(rigid_bodies{body_index}, size(rigid_bodies{body_index}.position_history, 1));
    for k = 1:size(coords_i, 1)
        penetration_left = -box_size/2 - coords_i(k, 1);
        penetration_right = coords_i(k, 1) - box_size/2;
        penetration_bottom = -box_size/2 - coords_i(k, 2);
        penetration_top = coords_i(k, 2) - box_size/2;

        if penetration_left > 0
            force_magnitude = k_inter * penetration_left;
            force_total(1) = force_total(1) + force_magnitude;
            r = coords_i(k, :) - rigid_bodies{body_index}.center_of_mass;
            torque = cross([r, 0], [force_magnitude, 0, 0]);
            torque_total = torque_total + torque(3);
            collision_count = collision_count + 1;
        elseif penetration_right > 0
            force_magnitude = k_inter * penetration_right;
            force_total(1) = force_total(1) - force_magnitude;
            r = coords_i(k, :) - rigid_bodies{body_index}.center_of_mass;
            torque = cross([r, 0], [-force_magnitude, 0, 0]);
            torque_total = torque_total + torque(3);
            collision_count = collision_count + 1;
        end

        if penetration_bottom > 0
            force_magnitude = k_inter * penetration_bottom;
            force_total(2) = force_total(2) + force_magnitude;
            r = coords_i(k, :) - rigid_bodies{body_index}.center_of_mass;
            torque = cross([r, 0], [0, force_magnitude, 0]);
            torque_total = torque_total + torque(3);
            collision_count = collision_count + 1;
        elseif penetration_top > 0
            force_magnitude = k_inter * penetration_top;
            force_total(2) = force_total(2) - force_magnitude;
            r = coords_i(k, :) - rigid_bodies{body_index}.center_of_mass;
            torque = cross([r, 0], [0, -force_magnitude, 0]);
            torque_total = torque_total + torque(3);
            collision_count = collision_count + 1;
        end
    end
end

function rigid_body = apply_force_torque(rigid_body, force_total, torque_total, k_damping, dt)
    % Apply total force and torque to the rigid body
    rigid_body.velocity_history(end+1, :) = rigid_body.velocity_history(end, :) + (force_total' / sum(rigid_body.mass_vector)) - k_damping * rigid_body.velocity_history(end, :) * dt;
    rigid_body.angular_velocity_history(end+1) = rigid_body.angular_velocity_history(end) + torque_total / rigid_body.moment_of_inertia - k_damping * rigid_body.angular_velocity_history(end) * dt;
end

function pressure = calculate_pressure(rigid_bodies, box_size, k_inter, sigma_micro)
    num_bodies = length(rigid_bodies);
    total_force = 0;
    
    for i = 1:num_bodies
        coords_i = get_particle_coordinates(rigid_bodies{i}, size(rigid_bodies{i}.position_history, 1));
        
        % Wall forces
        for k = 1:size(coords_i, 1)
            penetration_left = -box_size/2 - coords_i(k, 1);
            penetration_right = coords_i(k, 1) - box_size/2;
            penetration_bottom = -box_size/2 - coords_i(k, 2);
            penetration_top = coords_i(k, 2) - box_size/2;
            
            if penetration_left > 0
                total_force = total_force + k_inter * penetration_left;
            end
            if penetration_right > 0
                total_force = total_force + k_inter * penetration_right;
            end
            if penetration_bottom > 0
                total_force = total_force + k_inter * penetration_bottom;
            end
            if penetration_top > 0
                total_force = total_force + k_inter * penetration_top;
            end
        end
        
        % Particle-particle interactions
        for j = i+1:num_bodies
            coords_j = get_particle_coordinates(rigid_bodies{j}, size(rigid_bodies{j}.position_history, 1));
            
            for k = 1:size(coords_i, 1)
                for l = 1:size(coords_j, 1)
                    diff = coords_i(k, :)' - coords_j(l, :)';
                    distance = norm(diff);
                    
                    if distance <= sigma_micro
                        total_force = total_force + k_inter * (sigma_micro - distance);
                    end
                end
            end
        end
    end
    
    % Calculate pressure (force per unit length in 2D)
    pressure = total_force / (4 * box_size);
end

function visualize_system(fig, rigid_bodies, box_size)
    clf(fig);
    hold on;
    axis equal;
    xlim([-box_size/2, box_size/2]);
    ylim([-box_size/2, box_size/2]);
    
    for i = 1:length(rigid_bodies)
        coords = get_particle_coordinates(rigid_bodies{i}, size(rigid_bodies{i}.position_history, 1));
        plot(coords(:, 1), coords(:, 2), 'o');
    end
    
    % Draw box
    rectangle('Position', [-box_size/2, -box_size/2, box_size, box_size], 'EdgeColor', 'r', 'LineWidth', 2);
end