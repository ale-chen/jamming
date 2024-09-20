% Main simulation code
% Simulation parameters
num_bodies = 1;
total_time = 1;
dt = 0.0001;
k_inter = 1000;
sigma_micro = 1;
k_damping = 0.9;  % Damping coefficient
pressure_tolerance = 0.01;  % Tolerance for pressure convergence

% Initialize rigid bodies with randomized initial positions
rigid_bodies = cell(num_bodies, 1);
for i = 1:num_bodies
    fprintf('Drawing Rigid Body %d of %d:\n', i, num_bodies);
    rigid_bodies{i} = create_rigid_body();
    close;  % Close the figure after each rigid body is drawn
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



function rigid_bodies = initialize_positions(rigid_bodies, sigma_micro)
    num_bodies = length(rigid_bodies);
    min_distance = 10 * sigma_micro;  % Adjust the minimum distance based on macro particle size
    
    for i = 1:num_bodies
        while true
            % Generate random initial position
            pos = (rand(2, 1) - 0.5) * 100;
            rigid_bodies{i}.position_history(1, :) = pos';
            
            % Check for overlap with other macro particles
            overlap = false;
            for j = 1:i-1
                dist = norm(pos - rigid_bodies{j}.position_history(1, :)');
                if dist < min_distance
                    overlap = true;
                    break;
                end
            end
            
            if ~overlap
                break;
            end
        end
    end
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