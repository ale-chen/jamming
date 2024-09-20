classdef RigidBody
    properties
        particle_centers
        mass_vector
        center_of_mass
        moment_of_inertia
        position_history
        orientation_history
        velocity_history
        angular_velocity_history
    end
    
    methods
        function obj = RigidBody(particle_centers, mass_vector)
            obj.particle_centers = particle_centers;
            obj.mass_vector = mass_vector;
            obj.center_of_mass = sum(particle_centers .* mass_vector, 1) / sum(mass_vector);
            obj.moment_of_inertia = calculate_moment_of_inertia(particle_centers - repmat(obj.center_of_mass, size(particle_centers, 1), 1), mass_vector);
            obj.position_history = zeros(1, 2);
            obj.orientation_history = 0;
            obj.velocity_history = zeros(1, 2);
            obj.angular_velocity_history = 0;
        end
        
        function update_position(obj, dt)
            obj.position_history(end+1, :) = obj.position_history(end, :) + obj.velocity_history(end, :) * dt;
            obj.orientation_history(end+1) = obj.orientation_history(end) + obj.angular_velocity_history(end) * dt;
        end
        
        function apply_force(obj, force, application_point)
            obj.velocity_history(end+1, :) = obj.velocity_history(end, :) + force / sum(obj.mass_vector);
            r = application_point - obj.center_of_mass;
            torque = cross([r, 0], [force, 0]);  % Ensure 3D vectors for cross product
            obj.angular_velocity_history(end+1) = obj.angular_velocity_history(end) + torque(3) / obj.moment_of_inertia;
        end
        
        function coords = get_particle_coordinates(obj, index)
            rotation_matrix = [cos(obj.orientation_history(index)), -sin(obj.orientation_history(index)); 
                                sin(obj.orientation_history(index)), cos(obj.orientation_history(index))];
            relative_positions = obj.particle_centers - repmat(obj.center_of_mass, size(obj.particle_centers, 1), 1);
            rotated_positions = (rotation_matrix * relative_positions')';
            coords = rotated_positions + repmat(obj.position_history(index, :), size(obj.particle_centers, 1), 1);
        end
    end
end

function moment_of_inertia = calculate_moment_of_inertia(centers, mass_vector)
    moment_of_inertia = sum(mass_vector .* sum(centers.^2, 2));
end