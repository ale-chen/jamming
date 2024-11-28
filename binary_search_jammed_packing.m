% Parameters
num_polygons = 10;
sigma = 1;
sides = 3;
particles_per_side = 2;
total_particles = sides * particles_per_side;
m = ones(1, total_particles);
k_int = 10000;
b_trans = 100;
b_ang = 100;
dt = 0.001;
max_t = 20;

% Binary search parameters
L = 7.9;
L_old = L;
r_scale = 0.999; 
P_t = 10; 
P_l = 0.95 * P_t;  
P_h = 1.05 * P_t;
L_l = -1;
L_h = -1;
L_tol = 0.001;
E_thresh = 1e-4;
num_minimization_steps = 5000;

% Initialize system
polygons = generate_initial_polygons(num_polygons, sigma, sides, particles_per_side, L, m);
box = periodic_box_polygon(polygons, k_int, b_trans, b_ang, -L/2, L/2, -L/2, L/2);
old_state = box.copy();

% Create data directory if it doesn't exist
if ~exist('./data', 'dir')
    mkdir('./data');
end

% Initialize visualization
fig = figure('Position', [100, 100, 800, 1000]);
subplot(2,1,1);
hold on;
axis equal;
buffer = particles_per_side * sigma;
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
set(gca, 'YScale', 'log');

colors = lines(num_polygons);
total_t = 0;

% Binary search for jammed state
fprintf('Starting binary search...\n');
iteration = 0;
success = false;


while true
    iteration = iteration + 1;
    fprintf('Binary search iteration: %d\n', iteration);
    
    % Energy minimization
    [box, converged_energy] = energy_minimize(box, E_thresh, dt, num_minimization_steps, colors, ...
        line_kinetic_trans, line_kinetic_rot, line_potential, line_total);
    fprintf('Energy after minimization: %.6e\n', converged_energy);
    
    % Calculate pressure after minimization
    [P, overlap_count] = calculate_pressure(box);
    fprintf('Current L: %.6f, Current P: %.6f (Target: %.6f), Contacts: %d\n', L, P, P_t, overlap_count);
    
    % Binary search logic with explicit state tracking
    if P < P_l
        fprintf('  P < P_l: System is UNJAMMED at L = %.6f\n', L);
        old_state = box.copy();  % Save unjammed state
        L_old = L;              % Save current length for rescaling
        L_h = L;                % This length becomes upper bound (unjammed)
        
        if L_l > 0
            fprintf('    Known jammed state exists at L_l = %.6f\n', L_l);
            L = (L_h + L_l)/2;  % Try halfway
            fprintf('    Trying halfway point L = %.6f\n', L);
            L_l = -1;           % Reset L_l as rearrangement possible
        else
            fprintf('    No known jammed state, shrinking by r_scale\n');
            L = L * r_scale;    % Shrink by fixed amount
        end
    else
        if P > P_h
            fprintf('  P > P_h: System is JAMMED at L = %.6f\n', L);
            L_l = L;                % This length becomes lower bound (jammed)
            box = old_state.copy();  % Reset to last unjammed state
            L = (L_h + L_l)/2;      % Try halfway between jammed and unjammed
            fprintf('    Reset to unjammed state at L_h = %.6f\n', L_h);
            fprintf('    Trying halfway point L = %.6f\n', L);
        else
            fprintf('  P_l ≤ P ≤ P_h: Target pressure achieved!\n');
            success = true;
            break;
        end
    end
    
    % Update box size and rescale system
    scale_factor = L/L_old;
    fprintf('  Rescaling system by factor: %.6f\n', scale_factor);
    
    for i = 1:numel(box.polygons)
        box.polygons(i).q = box.polygons(i).q * scale_factor;
        box.polygons(i).v = box.polygons(i).v * sqrt(scale_factor);
    end
    box.xmin = -L/2;
    box.xmax = L/2;
    box.ymin = -L/2;
    box.ymax = L/2;
    
    % Check for convergence failure with detailed state info
    if L_l > 0 && L_h > 0 && abs((L_h / L_l) - 1) < L_tol
        fprintf('\nBinary search failed to find jammed state:\n');
        fprintf('L_h = %.6f (unjammed, P < %.6f)\n', L_h, P_l);
        fprintf('L_l = %.6f (jammed, P > %.6f)\n', L_l, P_h);
        fprintf('(L_h/L_l - 1) = %.6e (tolerance: %.6e)\n', abs(L_h/L_l - 1), L_tol);
        error('No jammed packing found in this pressure window.');
    end
end
if success
    fprintf('Jammed state found! Saving final configuration...\n');
    save_jammed_state(box, L, P);
end

% Support functions
function [box, final_energy] = energy_minimize(box, E_thresh, dt, max_steps, colors, ...
    line_kinetic_trans, line_kinetic_rot, line_potential, line_total)
    fprintf('Starting energy minimization...\n');
    
    last_energy = inf;
    t = 0;
    
    % Get handle to total_t from base workspace
    total_t = evalin('base', 'total_t');
    
    for step = 1:max_steps
        % Perform time step
        box.iterate_time(dt);
        t = t + dt;
        
        % Calculate energies and pressure
        [K_trans, K_rot] = box.get_kinetic();
        V = box.get_potential();
        current_energy = K_trans + K_rot + V;
        [virial_pressure, overlap_count] = calculate_pressure(box);
        
        % Update visualization every 10 steps
        if mod(step, 10) == 0
            % Clear and update system plot
            subplot(2,1,1);
            cla;
            
            % Draw simulation box
            L = box.xmax - box.xmin;
            rectangle('Position', [box.xmin, box.ymin, L, L], 'EdgeColor', 'k', 'LineWidth', 1);
            
            % Draw polygons
            for i = 1:numel(box.polygons)
                poly = box.polygons(i);
                % Wrap q using PBC
                wrapped_q = box.pbc_vector_difference(poly.q, [0, 0]);
                vertices = poly.vertices_relative + wrapped_q;
                
                % Draw edges
                plot([vertices(:,1); vertices(1,1)], [vertices(:,2); vertices(1,2)], ...
                    '-', 'Color', colors(i,:), 'LineWidth', 2);
                
                % Draw vertices as circles
                for j = 1:size(vertices, 1)
                    rectangle('Position', [vertices(j,1)-poly.sigma/2, vertices(j,2)-poly.sigma/2, ...
                        poly.sigma, poly.sigma], 'Curvature', [1, 1], ...
                        'FaceColor', colors(i,:), 'EdgeColor', 'none');
                end
                
                % Draw center
                plot(wrapped_q(1), wrapped_q(2), 'k.', 'MarkerSize', 10);
            end
            
            title(sprintf('t = %.2f, Virial P = %.2f, Contacts = %d', total_t + t, virial_pressure, overlap_count));
            
            % Update energy plot with continuous time
            subplot(2,1,2);
            addpoints(line_kinetic_trans, total_t + t, K_trans);
            addpoints(line_kinetic_rot, total_t + t, K_rot);
            addpoints(line_potential, total_t + t, V);
            addpoints(line_total, total_t + t, current_energy);
            
            drawnow;
        end
        
        % Check for convergence
        if abs(current_energy - last_energy) < E_thresh
            fprintf('Converged at step %d with energy %.6e\n', step, current_energy);
            break;
        end
        
        % Progress update
        if mod(step, 500) == 0
            fprintf('Step %d: Energy = %.6e\n', step, current_energy);
        end
        
        last_energy = current_energy;
    end
    
    % Update total_t in base workspace
    assignin('base', 'total_t', total_t + t);
    
    final_energy = current_energy;
end

function [pressure, overlap_count] = calculate_pressure(box)
    % Calculate virial pressure using center-to-center distances
    virial = 0;
    overlap_count = 0;
    
    % Sum over all polygon pairs
    for i = 1:numel(box.polygons)-1
        for j = i+1:numel(box.polygons)
            poly_i = box.polygons(i);
            poly_j = box.polygons(j);
            
            % Get center-to-center distance
            Rij = box.pbc_vector_difference(poly_i.q, poly_j.q);
            
            vertices_i = poly_i.get_vertices();
            vertices_j = poly_j.get_vertices();
            
            % Sum over all vertex pairs
            for ii = 1:size(vertices_i, 1)
                for jj = 1:size(vertices_j, 1)
                    rij = box.pbc_vector_difference(vertices_i(ii, :), vertices_j(jj, :));
                    distance = norm(rij);
                    
                    sigma = (poly_i.sigma + poly_j.sigma) / 2;
                    
                    if distance < sigma
                        overlap_count = overlap_count + 1;
                        force = box.k_int * (sigma - distance) * rij / distance;
                        % Use Rij (center-to-center) in virial calculation
                        virial = virial + dot(Rij, force);
                    end
                end
            end
        end
    end
    
    pressure = virial / (2 * box.get_area);
end

function save_jammed_state(box, L, P)
    % Create data structure for saving
    data = struct();
    data.L = L;
    data.P = P;
    data.num_polygons = numel(box.polygons);
    data.positions = zeros(data.num_polygons, 2);
    data.vertices = cell(data.num_polygons, 1);
    data.velocities = zeros(data.num_polygons, 2);
    data.angles = zeros(data.num_polygons, 1);
    data.angular_velocities = zeros(data.num_polygons, 1);
    
    % Store polygon data with PBC wrapping
    for i = 1:data.num_polygons
        % Wrap positions using PBC
        data.positions(i,:) = box.pbc_vector_difference(box.polygons(i).q, [0,0]);
        data.velocities(i,:) = box.polygons(i).v;
        data.angles(i) = box.polygons(i).theta;
        data.angular_velocities(i) = box.polygons(i).w;
        data.vertices{i} = box.polygons(i).vertices_relative + data.positions(i,:);
    end
    
    % Save data
    save('./data/jammed_state.mat', 'data');
    
    % Create final visualization
    figure('Position', [100, 100, 800, 800]);
    hold on;
    axis equal;
    
    % Draw polygons
    colors = lines(data.num_polygons);
    buffer = box.polygons(1).sigma * box.polygons(1).particles_per_side;
    
    for i = 1:data.num_polygons
        vertices = data.vertices{i};
        
        % Draw edges
        plot([vertices(:,1); vertices(1,1)], [vertices(:,2); vertices(1,2)], ...
            '-', 'Color', colors(i,:), 'LineWidth', 2);
        
        % Draw center
        plot(data.positions(i,1), data.positions(i,2), 'k.', 'MarkerSize', 10);
        
        % Draw vertices as circles
        for j = 1:size(vertices, 1)
            rectangle('Position', [vertices(j,1)-box.polygons(1).sigma/2, ...
                vertices(j,2)-box.polygons(1).sigma/2, ...
                box.polygons(1).sigma, box.polygons(1).sigma], ...
                'Curvature', [1, 1], 'FaceColor', colors(i,:), 'EdgeColor', 'none');
        end
    end
    
    title(sprintf('Jammed State (L = %.3f, P = %.3f)', data.L, data.P));
    xlabel('x');
    ylabel('y');
    axis([-L/2-buffer L/2+buffer -L/2-buffer L/2+buffer]);
    
    % Save figures
    saveas(gcf, './data/jammed_state.fig');
    saveas(gcf, './data/jammed_state.png');
    close(gcf);
end

function polygons = generate_initial_polygons(num_polygons, sigma, sides, particles_per_side, L, m)
    polygons = cell(1, num_polygons);
    
    grid_size = ceil(sqrt(num_polygons));
    polygon_size = sigma * particles_per_side;
    spacing = (L - grid_size * polygon_size) / (grid_size + 1);
    
    for i = 1:num_polygons
        row = floor((i - 1) / grid_size);
        col = mod(i - 1, grid_size);
        q = [-L/2 + (col + 1) * spacing + col * polygon_size, ...
             -L/2 + (row + 1) * spacing + row * polygon_size];
        v = rand(1,2) - 0.5;
        w = 0;
        polygons{i} = regular_polygon(sigma, sides, m, q, v, w, particles_per_side);
    end
end