
load('jammed_packings/iteration48L8.341P123.4472.mat');

jammed_box = box;

function V = total_potential_energy(box)
    V = box.get_potential();
end

% Extract the box object and relevant parameters
box = jammed_box;
dt = 0.001; % Adjust the time step as needed
max_t = 20; % Adjust the maximum time as needed
E_thresh = 1e-6; % Adjust the energy threshold as needed

% Energy minimization parameters
min_thresh = 1e-16; % Minimum potential energy per particle threshold
max_steps = 10000; % Maximum number of minimization steps

% Energy minimization loop
for step = 1:max_steps
    % Calculate current potential energy per particle
    V = box.get_potential();
    disp(box.get_potential());
    V_per_particle = V / numel(box.polygons);

    % Check if the potential energy is below the threshold
    if V_per_particle < min_thresh
        fprintf('Energy minimization converged after %d steps.\n', step);
        break;
    end

    % Perform one step of the relaxation routine
    box.iterate_time(dt);
end

% Update the jammed_box variable with the minimized box
jammed_box = box;

dt = 0.001;
steps = 1000;
epsilon = 0.05;

[dUdV, strained_box] = apply_isometric_strain(jammed_box, dt, steps, epsilon);

disp(dUdV);

function [dUdV, strained_box] = apply_isometric_strain(box, dt, steps, epsilon)
    % Apply an isometric strain. x_i = x_i * (1 + \epsilon),
    % y_i = y_i *(1+ \epsilon), L_x = L_x *(1+\epsilon), L_x = L_x *(1+\epsilon).
    % Measure U and V before and after the strain. 

    strained_box = box.copy();
    
    % Calculate initial potential energy and volume
    U_initial = strained_box.get_potential();
    V_initial = (strained_box.xmax - strained_box.xmin) * (strained_box.ymax - strained_box.ymin);
    
    % Apply isometric strain
    for polygon = 1:numel(strained_box.polygons)
        strained_box.polygons(polygon).q = strained_box.polygons(polygon).q * (1 + epsilon);
        strained_box.polygons(polygon).vertices_relative = strained_box.polygons(polygon).vertices_relative * (1 + epsilon);
    end
    strained_box.xmin = strained_box.xmin * (1 + epsilon);
    strained_box.xmax = strained_box.xmax * (1 + epsilon);
    strained_box.ymin = strained_box.ymin * (1 + epsilon);
    strained_box.ymax = strained_box.ymax * (1 + epsilon);
    
    % Visualization setup
    num_polygons = numel(strained_box.polygons);
    colors = lines(num_polygons);
    figure;
    hold on;
    axis equal;
    buffer = strained_box.polygons(1).particles_per_side * strained_box.polygons(1).sigma;
    xlim([(1 + epsilon) * strained_box.xmin - buffer, (1 + epsilon) * strained_box.xmax + buffer]);
    ylim([(1 + epsilon) * strained_box.ymin - buffer, (1 + epsilon) * strained_box.ymax + buffer]);
    
    % Relax the system and visualize
    for step = 1:steps
        strained_box.iterate_time(dt);
        
        if mod(step, 10) == 0
            cla;
            rectangle('Position', [strained_box.xmin, strained_box.ymin, strained_box.xmax - strained_box.xmin, strained_box.ymax - strained_box.ymin], ...
                'EdgeColor', 'k', 'LineWidth', 2);
            
            for polygon = 1:num_polygons
                sigma = strained_box.polygons(polygon).sigma;
                vertices = strained_box.polygons(polygon).vertices_relative + strained_box.polygons(polygon).q;
                
                % Apply PBC to the vertices for visualization
                vertices_pbc = strained_box.apply_pbc2d(vertices);
                
                plot([vertices_pbc(:,1); vertices_pbc(1,1)], [vertices_pbc(:,2); vertices_pbc(1,2)], ...
                    '-', 'Color', colors(polygon,:), 'LineWidth', 2);
                
                for j = 1:size(vertices_pbc, 1)
                    rectangle('Position', [vertices_pbc(j,1)-sigma/2, vertices_pbc(j,2)-sigma/2, sigma, sigma], ...
                        'Curvature', [1, 1], 'FaceColor', colors(polygon,:));
                end
            end
            
            title(sprintf('Step: %d', step));
            drawnow;
        end
    end
    
    % Calculate new potential energy and volume
    U_final = strained_box.get_potential();
    V_final = (strained_box.xmax - strained_box.xmin) * (strained_box.ymax - strained_box.ymin);
    
    % Calculate change in potential energy per change in volume
    dUdV = (U_final - U_initial) / (V_final - V_initial);
    
    close all;
end