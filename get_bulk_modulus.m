load('jammed_packings/iteration48L8.341P123.4472.mat');

jammed_box = box;

function V = total_potential_energy(box)
    V = box.get_potential();
end

box = jammed_box;
box.b_trans = 100;
box.b_ang = 100;
dt = 0.001;
E_thresh = 1e-8;

min_thresh = 1e-16; % min potential energy per particle threshold
max_steps = 5000; % max number of minimization steps

% seems like the box is too large slightly
epsilon = -.0140568069;
for polygon = 1:numel(box.polygons)
    box.polygons(polygon).q = box.polygons(polygon).q * (1 + epsilon);
    % box.polygons(polygon).vertices_relative = box.polygons(polygon).vertices_relative * (1 + epsilon);
end
box.xmin = box.xmin * (1 + epsilon);
box.xmax = box.xmax * (1 + epsilon);
box.ymin = box.ymin * (1 + epsilon);
box.ymax = box.ymax * (1 + epsilon);

% energy minimization loop
for step = 1:max_steps
    % calculate current potential energy per particle
    V = box.get_potential();
    disp(box.get_potential());
    V_per_particle = V / numel(box.polygons);

    % check if the potential energy is below the threshold
    if V_per_particle < min_thresh
        fprintf('Energy minimization converged after %d steps.\n', step);
        break;
    end

    box.iterate_time(dt);
end

jammed_box = box;
dt = 0.001; % Time step for relaxation
steps = 1000; % Number of relaxation steps

% define the range of strain increments
strain_increments = -logspace(-8, -1, 5);

% store the results
dUdV_values = zeros(size(strain_increments));
virial_pressure_initial = zeros(size(strain_increments));
virial_pressure_final = zeros(size(strain_increments));
packing_fraction_initial = zeros(size(strain_increments));
packing_fraction_final = zeros(size(strain_increments));

for i = 1:length(strain_increments)
    epsilon = strain_increments(i);
    
    % apply isometric strain and calculate dU/dV, virial pressures, and packing fractions
    [dUdV_values(i), virial_pressure_initial(i), virial_pressure_final(i),...
        contact_number_initial, contact_number_final, ~, packing_fraction_initial(i), packing_fraction_final(i)] = ...
        apply_isometric_strain(jammed_box, dt, steps, epsilon);
end

% plotting strain vs virial pressure (initial)
figure;
plot(strain_increments, virial_pressure_initial, 'o-');
xlabel('Strain');
ylabel('Virial Pressure (Initial)');
title('Strain vs Virial Pressure (Initial)');
grid on;
drawnow;

% plotting strain vs virial pressure (final)
figure;
plot(strain_increments, virial_pressure_final, 'o-');
xlabel('Strain');
ylabel('Virial Pressure (Final)');
title('Strain vs Virial Pressure (Final)');
grid on;
drawnow;

% plotting strain vs dU/dV
figure;
plot(strain_increments, dUdV_values, 'o-');
xlabel('Strain');
ylabel('dU/dV');
title('Strain vs dU/dV');
grid on;
drawnow;

% plotting virial pressure (initial) vs dU/dV
figure;
plot(virial_pressure_initial, dUdV_values, 'o-');
xlabel('Virial Pressure (Initial)');
ylabel('dU/dV');
title('Virial Pressure (Initial) vs dU/dV');
grid on;
drawnow;

% plotting virial pressure (final) vs dU/dV
figure;
plot(virial_pressure_final, dUdV_values, 'o-');
xlabel('Virial Pressure (Final)');
ylabel('dU/dV');
title('Virial Pressure (Final) vs dU/dV');
grid on;
drawnow;

diff_dUdV_virial = -dUdV_values - virial_pressure_final;

figure;
loglog(abs(strain_increments), abs(diff_dUdV_virial), 'o-');
xlabel('log(Change in Strain)');
ylabel('log(-dU/dV - P_{virial})');
title('log(-dU/dV - P_{virial}) vs log(Change in Strain)');
grid on;
drawnow;

% plotting strain vs packing fraction (initial)
figure;
plot(strain_increments, packing_fraction_initial, 'o-');
xlabel('Strain');
ylabel('Packing Fraction (Initial)');
title('Strain vs Packing Fraction (Initial)');
grid on;
drawnow;

% plotting strain vs packing fraction (final)
figure;
plot(strain_increments, packing_fraction_final, 'o-');
xlabel('Strain');
ylabel('Packing Fraction (Final)');
title('Strain vs Packing Fraction (Final)');
grid on;
drawnow;

function [dUdV, virial_pressure_initial, virial_pressure_final, contact_number_initial, contact_number_final, strained_box, packing_fraction_initial, packing_fraction_final] = apply_isometric_strain(box, dt, steps, epsilon)
    % Apply an isometric strain. x_i = x_i * (1 + \epsilon),
    % y_i = y_i *(1+ \epsilon), L_x = L_x *(1+\epsilon), L_x = L_x *(1+\epsilon).
    % measure U and V before and after the strain.
    strained_box = box.copy();
    
    % calculate initial potential energy, volume, and virial pressure
    U_initial = strained_box.get_potential();
    V_initial = (strained_box.xmax - strained_box.xmin) * (strained_box.ymax - strained_box.ymin);
    [virial_pressure_initial,contact_number_initial] = calculate_virial_pressure(strained_box);
    packing_fraction_initial = calculate_packing_fraction(strained_box);

    
    % apply isometric strain
    for polygon = 1:numel(strained_box.polygons)
        strained_box.polygons(polygon).q = strained_box.polygons(polygon).q * (1 + epsilon);
        % UNCOMMENT LINE BELOW FOR RELATIVE BUMP RESCALING
        % strained_box.polygons(polygon).vertices_relative = strained_box.polygons(polygon).vertices_relative * (1 + epsilon);
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
    
    for step = 1:steps
        strained_box.iterate_time(dt);
        
        if mod(step, 10) == 0
            cla;
            rectangle('Position', [strained_box.xmin, strained_box.ymin, strained_box.xmax - strained_box.xmin, strained_box.ymax - strained_box.ymin], ...
                'EdgeColor', 'k', 'LineWidth', 2);
            
            for polygon = 1:num_polygons
                sigma = strained_box.polygons(polygon).sigma;
                % PBC to the polygon position for visualization
                q = strained_box.apply_pbc2d(strained_box.polygons(polygon).q);
                vertices = strained_box.polygons(polygon).vertices_relative + q;
                
                plot([vertices(:,1); vertices(1,1)], [vertices(:,2); vertices(1,2)], ...
                    '-', 'Color', colors(polygon,:), 'LineWidth', 2);
                
                for j = 1:size(vertices, 1)
                    rectangle('Position', [vertices(j,1)-sigma/2, vertices(j,2)-sigma/2, sigma, sigma], ...
                        'Curvature', [1, 1], 'FaceColor', colors(polygon,:));
                end
            end
            
            % Visualize contacts based on PBC-wrapped bump positions
            for i = 1:numel(strained_box.polygons)-1
                for j = i+1:numel(strained_box.polygons)
                    poly_i = strained_box.polygons(i);
                    poly_j = strained_box.polygons(j);
                    
                    vertices_i = poly_i.get_vertices();
                    vertices_j = poly_j.get_vertices();
                    
                    for ii = 1:size(vertices_i, 1)
                        for jj = 1:size(vertices_j, 1)
                            vertex_i_pbc = strained_box.apply_pbc2d(vertices_i(ii, :));
                            vertex_j_pbc = strained_box.apply_pbc2d(vertices_j(jj, :));
                            
                            rij = vertex_j_pbc - vertex_i_pbc;
                            distance = norm(rij);
                            
                            sigma = (poly_i.sigma + poly_j.sigma) / 2;
                            
                            if distance < sigma
                                line([vertex_i_pbc(1), vertex_j_pbc(1)], [vertex_i_pbc(2), vertex_j_pbc(2)], ...
                                    'Color', 'r', 'LineWidth', 1.5);
                                
                                % Check if both vertices are being wrapped
                                if any(abs(vertex_i_pbc - vertex_j_pbc) > [strained_box.xmax - strained_box.xmin, strained_box.ymax - strained_box.ymin] / 2)
                                    % Calculate the wrapped positions of both vertices
                                    vertex_i_wrapped = vertex_i_pbc;
                                    vertex_j_wrapped = vertex_j_pbc;
                                    
                                    if vertex_i_pbc(1) < vertex_j_pbc(1)
                                        vertex_i_wrapped(1) = vertex_i_wrapped(1) + (strained_box.xmax - strained_box.xmin);
                                    else
                                        vertex_j_wrapped(1) = vertex_j_wrapped(1) + (strained_box.xmax - strained_box.xmin);
                                    end
                                    
                                    if vertex_i_pbc(2) < vertex_j_pbc(2)
                                        vertex_i_wrapped(2) = vertex_i_wrapped(2) + (strained_box.ymax - strained_box.ymin);
                                    else
                                        vertex_j_wrapped(2) = vertex_j_wrapped(2) + (strained_box.ymax - strained_box.ymin);
                                    end
                                    
                                    % Draw the additional contact line for the wrapped case
                                    line([vertex_i_wrapped(1), vertex_j_wrapped(1)], [vertex_i_wrapped(2), vertex_j_wrapped(2)], ...
                                        'Color', 'r', 'LineWidth', 1.5);
                                end
                            end
                        end
                    end
                end
            end
            
            title(sprintf('Step: %d', step));

            dof = 3 * num_polygons; % 2 translational + 1 rotational DOF per polygon
            [~, contact_number] = calculate_virial_pressure(strained_box);
            
            text(strained_box.xmin, strained_box.ymax + buffer, sprintf('DOF: %d, Contacts: %d', dof, contact_number), ...
                'FontSize', 12, 'VerticalAlignment', 'top');
            
            drawnow;
        end
    end
    
    % calculate new potential energy, volume, and virial pressure
    U_final = strained_box.get_potential();
    V_final = (strained_box.xmax - strained_box.xmin) * (strained_box.ymax - strained_box.ymin);
    [virial_pressure_final, contact_number_final] = calculate_virial_pressure(strained_box);
    packing_fraction_final = calculate_packing_fraction(strained_box);
    
    % calculate change in potential energy per change in volume
    dUdV = (U_final - U_initial) / (V_final - V_initial);
    disp(['Pressure right before return = ', num2str(virial_pressure_final)]);
end

function [pressure, overlap_count] = calculate_virial_pressure(obj)
% can be made to return contact number as well
    virial = 0;
    overlap_count = 0;
    
    % over all unique pairs of polygons
    for i = 1:numel(obj.polygons)-1
        for j = i+1:numel(obj.polygons)
            poly_i = obj.polygons(i);
            poly_j = obj.polygons(j);

            % distance between the polygons
            Rij = obj.pbc_vector_difference(poly_i.q, poly_j.q);
            
            vertices_i = poly_i.get_vertices();
            vertices_j = poly_j.get_vertices();
            
            % over all pairs of vertices between the two polygons
            for ii = 1:size(vertices_i, 1)
                for jj = 1:size(vertices_j, 1)
                    % distance between the two vertices
                    rij = obj.pbc_vector_difference(vertices_i(ii, :), vertices_j(jj, :));
                    distance = norm(rij);
                    
                    sigma = (poly_i.sigma + poly_j.sigma) / 2;
                    
                    % is there is an overlap between the vertices
                    if distance < sigma
                        % contact number
                        overlap_count = overlap_count + 1;

                        % calculate the force between the vertices
                        force = obj.k_int * (sigma - distance) * rij / distance;
                        
                        % add the contribution to the virial
                        % change Rij to rij for old configuration
                        virial = virial + dot(Rij, force);
                        disp('virial:')
                        disp(virial);
                    end
                end
            end
        end
    end
    
    % virial theorem
    area = obj.get_area();
    disp(['Area = ', num2str(area)]);
    pressure = virial / (2 * area);
    disp(['Final calculated pressure = ', num2str(pressure)]);
    return
end

function packing_fraction = calculate_packing_fraction(box)
    total_particle_area = 0;
    total_overlap_area = 0;
    
    for i = 1:numel(box.polygons)
        poly = box.polygons(i);
        vertices = poly.get_vertices();
        
        % calculate the area of the particle
        particle_area = 0;
        for j = 1:size(vertices, 1)
            particle_area = particle_area + pi * (poly.sigma / 2)^2;
        end
        total_particle_area = total_particle_area + particle_area;
        
        % calculate the overlap area with other particles
        for j = i+1:numel(box.polygons)
            other_poly = box.polygons(j);
            other_vertices = other_poly.get_vertices();
            
            for ii = 1:size(vertices, 1)
                for jj = 1:size(other_vertices, 1)
                    rij = box.pbc_vector_difference(vertices(ii, :), other_vertices(jj, :));
                    distance = norm(rij);
                    
                    k1 = poly.sigma / 2;
                    k2 = other_poly.sigma / 2;
                    
                    if distance < k1 + k2
                        gamma = k1 + k2 - distance;
                        alpha = (k1 * gamma) / (k1 + k2);
                        beta = (k2 * gamma) / (k1 + k2);
                        
                        overlap = (2 * acos(alpha / k1) / (2 * pi)) * pi * k1^2 - sqrt(k1^2 - alpha^2) * alpha + ...
                                  (2 * acos(beta / k2) / (2 * pi)) * pi * k2^2 - sqrt(k2^2 - beta^2) * beta;
                        
                        total_overlap_area = total_overlap_area + overlap;
                    end
                end
            end
        end
    end
    
    box_area = (box.xmax - box.xmin) * (box.ymax - box.ymin);
    packing_fraction = (total_particle_area - total_overlap_area) / box_area;
end