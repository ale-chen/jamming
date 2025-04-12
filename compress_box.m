%%% NEED TO CHANGE MAX_ENERGY_ITER TO BE MUCH LARGER
%%% ENERGY THRESH MIGHT HAVE TO BE 10 ORDERS OF MAGNITUDE LOWER THAN THE
%%% PRESSURE TARGET
%%% SYSTEM STABILITY CHECK NEEDS TO LOWER THE KE THRESH 


function compress_box(...
    param_json,...
    box_dir,...
    box_file,...
    output_dir,...
    do_visual...
    )
    arguments
        param_json string
        box_dir string
        box_file string
        output_dir string
        do_visual string
    end

    is_json = strcmp(param_json{1}(...
        strlength(param_json)-4:strlength(param_json)),...
        ".json");
    if ~is_json
        ME = MException('Jamming:NotAJson',...
            "%s does not end with .json", param_json{1});
        throw(ME);
    end
    params = jsondecode(fileread(fullfile(box_dir,param_json)));
    loaded_data = load(fullfile(box_dir,box_file));
    box = loaded_data.box;

    % EXAMPLE PARAMS
    % {
    % "N": 500,
    % "sigma": 1,
    % "sides": 3,
    % "bumps_per_side": 2,
    % "k": 1,
    % "dt": 0.001,
    % "max_energy_iter": 1000000,
    % "E_thresh": 0.0000000001, 
    % "P_t": 0.00001,
    % "P_tol_log": 0.001,
    % "r_scale": 0.999
    % }
    
    if strcmpi(do_visual, "true")
        global fig line_kinetic_trans line_kinetic_rot line_potential line_total;
        
        fig = figure('Position', [100, 100, 800, 1000]);
        subplot(2,1,1);
        hold on;
        axis equal;
        buffer = ceil(params.bumps_per_side * sqrt(2)) * params.sigma;
        xlim([box.xmin-buffer, box.xmax+buffer]);
        ylim([box.ymin-buffer, box.ymax+buffer]);
        
        subplot(2,1,2);
        line_kinetic_trans = animatedline('Color', 'r');
        line_kinetic_rot = animatedline('Color', 'g');
        line_potential = animatedline('Color', 'b');
        line_total = animatedline('Color', 'k');
        title('Energy');
        xlabel('Time');
        ylabel('Energy');
        legend('Kinetic (Trans)', 'Kinetic (Rot)', 'Potential', 'Total');
        legend('Location','southwest');
        set(gca, 'YScale', 'log');
    end
    fprintf('Starting binary search...\n N=%d, sides=%d, bumps/side=%d\n',...
                        length(box.polygons), params.sides, params.bumps_per_side');
    if ~exist('total_t', 'var')
        total_t = 0;
    end
    dt = params.dt;
    max_energy_iter = params.max_energy_iter;
    E_thresh = params.E_thresh;
    P_t = params.P_t;
    P_l = 10^(log10(params.P_t) - params.P_tol_log);
    P_h = 10^(log10(params.P_t) + params.P_tol_log);

    P_old = 0;
    
    L = get_L(box);
    L_old = L;

    iteration = 1;

    backup_box = copy(box);
    current_r_scale = params.r_scale;

    max_r_scale = 0.999999999999999999; % WE'RE GOING TOO FINELY
    MAX_BOX_TRIES = 100000;

    success = false;
    give_up = false;
    for box_try = 1:1:MAX_BOX_TRIES
        iteration = iteration + 1;
        energy_minimize(box,E_thresh,dt,max_energy_iter,do_visual);
        [P,overlap_count] = calculate_pressure(box);

        fprintf('Iteration: %d\nCurrent L: %.6f, Current P: %.6f (Target: %.6f), Contacts: %d\nr_scale: %d\n'...
            ,box_try, L, P, P_t, overlap_count, current_r_scale);
            % Search Logic
        if P < P_l % PRESSURE IS TOO LOW, OR UNDERJAMMED
            fprintf('  P < P_l: System is UNDERJAMMED at L = %.6f\n', L);
            backup_box = copy(box); % PRESERVE UNDERJAMMED STATE
            L_old = L;

            L = current_r_scale * L;
        else % P > P_l: We're somewhere near
            if P > P_h || P < P_old % We went too far
                if max_r_scale < current_r_scale
                    % FAILURE, GIVE UP
                    break
                else
                    fprintf('  P > P_h: System is OVERJAMMED at L = %.6f\n', L);
                    current_r_scale = sqrt(current_r_scale); % GO EVEN MORE FINELY              <----- arithmetic mean
                    box = backup_box; % RESET: WE OVERSHOT

                    L = current_r_scale * L_old; % TRY AT NEW RATE
                end
            else
                success = true;
                break
            end
        end
        P_old = P;

        for i = 1:numel(box.polygons)
            box.polygons(i).q = box.polygons(i).q * current_r_scale;
            box.polygons(i).v = box.polygons(i).v * sqrt(current_r_scale);
        end
            box.xmin = -L/2;
            box.xmax = L/2;
            box.ymin = -L/2;
            box.ymax = L/2;

        if strcmpi(do_visual, "true")
            subplot(2,1,2);
            xline(total_t, 'r-', ['L = ' num2str(L, '%.3f')], 'LabelOrientation', 'aligned', 'HandleVisibility', 'off');
        end
        if max_r_scale < current_r_scale
            % FAILURE, GIVE UP
            break
        end

    end

    if success
        fprintf('Jammed state found! Saving final configuration...\n');
        subdirectory = strrep(strrep(datestr(now), ' ', '-'), ':', '-');
        mkdir(fullfile(output_dir, subdirectory));
        save(fullfile(output_dir,subdirectory, 'box.mat'),"box");
        copyfile(param_json,...
        fullfile(output_dir,subdirectory));
    else
        fprintf('Failed...\n');
    end
end

function energy_minimize(box, E_thresh, dt, max_energy_iter, do_visual)
    global line_kinetic_trans line_kinetic_rot line_potential line_total;
    try
        total_t = evalin('base', 'total_t');
    catch
        total_t = 0; % initialize to 0 if it doesn't exist
    end
    visualize = false;
    if strcmpi(do_visual, "true")
        visualize = true;
    end
    colors_small = [0.3, 0.6, 0.9]; % Blue for small polygons
    colors_large = [0.9, 0.3, 0.3]; % Red for large polygons
    N = length(box.polygons);
    colors = zeros(N, 3);
    small_count = ceil(N/2);
    fprintf('Starting energy minimization...\n');

    for i = 1:N
        if i <= small_count
            colors(i,:) = colors_small;
        else
            colors(i,:) = colors_large;
        end
    end
    
    t = 0;
    
    for step = 1:max_energy_iter
        box.iterate_time(dt);
        t = t + dt;
        
        [K_trans, K_rot] = box.get_kinetic();
        total_kinetic = K_trans + K_rot;
        if visualize
            if mod(step, 10) == 0
                V = box.get_potential();
                current_energy = K_trans + K_rot + V;
                [virial_pressure, overlap_count] = calculate_pressure(box);
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
                
                % After calculating other quantities
                %polygon_contacts = box.get_polygon_contacts();
                polygon_contacts = 0;                                                                % FIX THIS!!!!
                
                % Update title to include both contact types
                title(sprintf('t = %.2f, P ~ 10^^%.2f, Asperity Contacts = %d, Polygon Contacts = %d', ...
                    total_t + t, log10(virial_pressure), overlap_count, polygon_contacts));
                
                % Update energy plot with continuous time
                subplot(2,1,2);
                addpoints(line_kinetic_trans, total_t + t, K_trans);
                addpoints(line_kinetic_rot, total_t + t, K_rot);
                addpoints(line_potential, total_t + t, V);
                addpoints(line_total, total_t + t, current_energy);
                
                drawnow;
            end
        end
        
        % check for convergence
        if total_kinetic < E_thresh
            fprintf('Converged at step %d with kinetic energy %.6e\n', step, total_kinetic);
            break;
        end
    end
    assignin('base', 'total_t', total_t + t);
end
function energy = get_energy(box)
    [K_trans, K_rot] = box.get_kinetic();
    V = box.get_potential();
    energy = K_trans + K_rot + V;
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
function L = get_L(box)
    L_x = box.xmax - box.xmin;
    L_y = box.ymax - box.ymin;
    L = (L_x + L_y) / 2;
end
function halfway = geo_mean(L,H)
    halfway = sqrt(L * H);
end