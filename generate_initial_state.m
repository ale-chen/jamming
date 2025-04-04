function generate_initial_state(...
    output_dir,...
    params_dir,...
    params_json,...
    do_visual...
    )
    arguments
        output_dir string
        params_dir string
        params_json string
        do_visual string
    end
    is_json = strcmp(params_json{1}(...
        strlength(params_json)-4:strlength(params_json)),...
        ".json");
    if ~is_json
        ME = MException('Jamming:NotAJson',...
            "%s does not end with .json", params_json{1});
        throw(ME);
    end
    
    params = jsondecode(fileread(fullfile(params_dir,params_json)));

    %{
    EXAMPLE PARAMS JSON
    {
        "N": 1000,
        "sigma": 1,
        "sides": 3,
        "bumps_per_side": 2,
        "k": 1,
        "dt": 0.001,
        "max_energy_iter": 1000000,
        "E_thresh": 
        "P_t": 0.001,
        "P_tol_log": 3,
        "r_scale": 0.999
    %}
    % x'' + kx/m + 2\zeta\sqrt{k/m} x' = 0 would be ideal,
    % but damping is done on each larger particle, not asperity.
    % however k is done per asperity, so k/m should be asperity m'
    % so m' = m/(bumps_per_side*sides) --> nice that reg_poly has obj.m = m'
    % b_trans should be in terms of total m.
    % b_ang should be in terms of total moment of inertia, however b_ang is for
    % all particles, regardless of actual polygon.mofi
    % bumps_per_side is the smaller one
    % Pressure threshold is on the log scale (arithmetic mean is P_t)
    % e^-0.001 gives 0.999 approximately.

    N = params.N;
    sigma = params.sigma;
    sides = params.sides;
    n1 = params.bumps_per_side;
    n2 = ceil(n1 * sqrt(2)); % sqrt(2) ratio
    m_tot = 1;
    
    zeta = 1;
    
    k = params.k;
    % We want to be close to zeta = 1 (critical). m per asperity =
    % 1/(sides*particles_per_side)
    % 2 * zeta * sqrt(k/m)
    b_trans = 2*zeta*sqrt(k*m_tot);
    % This is tougher, as m of i is different so we just use geom. mean
    m_vec_small = (m_tot/(n1*sides)).* ones(1,n1*sides);
    m_vec_big = (m_tot/(n2*sides)).* ones(1,n2*sides);

    small_poly = regular_polygon(sigma,sides,m_vec_small,[0.0,0.0],[0.0,0.0],0.0,n1);
    big_poly = regular_polygon(sigma,sides,m_vec_big,[0.0,0.0],[0.0,0.0],0.0,n2);
    b_ang = 2*zeta*sqrt(k*sqrt(small_poly.mofi*big_poly.mofi)/2);
    
    [init_polygons_cell,L] = initial_config(N, sigma, sides, n1, m_tot);
    bounds = [-L/2,L/2,-L/2,L/2]; % (X)
        [xmin,xmax,ymin,ymax] = deal(bounds(1),bounds(2),bounds(3),bounds(4));
    
    if ~(size(init_polygons_cell)==N)
        ME2 = MException('Jamming:Initial_Config_Failed',...
            'Failed to find inital state with %d particles, L = %d',N,L);
        throw(ME2);
    end

    box = periodic_box_polygon(init_polygons_cell, k, b_trans, b_ang, xmin, xmax, ymin, ymax);

    if strcmpi(do_visual, "true")
        visualize_system(box,params);
    end
    subdirectory = strrep(strrep(datestr(now), ' ', '-'), ':', '-');
    mkdir(fullfile(output_dir, subdirectory));
    save(fullfile(output_dir,subdirectory, 'box.mat'),"box");
    copyfile(fullfile(params_dir,params_json),...
        fullfile(output_dir,subdirectory));
end
function [polygon_cell, L]...
    = initial_config(N, sigma, sides, bumps_per_side_small, m_tot)
    MAX_ITER = 1000;
    rand_q =@(l) (rand([1,2])-(1/2)).*(l);

    bumps_per_side_big = ceil(bumps_per_side_small * sqrt(2)); % sqrt(2) ratio
    
    m_vec_small = (m_tot/(bumps_per_side_small*sides))...
        .* ones(1,bumps_per_side_small*sides);
    m_vec_big = (m_tot/(bumps_per_side_big*sides))...
        .* ones(1,bumps_per_side_big*sides);
    
    spoke_length = sin(pi*(0.5 - (1/sides)))/sin(2*pi/sides);% in terms of side length
    failure = true;
    for k = 2:MAX_ITER
        polygons = [];
        
        % replace polygons each time
        % k * 2 * spoke length * (1/2) * (large s.l. + small s.l.)
        L = k * sigma * (bumps_per_side_big + bumps_per_side_small)...
            * spoke_length*sqrt(N);
    
        bounds = [-L/2,L/2,-L/2,L/2]; % (X)
        [xmin,xmax,ymin,ymax] = deal(bounds(1),bounds(2),bounds(3),bounds(4));
    
        test_box = periodic_box_polygon({},1,1,1,xmin,xmax,ymin,ymax);
    
        % IS ITERATIVE PLACEMENT LIKE THIS ACTUALLY RANDOM?? idt so...
    
        % small polygons
        for i = 1:ceil(N/2)
            new_polygon = regular_polygon(...
                sigma,sides,m_vec_small,rand_q(L),[0.0,0.0],0.0,bumps_per_side_small...
                );
            new_polygon.rotate_vertices(pi*rand);
            polygons = [polygons, new_polygon];
            for j = 1:MAX_ITER
                if check_potential_overlap(test_box,polygons)
                    polygons(end).q = rand_q(L);
                    new_polygon.rotate_vertices(pi*rand);
                else
                    fprintf(...
                        '\n%d/%d polygons (small) placed in box size %f'...
                        , i, N, L);
                    failure = false;
                    break
                end
            end
            if failure
                break % we failed, break out
            end
        end
        if failure
            continue % we failed, try next k
        end

        
        % big polygons
        for i = ceil(N/2)+1:N
            failure = true;
            new_polygon = regular_polygon(...
                sigma,sides,m_vec_big,rand_q(L),[0.0,0.0],0.0,bumps_per_side_big...
                );
            new_polygon.rotate_vertices(pi*rand);
            polygons = [polygons, new_polygon];
            for j = 1:MAX_ITER
                if check_midpoint_overlap_big(...
                        polygons(end).q,polygons,bounds,N)
                    polygons(end).q = rand_q(L);
                    new_polygon.rotate_vertices(pi*rand);
                elseif check_potential_overlap(test_box,polygons)
                    polygons(end).q = rand_q(L);
                    new_polygon.rotate_vertices(pi*rand);
                else
                    fprintf(...
                        '\n%d/%d polygons (large) placed in box size %f'...
                        , i, N, L);
                    failure = false;
                    break
                end
            end
            if failure
                break % we failed, break out.
            end
        end
        if failure
            continue % we failed, try next k
        else
            fprintf('\nSuccess! %d polygons placed!', N);
            polygon_cell = num2cell(polygons); % we succeeded, return
            return
        end
    end
    polygon_cell = num2cell(polygons); % we failed, return
end

function overlapping = check_potential_overlap(box,polygons)
    overlapping = false;
    box.polygons = polygons;
    [~,~,~,V] = box.iterate_time(0.000000000001);
    if V > 0
        overlapping = true;
    end
end

function overlapping =...
    check_midpoint_overlap_big(position, polygons, bounds, N)
    % once all small polygons have been placed, we run into danger of
    % sqrt(2) s.l. large polygons fully circumscribing small polygon in
    % this case, the minimum permissible midpoint distance between a big
    % polygon and a small polygon is face-to-face small s.l. (2sqrt(2)-2)/3
    % (X)
    n = polygons(1).particles_per_side; % particles per side for small
    % minimum possible non-overlapping distance: base to base
    threshold_small = polygons(1).sigma*(n + ceil(n * sqrt(2)));
    threshold_big = polygons(1).sigma*(2*ceil(n * sqrt(2)));
    
    valid_mask = false(1, length(polygons)-1);
    for i = 1:length(polygons)-1
        valid_mask(i) = ~isempty(polygons(i));
    end
    valid_indices = find(valid_mask);
    num_valid = length(valid_indices);

    positions = zeros(num_valid, length(position));
    distances = zeros(num_valid, 1);

    for i = 1:num_valid
        idx = valid_indices(i);
        positions(i,:) = polygons(idx).q;
    end
    
    for i = 1:num_valid
        distances(i) =...
            norm(pbc_vector_difference(positions(i,:),position,bounds));
    end
    overlapping = any(distances(1:ceil(N/2))<threshold_small);
    if num_valid > (N/2)+1
        overlapping = overlapping...
            || any(distances(ceil(N/2)+1:num_valid)<threshold_big);
    end 
end
function d12 = pbc_vector_difference(q1, q2, bounds)
    % q1 - q2 --> produce vector q2 to q1
    % assume q1 and q2 are 1x2 row vectors
    [xmin,xmax,ymin,ymax] = deal(bounds(1),bounds(2),bounds(3),bounds(4));
    L = [xmax - xmin, ymax - ymin];
    d12 = q1 - q2;
    d12 = d12 - L .* round(d12./L);
end
function visualize_system(box, params)
% (X) totally
    figure('Name', 'Polygon Configuration', 'Position', [100, 100, 600, 600]);
    colors_small = [0.3, 0.6, 0.9]; % Blue for small polygons
    colors_large = [0.9, 0.3, 0.3]; % Red for large polygons

    N = length(box.polygons);
    small_polygons = box.polygons(1:ceil(N/2));
    big_polygons = box.polygons(ceil(N/2)+1:N);

    plot([box.xmin, box.xmax, box.xmax, box.xmin, box.xmin], ...
         [box.ymin, box.ymin, box.ymax, box.ymax, box.ymin], ...
         'k--', 'LineWidth', 1.5);
    color = colors_small;
    for i = 1:N
        polygon = box.polygons(i);
        vertices = polygon.get_vertices();

        if i > ceil(N/2)
            color = colors_large;
        end

        for j = 1:size(vertices, 1)
            r = polygon.sigma/2;
            pos = [vertices(j,1)-r, vertices(j,2)-r, 2*r, 2*r];
            rectangle('Position', pos, 'Curvature', [1,1], 'FaceColor', color, 'EdgeColor', 'none');
        end

        text(polygon.q(1), polygon.q(2), num2str(i), 'FontSize', 8, 'HorizontalAlignment', 'center');
    end

    axis equal;
    title(['Initial Configuration: ' num2str(N) ' ' num2str(params.sides) '-sided Polygons']);
    xlabel('X');
    ylabel('Y');
    grid on;

    info_text = sprintf('N=%d, sides=%d, bumps/side=%d', ...
                        N, params.sides, params.bumps_per_side);
    annotation('textbox', [0.15, 0.01, 0.7, 0.05], ...
               'String', info_text, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center');
    
    hold off;
end