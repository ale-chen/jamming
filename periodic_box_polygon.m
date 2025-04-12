classdef periodic_box_polygon < handle & matlab.mixin.Copyable
  properties
    % dynamic: will be updated with each time step
    t double % current time
    polygons % list of polygons themselves
    
    % static: stays constant
    k_int double % interaction spring constant
    b_trans double % damping term for translational velocity
    b_ang double % damping term for angular velocity

    dim % spatial dimensions of the system
    xmin % minimum x value
    xmax % maximum x value
    ymin % minimum y value
    ymax % maximum y value

    polygon_contact_pairs
  end

  methods
    % INIT
    function obj = periodic_box_polygon(polygons, k_int, b_trans, b_ang,...
            xmin, xmax, ymin, ymax)
      obj.t = 0;
      obj.polygons = [polygons{:}];
      obj.k_int = k_int;
      obj.b_trans = b_trans;
      obj.b_ang = b_ang;
      obj.xmin = xmin;
      obj.xmax = xmax;
      obj.ymin = ymin;
      obj.ymax = ymax;
    end

    function [pressure, T_trans, T_rot, V] = iterate_time(obj, dt)
        % updates positions of all polygons cofm so they are modded to box
        % position and angle displacement calculations, with damping terms
        % scaled to sum(polygon.m) and i
        % 
        % Verlet Integration Review:
        % x(t+dt) := x(t) + v(t)dt + \frac{1}{2m}F(x(t))dt^2 - bv(t)dt^2
        % v(t + dt) := v(t) + \frac{dt}{2m} (F(x(t)) + F(x(t+dt))) -
        % (b/m)v(t)dt
        % NEEDS TO RETURN PRESSURE CALCULATION
        
        % PERSISTENT VARIABLES THROUGH STEPS/STAGES

        force_rad_old = zeros(size(obj.polygons,2),2);
        torque_old = zeros([size(obj.polygons,2),1]);

        force_rad_new = zeros(size(obj.polygons,2),2);
        torque_new = zeros([size(obj.polygons,2),1]);

        pressure = 0;
        
        % FIRST STEP:
        for i = 1:size(obj.polygons, 2)
            % in the case that this is the first itreation,
            % calculate forces. Otherwise, grab from previous iteration.
            if obj.t == 0
                [force_rad_old(i,:), torque_old(i), ~, ~] = obj.total_force_torque(i);
                force_rad_old(i,:) = force_rad_old(i,:) - (obj.b_trans * obj.polygons(i).v);
                torque_old(i) = torque_old(i,:) - (obj.b_ang * obj.polygons(i).w);
            else
                force_rad_old(i,:) = obj.polygons(i).F_rad;
                torque_old(i) = obj.polygons(i).T;
            end

            obj.polygons(i).q...
                = obj.polygons(i).q... % x_0
                + obj.polygons(i).v * dt... % v_0
                + ((1/(2*obj.polygons(i).mass_tot)) * obj.polygons(i).F_rad) * (dt^2); % APPLY PBC !!NOT DONE YET!! FIGURE OUT IF I WANT TO SCALE TO MASS (ask mark)
            
            changetheta =...
                + obj.polygons(i).w * dt... % w_0
                + ((1/(2*(obj.polygons(i).mofi))) * obj.polygons(i).T) * (dt^2);
            obj.polygons(i).rotate_vertices(changetheta)
            obj.polygons(i).theta = obj.polygons(i).theta + changetheta;
        end
        
        % SECOND STEP:
        for i = 1:size(obj.polygons,2)
            [new_force, new_torque, P_cont, V_check] = obj.total_force_torque(i);
            
            new_force = new_force - (obj.b_trans * obj.polygons(i).v);
            new_torque = new_torque - (obj.b_ang * obj.polygons(i).w);

            obj.polygons(i).F_rad = new_force;
            obj.polygons(i).T = new_torque;
            % now polygons' internal force/torque are F_1 and T_1

            force_rad_new(i,:) = new_force;
            torque_new(i) = new_torque;
            pressure = pressure + P_cont;
        end
        
        V = obj.get_potential();
        
        % STEP 3
        for i = 1:size(obj.polygons,2)
            radial_force_t = force_rad_new(i,:);
            torque_t = torque_new(i);

            obj.polygons(i).v =  obj.polygons(i).v... % v_0
                + (dt/(2*obj.polygons(i).mass_tot))...
                    * (radial_force_t + force_rad_old(i,:)); % F_0, F_1

            
            obj.polygons(i).w = obj.polygons(i).w... % w_0
                + (dt/(2*obj.polygons(i).mofi))...
                    * (torque_t + torque_old(i)); % T_0, T_1
            % now polygons' internal force/torque are v_1, w_1
        end
        % DEBUG
        % disp('=====================')
        % disp(['Snapshot at ', num2str(obj.t)]) 
        % disp(['Polygon ', num2str(i), ': x:', num2str(newq(1)), ', y: ', num2str(newq(2))])
        % disp(['v_x: ', num2str(obj.polygons(i).v(1)), ', v_y: ', num2str(obj.polygons(i).v(2))])
        % disp(['Radial Force: ', num2str(radial_force_tdt(1)), ' ', num2str(radial_force_tdt(2))])
        % disp(['Torque: ', num2str(torque_tdt)])
        % disp('=====================')
        [T_trans, T_rot] = obj.get_kinetic();
        obj.t = obj.t + dt;
        pressure = pressure*(1/(4*obj.get_area));
    end

    function [F_tot, torque, P_cont, V] = total_force_torque(obj, polygon_id)
        % Force is vector valued [x,y]
        % torque: scalar valued (counterclockwise via rhr)
        
        V = 0;
        F_tot = [0,0];
        torque = 0; % return vector torque (shear) with positive/negative contributions
        P_cont = 0;
        % don't iterate over self
        main_polygon = obj.polygons(polygon_id);
        % disp(main_polygon)
        % disp(main_polygon.get_vertices())
        % mask = obj.polygons ~= main_polygon;
        main_vertices = main_polygon.get_vertices();

        % HORRIBLY INEFFICIENT, DOUBLE-COUNTS INTERACTIONS, NOT VECTORIZED
        % forces = zeros(size(main_vertices,1),2);

        for vertex = 1:size(main_vertices, 1)
            r = main_polygon.vertices_relative(vertex, :);
            main_vertex = main_vertices(vertex,:);
            vertex_tot = [0,0]; % don't split force/torque yet
            for j = 1:numel(obj.polygons)
                if j ~= polygon_id

                    rij = obj.pbc_vector_difference(obj.polygons(polygon_id).q...
                        ,obj.polygons(j).q);
                    
                    other_polygon = obj.polygons(j);
                    other_vertices = other_polygon.get_vertices();
            
                    for k = 1:size(other_vertices,1) % see if this can be vectorized
                        other_vertex = other_vertices(k,:);
                        [F12, E12] = obj.particle_vector_force...
                            (other_polygon.sigma,main_polygon.sigma,...
                            other_vertex, main_vertex);
                        vertex_tot = vertex_tot + F12;
                        V = V + E12;
                        P_cont = P_cont + dot(F12,rij);
                    end
                end
            end
            % ADD HERE TO PERPENDICULAR FORCE AND TORQUE
            % calculate torque (right hand rule)
            torque = torque + (r(1) * vertex_tot(2)) - (r(2) * vertex_tot(1)); % can be done better --think about again
            
            % % calculate radial force contribution
            % r_unit = r / norm(r);
            % force_radial_mag = dot(vertex_tot, r_unit);
            % force_radial = r_unit * force_radial_mag;
            % F_rad_tot = F_rad_tot + force_radial;

            % total force instead
            F_tot = F_tot + vertex_tot;

            % DEBUG
            % disp(['vertex ', num2str(vertex),' r_unit: ',num2str(r_unit)])
            % disp(['vertex ', num2str(vertex),' vertex_tot: ',num2str(vertex_tot)])
            % disp(['vertex ', num2str(vertex),' force_radial_mag: ',num2str(force_radial_mag)])
            % disp(['vertex ', num2str(vertex),' force_radial: ',num2str(force_radial)])
            % 
            % disp(['vertex ', num2str(vertex),' F_tot: ',num2str(F_tot)])

            % forces(i,:) = vertex_tot;


            % torque = torque...
            %     + cross(main_polygon.vertices_relative(i,:), vertex_tot); % fix torque calculation
            % vertex_unit_vec = vertex_tot / norm(vertex_tot);
            % F_tot = F_tot + (dot(vertex_tot,vertex_unit_vec) * vertex_unit_vec);  Not correct--should be purely radial
        end

        % for other_polygon = obj.polygons(mask)
        %     other_vertices = other_polygon.get_vertices;
        %     forces = zeros(size(main_vertices,1),2);
        %     for i = 1:size(main_vertices, 1)
        %         main_vertex = main_vertices(i, :);
        %         vertex_tot = [0,0]; % don't split force/torque yet
        % 
        %         for j = 1:size(other_vertices, 1)
        %             other_vertex = other_vertices(j, :);
        %             vertex_tot = vertex_tot...
        %                 + obj.particle_vector_force...
        %                 (other_polygon.sigma,main_polygon.sigma,...
        %                 other_vertex, main_vertex);
        %         end
        %         forces(i,:) = vertex_tot;
        %     end
        % 
        %     main_polygon.vertices_relative(i,:) % not sure if this dimension call is correct
        % end
        % for debug

    end
    function V = get_potential(obj)
        V = 0;
        % potential
        for i = 1:numel(obj.polygons)-1
            poly_i = obj.polygons(i);
            vertices_i = poly_i.get_vertices();
    
            for j = i+1:numel(obj.polygons)
                poly_j = obj.polygons(j);
                vertices_j = poly_j.get_vertices();
    
                for ii = 1:size(vertices_i, 1)
                    for jj = 1:size(vertices_j, 1)
                        rij = obj.pbc_vector_difference...
                            (vertices_i(ii, :), vertices_j(jj, :));
                        distance = norm(rij);
                        sigma = (poly_i.sigma + poly_j.sigma) / 2;
                        H = sigma - distance >= 0;
                        V = V...
                            + 0.5 * obj.k_int * (sigma - distance)^2 * H;
                    end
                end
            end
        end
    end
    function [kinetic_trans, kinetic_rot] = get_kinetic(obj)
        kinetic_trans = 0;
        kinetic_rot = 0;
    
        % kinetic
        for i = 1:numel(obj.polygons)
            poly = obj.polygons(i);
            kinetic_trans = kinetic_trans + 0.5 * sum(poly.m) * dot(poly.v, poly.v);
            kinetic_rot = kinetic_rot +  0.5 * poly.mofi * poly.w^2;
        end
    end

    function [F21, E21] = particle_vector_force(obj, sigma1, sigma2, q2, q1)
        % calculates individual pairwise particle vector force from two
        % positions, assuming both are of sigma diameter.
        % q2 ON q1 \alpha k * (sigma - (q1-q2))
        d12 = obj.pbc_vector_difference(q1,q2); % q1 - q2
        r = norm(d12);
        H = ((sigma1+sigma2)/2) - r >= 0;

        F21 = (d12./r) .*(((sigma1+sigma2)/2) - r) .* obj.k_int .* H;
        E21 = (((sigma1+sigma2)/2) - r)^2 .* obj.k_int .* H / 2;

        F21(isnan(F21) | isinf(F21)) = 0;
    end

    function d12 = pbc_vector_difference(obj, q1, q2)
        % q1 - q2 --> produce vector q2 to q1
        % assume q1 and q2 are 1x2 row vectors
        L = [obj.xmax - obj.xmin, obj.ymax - obj.ymin];
        d12 = q1 - q2;
        d12 = d12 - L .* round(d12./L);
    end

    function area = get_area(obj)
        area = (obj.xmax - obj.xmin) * (obj.ymax - obj.ymin);
    end

    function wrapped_qs = apply_pbc2d(obj, qs) % DON'T DO THIS
        
        % Extract x and y coordinates
        x = qs(:, 1);
        y = qs(:, 2);
        
        x = x - obj.xmin;
        x = mod(x, obj.xmax - obj.xmin);
        x = x + obj.xmin;
        y = y - obj.ymin;
        y = mod(y, obj.ymax - obj.ymin);
        y = y + obj.ymin;
        
        % Combine x and y
        wrapped_qs = [x, y];
    end
    function update_polygon_contacts(obj)
        obj.polygon_contact_pairs = [];  % Clear existing contacts
        
        for i = 1:numel(obj.polygons)-1
            for j = i+1:numel(obj.polygons)
                poly_i = obj.polygons(i);
                poly_j = obj.polygons(j);
                
                vertices_i = poly_i.get_vertices();
                vertices_j = poly_j.get_vertices();
                
                has_contact = false;
                for ii = 1:size(vertices_i, 1)
                    for jj = 1:size(vertices_j, 1)
                        rij = obj.pbc_vector_difference(vertices_i(ii, :), vertices_j(jj, :));
                        distance = norm(rij);
                        sigma = (poly_i.sigma + poly_j.sigma) / 2;
                        
                        if distance < sigma
                            has_contact = true;
                            break;
                        end
                    end
                    if has_contact
                        break;
                    end
                end
                
                if has_contact
                    obj.polygon_contact_pairs = [obj.polygon_contact_pairs; i j];
                end
            end
        end
    end
    function polygon_contacts = get_polygon_contacts(obj)
        obj.update_polygon_contacts();
        polygon_contacts = obj.polygon_contact_pairs;
    end
  end

  methods (Access = protected)
        function newObj = copyElement(obj)
            % Call the superclass method to create a shallow copy
            newObj = copyElement@matlab.mixin.Copyable(obj);
            
            % Deep copy each polygon within the box
            newObj.polygons = cell(size(obj.polygons));
            for i = 1:numel(obj.polygons)
                newObj.polygons{i} = obj.polygons(i).copy();
            end
            newObj.polygons = [newObj.polygons{:}];
        end

  end
end


        % newq_store = zeros(size(obj.polygons,2),2);
        % newtheta_store = zeros(size(obj.polygons,2));
        % 
        % for i = 1:size(obj.polygons,2)
        %     % poly = obj.polygons(i);
        %     % It's essential to make sure that this is the object itself, not a pointer.
        % 
        %     % [radial_force_t, torque_t] = obj.total_force_torque(i); % ADD DAMPING DIRECTLY TO FORCE
        %     [radial_force_t,torque_t] = obj.total_force_torque(i);
        % 
        %     theta_0(i,:) = obj.polygons(i).theta;
        %     % disp('theta_prev')
        %     % disp(theta_prev)
        %     force_torque_0{i,1} = radial_force_t;
        %     force_torque_0{i,2} = torque_t;
        % 
        %     newq_store(i,:) = obj.polygons(i).q + obj.polygons(i).v * dt +...
        %         ((1/(2*obj.polygons(i).mass_tot)) * radial_force_t - obj.b_trans * obj.polygons(i).v) * (dt^2); % APPLY PBC !!NOT DONE YET!! FIGURE OUT IF I WANT TO SCALE TO MASS (ask mark)
        %     newtheta_store(i) = obj.polygons(i).theta + obj.polygons(i).w * dt +...
        %         ((1/(2*(obj.polygons(i).mofi))) * torque_t - obj.b_ang * obj.polygons(i).w) * (dt^2);
        % 
        %     % disp('new_theta')
        %     % disp(newtheta)
        % end
        % for i = 1:size(obj.polygons,2)
        %     obj.polygons(i).q = newq_store(i,:);
        %     % obj.polygons(i).q = obj.apply_pbc2d(newq);
        %     obj.polygons(i).theta = newtheta_store(i);
        % 
        %     [radial_force_tdt, torque_tdt] = obj.total_force_torque(i); % ADD DAMPING DIRECTLY TO FORCE
        %     radial_force_t = force_torque_0{i,1};
        %     torque_t = force_torque_0{i,2};
        %     theta_prev = theta_0(i);
        % 
        %     obj.polygons(i).v = obj.polygons(i).v...
        %         + (dt/(2*obj.polygons(i).mass_tot)) * (radial_force_t + radial_force_tdt)...
        %         - (obj.b_trans/obj.polygons(i).mass_tot) * obj.polygons(i).v .* dt;
        % 
        %     obj.polygons(i).w = obj.polygons(i).w...
        %         + (dt/(2*obj.polygons(i).mofi)) * (torque_t + torque_tdt)...
        %         - (obj.b_ang/obj.polygons(i).mofi) * obj.polygons(i).w * dt;
        %     obj.polygons(i).rotate_vertices(obj.polygons(i).theta - theta_prev);
        % 
        %     radial_force_magnitudes(i) = norm(radial_force_tdt);
        %     % DEBUG
        %     % disp('=====================')
        %     % disp(['Snapshot at ', num2str(obj.t)]) 
        %     % disp(['Polygon ', num2str(i), ': x:', num2str(newq(1)), ', y: ', num2str(newq(2))])
        %     % disp(['v_x: ', num2str(obj.polygons(i).v(1)), ', v_y: ', num2str(obj.polygons(i).v(2))])
        %     % disp(['Radial Force: ', num2str(radial_force_tdt(1)), ' ', num2str(radial_force_tdt(2))])
        %     % disp(['Torque: ', num2str(torque_tdt)])
        %     % disp('=====================')
        % end
% 
%     % UPDATE TIME
%     function iterate_verlet(obj, dt)
%       % calculate initial forces
%       forces = zeros(length(obj.polygons), 2);
%       for i = 1:length(obj.polygons)
%         forces(i,:) = obj.calculate_force(i);
%       end
% 
%       % update positions
%       for i = 1:length(obj.polygons)
%         obj.polygons(i).q = obj.polygons(i).q + dt * obj.polygons(i).v + ...
%           0.5 * dt^2 / sum(obj.polygons(i).m) * forces(i,:);
% 
%         % apply periodic boundary conditions
%         obj.polygons(i).q = obj.apply_periodic_boundary(obj.polygons(i).q);
% 
%         % Print the updated position of each polygon
%         fprintf('Polygon %d position: [%.2f, %.2f]\n', i, obj.polygons(i).q(1), obj.polygons(i).q(2));
%       end
% 
%       % update vertices based on new positions
%       for i = 1:length(obj.polygons)
%         obj.polygons(i).vertices = obj.polygons(i).get_vertices();
%       end
% 
%       % calculate new forces
%       new_forces = zeros(length(obj.polygons), 2);
%       for i = 1:length(obj.polygons)
%         new_forces(i,:) = obj.calculate_force(i);
%       end
% 
%       % update velocities
%       for i = 1:length(obj.polygons)
%         obj.polygons(i).v = obj.polygons(i).v + ...
%           0.5 * dt / sum(obj.polygons(i).m) * (forces(i,:) + new_forces(i,:));
% 
%         % print the updated velocity of each polygon
%         fprintf('Polygon %d velocity: [%.2f, %.2f]\n', i, obj.polygons(i).v(1), obj.polygons(i).v(2));
%       end
% 
%       obj.t = obj.t + dt;
%       fprintf('Time step: %.3f\n', obj.t);
%     end
% 
%     % MICROSTATE VARIABLES
%     function energy = get_internal_energy(obj)
%       energy = 0;
%       for i = 1:length(obj.polygons)
%         for j = i+1:length(obj.polygons)
%           r_ij = obj.distance(obj.polygons(i).q, obj.polygons(j).q);
%           if r_ij < (obj.polygons(i).sigma + obj.polygons(j).sigma) / 2
%             energy = energy + 0.5 * obj.k_int * ((obj.polygons(i).sigma + obj.polygons(j).sigma) / 2 - r_ij)^2;
%           end
%         end
%       end
%     end
% 
%     % MACROSTATE VARIABLES
%     function pressure = get_pressure(obj)
%       force = zeros(size(obj.box_bounds));
%       for i = 1:length(obj.polygons)
%         force = force + obj.calculate_force(i);
%       end
%       pressure = norm(force) / (norm(obj.box_bounds(2,:) - obj.box_bounds(1,:)));
% 
%       %% Fixed boundaries
%       virial = 0;
%       for i = 1:length(obj.polygons)
%         virial = virial + dot(r(i), f(i));
%       end
% 
%       %% PBC 
%       virial = 0;
%       for i = 1:length(obj.polygons)-1
%           for j = i:length(obj.polyogns);
%               virial = virial + dot(r_ij(i), f_ij(i)); % 
%           end
%       end
% 
%       %%
%     end
% 
%     % HELPER FUNCTIONS
%     function force = calculate_force(obj, polygon_id)
%       force = [0, 0];
% 
%       for i = 1:length(obj.polygons)
%         if i ~= polygon_id
%           r_ij = obj.distance(obj.polygons(i).q, obj.polygons(polygon_id).q);
% 
%           if r_ij < (obj.polygons(i).sigma + obj.polygons(polygon_id).sigma) / 2
%             direction = (obj.polygons(polygon_id).q - obj.polygons(i).q) / r_ij;
%             force = force + obj.k_int * ((obj.polygons(i).sigma + obj.polygons(polygon_id).sigma) / 2 - r_ij) * direction;
%           end
%         end
%       end
%     end
% 
%     function r_ij = distance(obj, q_i, q_j)
%       dq = q_i - q_j;
%       dq = obj.apply_periodic_boundary(dq);
%       r_ij = norm(dq);
%     end
% 
%     function v = apply_periodic_boundary(obj, v)
%       Lx = obj.box_bounds(2,1) - obj.box_bounds(1,1);
%       Ly = obj.box_bounds(2,2) - obj.box_bounds(1,2);
% 
%       v(1) = v(1) - Lx * round(v(1) / Lx);
%       v(2) = v(2) - Ly * round(v(2) / Ly);
%     end

    % function pressure = iterate_time(obj, dt)
    %     % updates positions of all polygons cofm so they are modded to box
    %     % position and angle displacement calculations, with damping terms
    %     % scaled to sum(polygon.m) and i
    %     % 
    %     % Verlet Integration Review:
    %     % x(t+dt) := x(t) + v(t)dt + \frac{1}{2m}F(x(t))dt^2 - bv(t)dt^2
    %     % v(t + dt) := v(t) + \frac{dt}{2m} (F(x(t)) + F(x(t+dt))) -
    %     % (b/m)v(t)dt
    %     % NEEDS TO RETURN PRESSURE CALCULATION
    %     radial_force_magnitudes = zeros(size(obj.polygons));
    % 
    %     for i = 1:size(obj.polygons,2)
    %         % poly = obj.polygons(i);
    %         % It's essential to make sure that this is the object itself, not a pointer.
    % 
    %         [radial_force_t, torque_t] = obj.total_force_torque(i); % ADD DAMPING DIRECTLY TO FORCE
    % 
    %         theta_prev = obj.polygons(i).theta;
    %         % disp('theta_prev')
    %         % disp(theta_prev)
    % 
    %         mass_tot = sum(obj.polygons(i).m);
    % 
    %         newq = obj.polygons(i).q + obj.polygons(i).v * dt +...
    %             ((1/(2*mass_tot)) * radial_force_t - obj.b_trans * obj.polygons(i).v) * (dt^2); % APPLY PBC !!NOT DONE YET!! FIGURE OUT IF I WANT TO SCALE TO MASS (ask mark)
    %         newtheta = obj.polygons(i).theta + obj.polygons(i).w * dt +...
    %             ((1/(2*(obj.polygons(i).mofi))) * torque_t - obj.b_ang * obj.polygons(i).w) * (dt^2);
    % 
    %         % disp('new_theta')
    %         % disp(newtheta)
    % 
    %         obj.polygons(i).q = newq;
    %         % obj.polygons(i).q = obj.apply_pbc2d(newq);
    %         obj.polygons(i).theta = newtheta;
    % 
    %         [radial_force_tdt, torque_tdt] = obj.total_force_torque(i); % ADD DAMPING DIRECTLY TO FORCE
    % 
    %         obj.polygons(i).v = obj.polygons(i).v...
    %             + (dt/(2*mass_tot)) * (radial_force_t + radial_force_tdt)...
    %             - (obj.b_trans/mass_tot) * obj.polygons(i).v .* dt;
    % 
    %         obj.polygons(i).w = obj.polygons(i).w...
    %             + (dt/(2*obj.polygons(i).mofi)) * (torque_t + torque_tdt)...
    %             - (obj.b_ang/obj.polygons(i).mofi) * obj.polygons(i).w * dt;
    %         obj.polygons(i).rotate_vertices(obj.polygons(i).theta - theta_prev);
    % 
    %         radial_force_magnitudes(i) = norm(radial_force_tdt);
    %         % DEBUG
    %         % disp('=====================')
    %         % disp(['Snapshot at ', num2str(obj.t)]) 
    %         % disp(['Polygon ', num2str(i), ': x:', num2str(newq(1)), ', y: ', num2str(newq(2))])
    %         % disp(['v_x: ', num2str(obj.polygons(i).v(1)), ', v_y: ', num2str(obj.polygons(i).v(2))])
    %         % disp(['Radial Force: ', num2str(radial_force_tdt(1)), ' ', num2str(radial_force_tdt(2))])
    %         % disp(['Torque: ', num2str(torque_tdt)])
    %         % disp('=====================')
    %     end
    % 
    % 
    %     obj.t = obj.t + dt;
    %     pressure = sum(radial_force_magnitudes)*(1/(2*obj.get_area));
    % end