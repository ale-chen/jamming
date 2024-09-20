classdef regular_polygon < handle & matlab.mixin.Copyable
    properties
    sigma % vertex diameter (float)
    sides % sides (int)
    k % vertex spring constant (float)
    m % vertex mass (N vector)
    q % position c of m (float)
    v % translational velocity (float)
    theta % angle
    w % angular velocity (float)
    mofi % moment of inertia
    vertices % individual vertex locations
    vertices_relative % individual vertex locations RELATIVE TO Q
    particles_per_side % self explanatory

    torque % total torque calculation
    force % total force calculation
    end

    methods
        function obj = regular_polygon(sigma, sides,...
            m, q, v, w, particles_per_side)
            obj.sigma = sigma;
            obj.sides = sides;
            obj.m = m;
            obj.q = q;
            obj.v = v;
            obj.theta = 0;
            obj.w = w;
            obj.vertices_relative =...
                obj.generate_vertices(sigma, sides, particles_per_side, m);
            obj.vertices = obj.get_vertices();
            obj.mofi = obj.calculate_moment_of_inertia();
            obj.particles_per_side = particles_per_side;
        end
        function vertices = get_vertices(obj)
            vertices = obj.vertices_relative + obj.q;
        end
        function p = get_translational_momentum(obj)
            p = obj.v * sum(obj.m);
        end
        function l = get_angular_momentum(obj)
            l = obj.mofi * obj.w;
        end

        function rotate_vertices(obj, theta)
            R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

            % apply rotation to all vertices
            obj.vertices_relative = obj.vertices_relative * R';
            
            % update the absolute vertex positions
            obj.vertices = obj.get_vertices();
        end

        function v = rotate2d(obj, theta,v)
            if (size(v)==size([1,1]))
                v = ([cos(theta),-sin(theta); sin(theta), cos(theta)] * v')';
            elseif(size(v)==size([1;1]))
                v = [cos(theta),-sin(theta); sin(theta), cos(theta)] * v;
            end
        end
    end
    methods (Access = protected)
        function newObj = copyElement(obj)
            % Call the superclass copy method to get a shallow copy
            newObj = copyElement@matlab.mixin.Copyable(obj);
        end
    end

    methods (Access = private)
        function mofi = calculate_moment_of_inertia(obj)
            % Ensure obj.m is a column vector
            if size(obj.m, 1) == 1
                obj.m = obj.m';
            end
            
            squared_distances = sum(obj.vertices_relative.^2, 2);
            mofi = sum(obj.m .* squared_distances);
        end
        function vertices = generate_vertices(obj, sigma, sides, particles_per_side, m)
            vertices = zeros(sides * particles_per_side, 2);
            r = sigma * particles_per_side / (2 * sin(pi / sides));
            

            % generate subsequent vertices

            for side = 1:sides
                theta = (2 * pi / sides); % NOT THE OBJ.THETA
                vertices((side-1)*particles_per_side + 1, :) = obj.rotate2d((side - 1) * theta, [r; 0]);
                phi = ((pi + theta)/2 ) + ((side - 1) * theta);
                v = obj.rotate2d(phi,[sigma;0]);
                % v = rotation * [sigma; 0];
                for particle = 2:particles_per_side
                    % vertices((side - 1)*particles_per_side + particle, :) = ...
                    %     vertices((side - 1)*particles_per_side + particle - 1, :) + v';
                    vertices((side-1)*particles_per_side + particle, :) = ...
                        vertices((side - 1)*particles_per_side + particle - 1, :) + v';
                end
            end

            % for side = 1:sides
            %     theta = (2 * pi / sides);
            %     vertices((side-1)*particles_per_side + 1, :) = obj.rotate2d((side - 1) * theta, [r; 0]);
            % 
            %     v = obj.rotate2d(((side+1) * theta) - (pi/2),[sigma;0]);
            %     % v = rotation * [sigma; 0];
            %     for particle = 2:particles_per_side % BROKEN RIGHT NOW
            %         % vertices((side - 1)*particles_per_side + particle, :) = ...
            %         %     vertices((side - 1)*particles_per_side + particle - 1, :) + v';
            %         vertices((side-1)*particles_per_side + particle, :) = ...
            %             vertices((side - 1)*particles_per_side + particle - 1, :) + v';
            %     end
            % end
            
            % Center by center of mass (1 x N) * (N x dim) ./ (1)
            center_of_mass = m*vertices ./ sum(m);
            vertices = vertices - center_of_mass;
        end

    end
end

% % Set the polygon parameters
% sigma = 1;
% sides = 3;
% particles_per_side = 1;
% % Generate the polygon vertices
% vertices = polygon(sigma, sides, particles_per_side);
% % Plot the polygon edges
% figure;
% plot(vertices(:,1), vertices(:,2), 'b-', 'LineWidth', 2);
% hold on;
% % Plot the particles as circles
% for i = 1:size(vertices, 1)
%  rectangle('Position', [vertices(i,1)-sigma/2, vertices(i,2)-sigma/2, sigma, sigma], ...
%  'Curvature', [1, 1], 'FaceColor', 'r');
% end
% % Set the axis limits and aspect ratio
% axis equal;
% xlim([-10, 10]);
% ylim([-10, 10]);
% % Add labels and title
% xlabel('X');
% ylabel('Y');
% title('Regular Polygon with Particles');
