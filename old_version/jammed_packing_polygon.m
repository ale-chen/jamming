classdef jammed_packing_polygon
    properties
        % immutable attributes
        polygon_attributes % matrix of polygon properties
        
        % dynamic: will be updated with each time step
        t % current time
        polygons % list of polygons themselves
        box_bounds % box coordinates ((top left x1, top left x2), (bottom right x1, bottom right x2))
        
    end
    methods
        % INIT
        

        % UPDATE TIME
        function iterate_verlet(dt)
        end

        % MICROSTATE VARIABLES
        function energy = get_internal_energy()
            energy;
        end
        
        % MACROSTATE VARIABLES
        function pressure = get_pressure()
            pressure;
        end

        % HELPER FUNCTIONS
        function [force,torque] = calculate_force_torque(polygon_id)
            % calculate total force and torque on a given polygon
        end

        function [x,y] = distance(v1,v2)
            % calculate distance over periodic box between two points
            % v1 - v2
            x;
            y;
        end
    end
end