function visualize_system(box)
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
    sides = box.polygons(i).sides;
    bumps_per_side = box.polygons(i).particles_per_side;
    axis equal;
    title(['Initial Configuration: ' num2str(N) ' ' num2str(sides) '-sided Polygons']);
    xlabel('X');
    ylabel('Y');
    grid on;

    info_text = sprintf('N=%d, sides=%d, bumps/side=%d', ...
                        N, sides, bumps_per_side);
    annotation('textbox', [0.15, 0.01, 0.7, 0.05], ...
               'String', info_text, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center');
    
    hold off;
end