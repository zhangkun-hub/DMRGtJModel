classdef SquareLattice
    properties
        Ly % Number of sites in the y-direction
        Lx % Number of sites in the x-direction
        a % lattice constant 
    end

    methods
        function obj = SquareLattice(Ly, Lx)
            obj.Ly = Ly;
            obj.Lx = Lx;
            obj.a = 1; 
        end

        function [x_coord, y_coord] = indexToCoord(obj, site_idx)
            % Transform site index to square lattice coordinates
            % site index count from 0
            x = fix(site_idx / obj.Ly) + 1;
            y = mod(site_idx, obj.Ly) + 1;
            x_coord = x * obj.a; % Offset every other row
            y_coord = y * obj.a;
        end

        function site_idx = coordToIndex(obj, x, y)
            site_idx = (x /obj.a - 1) * obj.Ly + (y /obj.a - 1);
        end

        function drawBonds(obj, bond_color, bond_width)

            hold on;


            for x = 1:obj.Lx
                for y = 1:obj.Ly
                    % Coordinates for the current site with staggered rows
                    x_coord = x; % Offset every other row
                    y_coord = y * obj.a;

                    % Horizontal bonds
                    if x < obj.Lx 
                        line([x_coord, x_coord + obj.a], [y_coord, y_coord], 'Color', bond_color, 'LineWidth', bond_width);
                    end

                    % \
                    if y < obj.Ly
                        x_vert = x_coord;
                        y_vert = y_coord+obj.a;
                        line([x_coord, x_vert], [y_coord, y_vert], 'Color', bond_color, 'LineWidth', bond_width);
                    end
                end
            end
            hold off;
        end

        function drawPartBonds(obj, bond_color, bond_width, x_start, x_end)
            hold on;
            for x = x_start+1:x_end+1
                for y= 1:obj.Ly-1
                    line([x,x],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
                    if x == x_end+1
                        line([x+1,x+1],[y,y+1],'color',bond_color,'linewidth',bond_width);  hold on;
                    end
                end
                for y = 1:obj.Ly
                    line([x,x+1],[y,y],'color',bond_color,'linewidth',bond_width);  hold on;
                end
            end
            hold off;
        end


        function drawOnsiteData(obj, data, circle_scale, positive_color, negative_color)
            % Draw onsite data as circles on the lattice
            hold on;
            for j = 1:size(data, 1)
                site_idx = data(j, 1);
                value = data(j, 2);
                [x_coord, y_coord] = obj.indexToCoord(site_idx);

                % Define radius for correlation circles
                radius = sqrt(abs(value)) * circle_scale;

                % Draw circle based on value sign
                if real(value) >= 0
                    color = positive_color;
                else
                    color = negative_color;
                end

                rectangle('Position', [x_coord - radius, y_coord - radius, 2 * radius, 2 * radius], ...
                          'Curvature', [1, 1], 'FaceColor', color, 'EdgeColor', 'none');
            end
            hold off;
        end
    end
end