classdef Materials1gD1
    %MATERIALS_PROVA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        materials_vec double % Matrix with the Ind of the material present in each cell
        sizes_vec % Length of each material cell
        material_boundaries  % Cartesian coordinates of the boundaries between materials (includes origin)
        D % Vector with the D corresponding to each cell
        sigma_a
        nu_sigma_f
    end
    
    methods
        function obj = Materials1gD1(materials_vec, sizes_vec, D, sigma_a, nu_sigma_f)
            obj.materials_vec = materials_vec;
            obj.sizes_vec = sizes_vec;
            obj.D = D(materials_vec);
            obj.sigma_a = sigma_a(materials_vec);
            obj.nu_sigma_f = nu_sigma_f(materials_vec);

            obj.material_boundaries = [0; cumsum(obj.sizes_vec)];
        end

        function materials_vec_ind = GetMaterialsInd(obj, x_coord)
            % Takes the x_coord and obtains the corresponding index i 
            % of the material cell in which the point is located.
            % For an x_coord between two material cells, the output
            % will be the index of the material on the cell on the right.
            % If the x_coord corresponds to the boundary of the
            % reactor, the output will be empty.

            min_material_boundaries_ind = find(x_coord < obj.material_boundaries, 1);
            materials_vec_ind = min_material_boundaries_ind - 1;
        end

        function [D, sigma_a, nu_sigma_f] = GetProperties(obj, x_coord)
            % Given the x_coord of a point, gets the material properties
            materials_vec_ind = GetMaterialsInd(x_coord);
            mat_ind = obj.materials_vec(materials_vec_ind);
            D = obj.D(mat_ind);
            sigma_a = obj.sigma_a(mat_ind);
            nu_sigma_f = obj.nu_sigma_f(mat_ind);
        end

    end
end

