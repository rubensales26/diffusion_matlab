% CLASS: Materials1D1g
% Handles the physical properties for the 1 energy group problem
% Assigns to each cell of the mesh one material and stores its
% cross-sections (D, Sigma_a, Nu_Sigma_f).

classdef Materials1gD1
    properties
        % Input attributes
        materials_vec (:,1) double {mustBePositive,mustBeInteger} % vector with the Ind of the material present in each material region
        sizes_vec (:,1) double {mustBePositive} % Length of each material cell

        % Computed attributes
        material_boundaries (:,1)  % Cartesian coordinates of the boundaries between materials (includes origin)
        D (:,1) % Vector with the D corresponding to each cell
        sigma_a (:,1) % Vector with the sigma_a corresponding to each cell
        nu_sigma_f (:,1) % % Vector with the nu_sigma_f corresponding to each cell
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
    end
end

