% CLASS: Materials1D1g
% Handles the physical properties for the 1 energy group problem
% Assigns to each cell of the mesh one material and stores its
% cross-sections (D, Sigma_a, Nu_Sigma_f).

classdef Materials1D1g
    properties
        materials_vec % Stores which Material ID is in which material region
        sizes_vec     % Physical width of each material region
        material_boundaries % Cartesian coordinates of boundaries between regions
        D             % Diffusion coefficient mapped to each Region
        sigma_a       % Absorption mapped to each Region
        nu_sigma_f    % Fission mapped to each Region
    end
    
    methods
        function obj = Materials1D1g(materials_vec, sizes_vec, D_library, sigA_library, nuSigF_library)
            obj.materials_vec = materials_vec;
            obj.sizes_vec = sizes_vec;
            
            % Map the library values to the actual regions
            % If materials_vec is [1; 2; 1], D becomes [D(1); D(2); D(1)]
            obj.D = D_library(materials_vec);
            obj.sigma_a = sigA_library(materials_vec);
            obj.nu_sigma_f = nuSigF_library(materials_vec);
            
            % Create the vector of boundaries
            obj.material_boundaries = [0; cumsum(obj.sizes_vec)];
        end

        function [D, sigma_a, nu_sigma_f] = GetProperties(obj, x_coord)
            % Finds the cell index for a specific coordinate
            idx = obj.GetMaterialsInd(x_coord);
            
            % Access properties mapped to that cell
            D = obj.D(idx);
            sigma_a = obj.sigma_a(idx);
            nu_sigma_f = obj.nu_sigma_f(idx);
        end

        function idx = GetMaterialsInd(obj, x_coord)
            % Returns the index of the material region containing x_coord
            if x_coord < 0 || x_coord > obj.material_boundaries(end)
                idx = []; % Outside reactor
            else
                idx = find(x_coord < obj.material_boundaries, 1) - 1;
            end
        end
    end
end