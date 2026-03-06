classdef Materials
    %Materials Physical properties manager
    %   It is independent of the dimensions of the reactor
    properties
        % Input attributes
        region_materials (:,1) double   % Material ID assigned to each region
        D_lib (:,1) double              % Library of Diffusion coefficients 
        sigma_a_lib (:,1) double        % Library of Absorption cross-sections
        nu_sigma_f_lib (:,1) double     % Library of Fission cross-sections
    end
    
    methods
        function obj = Materials(region_materials, D_lib, sigma_a_lib, nu_sigma_f_lib)
            % Store inputs as column vectors
            obj.region_materials = region_materials;
            obj.D_lib = D_lib;
            obj.sigma_a_lib = sigma_a_lib;
            obj.nu_sigma_f_lib = nu_sigma_f_lib;
        end

        function displayMaterials(obj)
            fprintf('==========================\n');
            fprintf('      MATERIAL DATA       \n');
            fprintf('==========================\n\n');

            % 1. Display the Library (Unique materials)
            fprintf('--- Material Library ---\n');
            T_lib = table((1:length(obj.D_lib))', obj.D_lib, obj.sigma_a_lib, obj.nu_sigma_f_lib, ...
                'VariableNames', {'ID', 'D', 'Sigma_a', 'NuSigma_f'});
            disp(T_lib);
        end
    end
end