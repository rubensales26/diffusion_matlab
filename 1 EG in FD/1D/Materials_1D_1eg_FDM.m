classdef Materials_1D_1eg_FDM
    properties
        % Input attributes
        region_materials (:,1) double   % Material ID assigned to each region
        D_lib (:,1) double              % Library of Diffusion coefficients 
        sigma_a_lib (:,1) double        % Library of Absorption cross-sections
        nu_sigma_f_lib (:,1) double     % Library of Fission cross-sections
        
        % Computed attributes (Properties mapped to regions)
        D_region (:,1) double           
        sigma_a_region (:,1) double     
        nu_sigma_f_region (:,1) double  
    end
    
    methods
        function obj = Materials_1D_1eg_FDM(region_materials, D_lib, sigma_a_lib, nu_sigma_f_lib)
            % Store inputs as column vectors
            obj.region_materials = region_materials;
            obj.D_lib = D_lib;
            obj.sigma_a_lib = sigma_a_lib;
            obj.nu_sigma_f_lib = nu_sigma_f_lib;
            
            % Map the library values to the actual regions
            obj.D_region = obj.D_lib(obj.region_materials);
            obj.sigma_a_region = obj.sigma_a_lib(obj.region_materials);
            obj.nu_sigma_f_region = obj.nu_sigma_f_lib(obj.region_materials);
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

            % 2. Display the Regional Mapping
            fprintf('--- Regional Mapping ---\n');
            T_reg = table((1:length(obj.region_materials))', obj.region_materials, ...
                obj.D_region, obj.sigma_a_region, obj.nu_sigma_f_region, ...
                'VariableNames', {'Region', 'Material_ID', 'D', 'Sigma_a', 'NuSigma_f'});
            disp(T_reg);
        end
    end
end