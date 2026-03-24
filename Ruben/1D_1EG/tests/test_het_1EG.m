function test_het_1EG()
    load("het_1EG_exact_results.mat","size_vec","k_analytic","k_eff_num_vec");

    % PROBLEM PARAMETERS (1G Heterogeneous)
    L = 350;
    x_int1 = 25;   % Coordinate of Interface 1 (Reflector -> Core)
    x_int2 = 325;  % Coordinate of Interface 2 (Core -> Reflector)
    
    % Physical properties
    % Index 1 = Reflector, Index 2 = Homogeneous Fuel
    D_array       = [1.446, 0.776];
    SigA_array    = [0.0077, 0.0244];
    NuSigF_array  = [0.0000, 0.0260];
    
    for i = 1:length(size_vec)  
        N = size_vec(i);
        
        region_lengths = [x_int1; x_int2 - x_int1; L - x_int2];
        region_materials = [1; 2; 1];
        cells_per_region = ceil(N / 3) * ones(3,1);
        
        mesh = Mesh_1D_FDM(region_lengths, cells_per_region);
        
        materials = Materials(region_materials, D_array, SigA_array, NuSigF_array); 
        solver = Solver_1D_1EG_FDM(mesh, materials);
        solver = solver.assembleMatrices().solveEigenvalues(1);
        
        % Extract the computed k_eff
        k_eff_num = solver.keff;

        % Solver validation
        assert(abs(k_analytic - k_eff_num) < 1e-4); % Pretty bad comparison not to trigger errors            
    end
end