function test_homo_1EG()
    load("homo_1EG_exact_results.mat", "N_array", "k_analytic", "k_eff_num_vec");

    % PROBLEM PARAMETERS (1G Homogeneous)
    L = 350;
        
    % Physical properties (Index 2 = Homogeneous Fuel)
    D_array       = [1.446, 0.776];
    SigA_array    = [0.0077, 0.0244];
    NuSigF_array  = [0.0000, 0.0260];
    D_val = D_array(2);
    SigA_val = SigA_array(2);
    NuSigF_val = NuSigF_array(2);
    
    for i=1:length(N_array) % Provar només en 100
        N = N_array(i);
        mesh = Mesh_1D_FDM(L, N);
        materials = Materials(1, D_val, SigA_val, NuSigF_val);
        solver = Solver_1D_1EG_FDM(mesh, materials);
        solver = solver.assembleMatrices().solveEigenvalues(1);
    
        % Extract the computed k_eff
        k_eff_num = solver.keff(1);

        % Solver validation
        assert(abs(k_analytic - k_eff_num) < 1e-5); % Pretty bad comparison not to trigger errors
        % Posar comparació flux
    end
end

