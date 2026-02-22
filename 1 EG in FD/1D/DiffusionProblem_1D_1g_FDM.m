% Implementation of the FDM for solving the 1D 1g neutron diffusion
% equation

classdef DiffusionProblem_1D_1g_FDM    
    properties
        mesh
        materials
        A
        Q
        phi
        keff
    end
    
    methods
        function obj = DiffusionProblem_1D_1g_FDM(mesh,materials)
            obj.mesh = mesh;
            obj.materials = materials;
        end

        function obj = assembleMatrices(obj)
            % 1. Parameter extraction
            nNodes = obj.mesh.nNodes;
            % This vector maps each node to its material region index
            region_idx = obj.mesh.material_region_indices; 
            
            % Map properties to every node using the region indices
            % obj.materials.D is mapped to regions
            D_vec = obj.materials.D(region_idx);
            SigA_vec = obj.materials.sigma_a(region_idx);
            NuSigF_vec = obj.materials.nu_sigma_f(region_idx);
            
            % Map the step size (delta_x) to each region to the nodes
            Delta_vec = obj.mesh.delta_region(region_idx);
            
            % 2. Calculate Coupling Coefficients (d)
            % d_i couples node i and i+1.
            % Formula: 2*D1*D2 / (dx2*D1 + dx1*D2)
            D_curr = D_vec(1:end-1);
            D_next = D_vec(2:end);
            dx_curr = Delta_vec(1:end-1);
            dx_next = Delta_vec(2:end);
            
            % Vector of coupling coefficients (size nNodes-1)
            d_vec = (2 .* D_curr .* D_next) ./ (dx_next .* D_curr + dx_curr .* D_next);
            
            % 3. Calculate Boundary Leakage Coefficients (L)
            % Mesh-Centered Zero Flux: L = 2*D / dx
            L_left  = (2 * D_vec(1)) / Delta_vec(1);
            L_right = (2 * D_vec(end)) / Delta_vec(end);
            
            % 4. Build Matrix A (Loss/Transport)
            % Diagonal: Absorption (SigA*dx) + Couplings (d_left + d_right) + Leakage
            diag_A = Delta_vec .* SigA_vec;
            diag_A = diag_A + [0; d_vec] + [d_vec; 0];
            
            % Apply boundary conditions to the first and last diagonal elements
            diag_A(1)   = diag_A(1) + L_left;
            diag_A(end) = diag_A(end) + L_right;
            
            % Off-Diagonals for the tridiagonal system
            upper_diag = [0; -d_vec];
            lower_diag = [-d_vec; 0];
            
            % d_vec has length nNodes-1. We pad it to nNodes.
            padded_d = [d_vec; 0]; 
            padded_upper = [0; -d_vec]; % This works, but let's be explicit:
            obj.A = spdiags([-padded_d, diag_A, padded_upper], -1:1, nNodes, nNodes);

            % 5. Build Matrix Q (Fission/Production)
            % Q_i = dx_i * nu * sigma_f,i
            diag_Q = Delta_vec .* NuSigF_vec;
            obj.Q = spdiags(diag_Q, 0, nNodes, nNodes);
        end

function obj = solveEigenvalues(obj, nm)
            % 1. Solve the generalized eigenvalue problem Q*phi = k*A*phi
            [V, D_eigen] = eigs(obj.Q, obj.A, nm, 'largestabs');
            
            % 2. Extract k-eff
            obj.keff = diag(D_eigen);
            
            % 3. Fix sign and Normalize (Vectorized)
            % Check the sum of the first 5 nodes for each mode
            sinais = sign(sum(V(1:min(5, end), :), 1)); 
            V = V .* sinais;
            
            % Divide each column by its own maximum
            obj.phi = V ./ max(V);
        end
    end
end

