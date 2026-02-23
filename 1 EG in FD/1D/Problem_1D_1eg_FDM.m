% Implementation of the FDM for solving the 1D 1g neutron diffusion
% equation
classdef Problem_1D_1eg_FDM    
    properties
        % Input attributes
        mesh
        materials
        
        % Computed attributes
        A (:,:) double
        Q (:,:) double
        phi (:,:) double
        keff (:,1) double
    end
    
    methods
        function obj = Problem_1D_1eg_FDM(mesh, materials)            
            obj.mesh = mesh;
            obj.materials = materials;
        end
        
        function obj = assembleMatrices(obj)
            % 1. Parameter extraction
            nNodes = obj.mesh.nNodes;
            
            % Map each node to its region
            reg_idx = obj.mesh.node_region_idx; 
            
            % Map physical properties from the Materials object to every node
            D_vec = obj.materials.D_region(reg_idx);
            SigA_vec = obj.materials.sigma_a_region(reg_idx);
            NuSigF_vec = obj.materials.nu_sigma_f_region(reg_idx);
            
            % Map the step size (delta_x) from the Mesh object to every node
            Delta_vec = obj.mesh.delta_region(reg_idx);
            
            % 2. Calculate Coupling Coefficients (d)
            % d_i couples node i and i+1.
            % Formula: 2*D1*D2 / (dx2*D1 + dx1*D2)
            D_curr = D_vec(1:end-1); 
            D_next = D_vec(2:end);   
            dx_curr = Delta_vec(1:end-1);
            dx_next = Delta_vec(2:end);
            
            % Vectorized coupling coefficients
            d_vec = (2 .* D_curr .* D_next) ./ (dx_next .* D_curr + dx_curr .* D_next);
            
            % 3. Calculate Boundary Leakage Coefficients (L)
            % Mesh-Centered Zero Flux: L = 2*D / dx
            L_left  = (2 * D_vec(1)) / Delta_vec(1);
            L_right = (2 * D_vec(end)) / Delta_vec(end);
            
            % 4. Build Matrix A (Loss/Transport)
            diag_A = Delta_vec .* SigA_vec;
            diag_A = diag_A + [0; d_vec] + [d_vec; 0];
            
            % Boundary conditions
            diag_A(1)   = diag_A(1) + L_left;
            diag_A(end) = diag_A(end) + L_right;
            
            % Tridiagonal setup
            padded_d = [d_vec; 0];
            padded_upper = [0; -d_vec];
            obj.A = spdiags([-padded_d, diag_A, padded_upper], -1:1, nNodes, nNodes);
            
            % 5. Build Matrix Q (Fission/Production)
            diag_Q = Delta_vec .* NuSigF_vec;
            obj.Q = spdiags(diag_Q, 0, nNodes, nNodes);
        end
        
        function obj = solveEigenvalues(obj, nm)
            % Solve Q*phi = k*A*phi
            [V, D_eigen] = eigs(obj.Q, obj.A, nm, 'largestabs');
            
            obj.keff = diag(D_eigen);
            
            % Fix sign (ensure flux is positive)
            sinais = sign(sum(V(1:min(5, end), :), 1)); 
            V = V .* sinais;
            
            % Normalize each mode to max = 1.0
            obj.phi = V ./ max(V);
        end
        
        function displayProblem(obj)
            fprintf('=================\n  PROBLEM DATA  \n=================\n\n');
            fprintf('K-effective (Fundamental): %.5f\n', obj.keff(1));
            if length(obj.keff) > 1
                fprintf('K-effective (Harmonics): ');
                fprintf('%.5f ', obj.keff(2:end));
                fprintf('\n');
            end
        end
    end
end