classdef Solver_1D_1EG_FDM
    %Solver_1D_1EG_FDM 1D 1EG diffusion equation solver
    %   An implementation based on Hérbert's mesh-centered finite
    %   difference method

    properties
        % Input attributes
        mesh
        materials
        
        % Computed attributes
        A (:,:) double                    % Removal/Transport Matrix
        Q (:,:) double                    % Fission Source Matrix
        phi (:,:) double                  % Flux profiles
        keff (:,1) double                 % Eigenvalues
        full_mesh_coordinates (:,1) double     % Coordinates for plotting (cells + boundaries)
    end
    
    methods
        function obj = Solver_1D_1EG_FDM(mesh, materials)            
            obj.mesh = mesh;
            obj.materials = materials;
            
            % Add the extreme boundary coordinates for plotting purposes
            left_boundary = obj.mesh.region_boundaries(1);
            right_boundary = obj.mesh.region_boundaries(end);
            obj.full_mesh_coordinates = [left_boundary; obj.mesh.cell_centers; right_boundary];
        end
        
        function obj = assembleMatrices(obj)
            nCells = obj.mesh.nCells;
            
            % 1. MAPPING: Link Geometry to Physics
            % Get the geometric region index for each cell
            geom_reg_idx = obj.mesh.cell_region_idx;
            
            % Ask the Materials class which Material ID corresponds to each region
            mat_id_vec = obj.materials.region_materials(geom_reg_idx);
            
            % Map physical properties to every cell
            D_vec = obj.materials.D_lib(mat_id_vec);
            SigA_vec = obj.materials.sigma_a_lib(mat_id_vec);
            NuSigF_vec = obj.materials.nu_sigma_f_lib(mat_id_vec);
            Delta_vec = obj.mesh.delta_x;
            
            % 2. Calculate Coupling Coefficients (d)
            % Harmonically averaged diffusion coefficient at the cell interfaces
            D_curr = D_vec(1:end-1);
            D_next = D_vec(2:end);   
            dx_curr = Delta_vec(1:end-1);
            dx_next = Delta_vec(2:end);
            
            d_vec = (2 .* D_curr .* D_next) ./ (dx_next .* D_curr + dx_curr .* D_next);
            
            % 3. Calculate Boundary Leakage Coefficients (L)
            % Cell-Centered Zero Flux at the outer edges
            L_left  = (2 * D_vec(1)) / Delta_vec(1);
            L_right = (2 * D_vec(end)) / Delta_vec(end);
            
            % 4. Build Matrix A (Loss/Transport)
            % Main diagonal (Removal) must be strictly POSITIVE
            diag_A = Delta_vec .* SigA_vec;
            diag_A = diag_A + [d_vec; 0] + [0; d_vec];
            
            % Add external leakage to the boundary cells
            diag_A(1)   = diag_A(1) + L_left;
            diag_A(end) = diag_A(end) + L_right;
            
            % Off-diagonals (Coupling) must be strictly NEGATIVE
            padded_upper = [0; -d_vec];
            padded_lower = [-d_vec; 0];
            obj.A = spdiags([padded_lower, diag_A, padded_upper], -1:1, nCells, nCells);
            
            % 5. Build Matrix Q (Fission/Production)
            % Fission production must be POSITIVE
            diag_Q = Delta_vec .* NuSigF_vec;
            obj.Q = spdiags(diag_Q, 0, nCells, nCells);
        end
        
        function obj = solveEigenvalues(obj, nm)
            % Solve the generalized eigenvalue problem: Q*Phi = keff * A*Phi
            [V, D_eigen] = eigs(obj.Q, obj.A, nm, 'largestabs');
            
            obj.keff = diag(D_eigen);
            
            % Fix sign (ensure fundamental mode flux is predominantly positive)
            sinais = sign(sum(V(1:min(5, end), :), 1)); 
            V = V .* sinais;
            
            % Normalize each mode to max = 1.0
            V_norm = V ./ max(abs(V));
            
            % Add the boundary condition points to each mode (0 flux at edges)
            obj.phi = zeros(size(V_norm) + [2, 0]); 
            obj.phi(1,:) = 0;
            obj.phi(end,:) = 0;
            obj.phi(2:end-1,:) = V_norm;
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
        
        function plotPhi(obj)
            figure
            hold on
            grid on
            plot(obj.full_mesh_coordinates, obj.phi, 'LineWidth', 1.5);
            xlabel('Position x (cm)');
            ylabel('Normalized Flux \Phi(x)');
            title('1D Reactor Flux Profile');
        end
    end
end