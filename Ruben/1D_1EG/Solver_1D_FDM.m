classdef Solver_1D_FDM
    %Solver_1D_1EG_FDM 1D 1EG diffusion equation solver
    %   An implementation based on Hérbert's mesh-centered finite
    %   difference method

    properties
        % Input attributes
        mesh
        materials
        xi_matrix % g x n matrix, n = number of material regions
        
        % Computed attributes
        n_energy_groups
        A (:,:) double                    % Removal/Transport Matrix
        B (:,:) double                    % Fission Source Matrix
        phi (:,:) double                  % Flux profiles
        keff (:,1) double                 % Eigenvalues
        full_mesh_coordinates (:,1) double     % Coordinates for plotting (cells + boundaries)
    end
    
    methods
        function obj = Solver_1D_1EG_FDM(mesh, materials,xi_matrix)
            %Solver_1D_1EG_FDM Constructor of the class
            
            % Check if the input leads to a valid problem
            if length(mesh.region_lengths) ~= length(materials.region_materials)
                error("Invalid input: Arguments region_lengths and region_materials must have the same dimension.");
            end

            obj.mesh = mesh;
            obj.materials = materials;

            g = dims(xi_matrix);
            obj.n_energy_groups = g(1);

            
            % Add the extreme boundary coordinates for plotting purposes
            left_boundary = obj.mesh.region_boundaries(1);
            right_boundary = obj.mesh.region_boundaries(end);
            obj.full_mesh_coordinates = [left_boundary; obj.mesh.cell_centers; right_boundary];
        end
        
        function obj = assembleMatrices(obj)
            %assembleMatrices Assemble the matrices of the system
            %   Vectorized construction of the FDM equations with the
            %   zero flux boundary conditions
            
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
            
            % FULL MATRIX BABY
            nCells = obj.mesh.nCells;
            A = zeros(obj.n_energy_groups*nCells);
            B = zeros(obj.n_energy_groups*nCells);
            
            for i=1:obj.n_energy_groups
                % 2. Calculate Coupling Coefficients (d)
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
                diag_A = Delta_vec .* SigA_vec;
                diag_A = diag_A + [d_vec; 0] + [0; d_vec];
                
                % Add external leakage to the boundary cells
                diag_A(1)   = diag_A(1) + L_left;
                diag_A(end) = diag_A(end) + L_right;
                
                % Off-diagonals (Coupling) must be strictly NEGATIVE
                padded_upper = [0; -d_vec];
                padded_lower = [-d_vec; 0];

                starting_coordinate = (i-1)* nCells + 1;
                A(starting_coordinate:nCells,starting_coordinate:nCells) = spdiags([padded_lower, diag_A, padded_upper], -1:1, nCells, nCells);
                
                diag_Q = Delta_vec .* NuSigF_vec;
                B(starting_coordinate:nCells,starting_coordinate:nCells) = spdiags(diag_Q, 0, nCells, nCells);

            end
            obj.A = A;
            obj.B = B;
        end
        
        function obj = solveEigenvalues(obj, nm)
            % Solve the generalized eigenvalue problem: Q*Phi = keff * A*Phi
            [V, D_eigen] = eigs(obj.Q, obj.A, nm, 'largestabs');
            
            obj.keff = diag(D_eigen);
                        
            % Normalize each mode to max = 1.0
            V_norm = V ./ max(abs(V));

            % If the first mode is mostly negative, flip the sign
            if sum(V_norm(:,1))<0
                V_norm(:,1) = -V_norm(:,1);
            end
            
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