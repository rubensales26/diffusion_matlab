classdef DifusionProblem_2D_1g_FD_O2
    properties
        mesh,
        materials,
        L_mat,
        M_mat,
        k_eff,
        Phi
    end
    
    methods
        function obj = DifusionProblem_2D_1g_FD_O2(mesh, materials)
            obj.mesh = mesh;
            obj.materials = materials;
        end
        
        function obj = AssembleSystem(obj)
            N = obj.mesh.nNodes;
            nx = obj.mesh.nx;
            ny = obj.mesh.ny;
            hx = obj.mesh.hx;
            hy = obj.mesh.hy;
            
            % Sparse triplets preallocation
            max_nz = 5 * N; % We are considering a 3 point stencil
            I = zeros(max_nz,1);
            J = zeros(max_nz,1);
            V = zeros(max_nz,1);
            M_diag = zeros(N,1);
            idx = 0;
            
            for j = 1:nx
                for i = 1:ny
                    k = obj.mesh.getIndexInVector(i, j);
                    
                    % 1. CENTER NODE PROPERTIES
                    x_C = obj.mesh.X(i,j);
                    y_C = obj.mesh.Y(i,j);
                    
                    [D_C, SigA_C, NuSigF_C] = obj.materials.GetProperties(x_C, y_C);
                    
                    % Start building diagonal term (Absorption + Leakage)
                    diag_val = SigA_C; 
                    
                    % --- EAST NEIGHBOR (j+1) ---
                    if j < nx
                        % Properties of Right Neighbor
                        x_E = obj.mesh.X(i,j+1);
                        [D_E, ~, ~] = obj.materials.GetProperties(x_E, y_C);
                        
                        % We represent D_{i+1/2} as the Harmonic Mean
                        D_int = (2 * D_C * D_E) / (D_C + D_E);
                        coeff = -D_int / (hx^2);
                        
                        % Add to matrix
                        idx=idx+1;
                        I(idx)=k; % We are building the equation for node k, so this term belongs in row k
                        J(idx)=obj.mesh.getIndexInVector(i, j+1); % index of the east neighbour
                        V(idx)=coeff;
                        diag_val = diag_val - coeff; % Subtract negative coeff (add positive)
                    else
                        % VACUUM BOUNDARY (Right Wall)
                        % Distance from Center to Wall = hx/2.
                        % Leakage = -D * (Phi_C - 0) / (h/2) -> Coeff = 2*D/h^2
                        diag_val = diag_val + (2*D_C)/(hx^2);
                    end
                    
                    % --- WEST NEIGHBOR (j-1) ---
                    if j > 1
                        x_W = obj.mesh.X(i,j-1);
                        [D_W, ~, ~] = obj.materials.GetProperties(x_W, y_C);
                        
                        D_int = (2 * D_C * D_W) / (D_C + D_W);
                        coeff = -D_int / (hx^2);
                        
                        idx=idx+1;
                        I(idx)=k;
                        J(idx)=obj.mesh.getIndexInVector(i, j-1);
                        V(idx)=coeff;
                        diag_val = diag_val - coeff;
                    else
                        % VACUUM BOUNDARY (Left Wall)
                        diag_val = diag_val + (2*D_C)/(hx^2);
                    end
                    
                    % --- NORTH / SOUTH (Same logic for Y) ---
                    % (Repeat the if/else logic for i-1 and i+1 using hy)
                    % ... (Use the same structure as East/West)

                    % --- NORTH NEIGHBOR (i+1) ---
                    if i < ny
                        y_N = obj.mesh.Y(i+1, j);
                        [D_N, ~, ~] = obj.materials.GetProperties(x_C, y_N);
                        
                        % Harmonic Mean
                        D_int = (2 * D_C * D_N) / (D_C + D_N);
                        coeff = -D_int / (hy^2);
                        
                        idx=idx+1;
                        I(idx)=k;
                        J(idx)=obj.mesh.getIndexInVector(i+1, j);
                        V(idx)=coeff;
                        diag_val = diag_val - coeff;
                    else
                        % VACUUM BOUNDARY (Top Wall)
                        % Distance is hy/2
                        diag_val = diag_val + (2*D_C)/(hy^2);
                    end
                    
                    % --- SOUTH NEIGHBOR (i-1) ---
                    if i > 1
                        y_S = obj.mesh.Y(i-1, j);
                        [D_S, ~, ~] = obj.materials.GetProperties(x_C, y_S);
                        
                        D_int = (2 * D_C * D_S) / (D_C + D_S);
                        coeff = -D_int / (hy^2);
                        
                        idx=idx+1;
                        I(idx)=k;
                        J(idx)=obj.mesh.getIndexInVector(i-1, j);
                        V(idx)=coeff;
                        diag_val = diag_val - coeff;
                    else
                        % VACUUM BOUNDARY (Bottom Wall)
                        diag_val = diag_val + (2*D_C)/(hy^2);
                    end                    
                    % 2. DIAGONAL ELEMENTS
                    idx=idx+1; I(idx)=k; J(idx)=k; V(idx)=diag_val;
                    M_diag(k) = NuSigF_C;
                end
            end
            
            obj.L_mat = sparse(I(1:idx), J(1:idx), V(1:idx), N, N);
            obj.M_mat = spdiags(M_diag, 0, N, N);
        end
        
        function obj = Solve(obj)
             obj = obj.AssembleSystem();
             [eig_vec, eig_val] = eigs(obj.M_mat, obj.L_mat, 1, 'largestabs');
             obj.k_eff = eig_val;
             obj.Phi = reshape(abs(eig_vec), obj.mesh.ny, obj.mesh.nx);
        end
    end
end