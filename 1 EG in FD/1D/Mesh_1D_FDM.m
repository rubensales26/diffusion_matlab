classdef Mesh_1D_FDM
    properties
        % Input attributes
        region_lengths (:,1) double     % Length of each region
        region_materials (:,1) double   % Material ID assigned to each region
        nodes_per_region (:,1) double   % Number of internal nodes in each region
        
        % Computed attributes
        nNodes (1,1) double             % Total number of nodes in the mesh
        nodes_vec (:,1) double          % Cartesian coordinates of the nodes
        node_region_idx (:,1) double    % The region index each node belongs to
        node_material_idx (:,1) double  % The material ID each node contains
        delta_region (:,1) double       % The step size (dx) for each region
        region_boundaries (:,1) double  % Cartesian coordinates of the region interfaces
    end
    
    methods
        function obj = Mesh_1D_FDM(region_lengths, region_materials, nodes_per_region)
            % Store inputs as column vectors
            obj.region_lengths = region_lengths;
            obj.region_materials = region_materials;
            obj.nodes_per_region = nodes_per_region;
            
            % Compute region boundaries and total nodes
            obj.region_boundaries = [0; cumsum(obj.region_lengths)];
            obj.nNodes = sum(obj.nodes_per_region);
            
            % Pre-allocate arrays
            obj.nodes_vec = zeros(obj.nNodes, 1);
            obj.node_region_idx = zeros(obj.nNodes, 1);
            obj.node_material_idx = zeros(obj.nNodes, 1);
            obj.delta_region = zeros(length(obj.region_lengths), 1);
            
            % Populate the nodes
            ind = 1; 
            for i = 1:length(obj.region_lengths)
                % Step size for a given region: length / (nodes + 1)
                obj.delta_region(i) = obj.region_lengths(i) / (obj.nodes_per_region(i) + 1);
                left_boundary = obj.region_boundaries(i);
                
                for j = 1:obj.nodes_per_region(i)
                    obj.nodes_vec(ind) = left_boundary + obj.delta_region(i) * j;
                    obj.node_region_idx(ind) = i; 
                    obj.node_material_idx(ind) = obj.region_materials(i); 
                    ind = ind + 1;
                end
            end
        end

        function displayMesh(obj)
            fprintf('==========================\n');
            fprintf('        MESH DATA         \n');
            fprintf('==========================\n\n');
            
            fprintf('Total Reactor Length: %.2f cm\n', obj.region_boundaries(end));
            fprintf('Total Number of Nodes: %d\n', obj.nNodes);
            fprintf('Number of Regions:     %d\n\n', length(obj.region_lengths));
            
            % Display Table of Regions
            fprintf('--- Region Breakdown ---\n');
            Region = (1:length(obj.region_lengths))';
            Material_ID = obj.region_materials;
            Length = obj.region_lengths;
            Nodes = obj.nodes_per_region;
            Delta_x = obj.delta_region;
            
            T = table(Region, Material_ID, Length, Nodes, Delta_x);
            disp(T);
            
            % Optional: Display first and last node coordinates to verify boundaries
            if obj.nNodes > 0
                fprintf('--- Nodal Extents ---\n');
                fprintf('First node position: %.4f cm\n', obj.nodes_vec(1));
                fprintf('Last node position:  %.4f cm\n', obj.nodes_vec(end));
            end
            fprintf('\n');
        end
    end
end