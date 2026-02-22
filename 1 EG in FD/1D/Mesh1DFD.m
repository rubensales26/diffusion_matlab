% This class implements a mesh for a squared-shaped reactor for FDM
classdef Mesh1DFD
    properties
        % Input attributes
        Materials Materials1gD1 % Representation of the grid of materials
        refinment (1,1) double {mustBePositive, mustBeInteger} = 1% Number of nodes in the smallest material region

        % Computed attributes
        nNodes (1,1) double % Number of nodes in the mesh
        nodes_vec (:,1) double % Cartesian coordinates of the nodes of the mesh
        material_region_indices (:,1) double % Vector with the material region index of each node in the mesh
        delta_region (:,1) double % Vector with the step sizes for every region in the mesh
    end
    
    methods
        function obj = Mesh1DFD(Materials, refinment)
            obj.Materials = Materials;
            obj.refinment = refinment;
            materials_vec = obj.Materials.materials_vec;
            material_boundaries = obj.Materials.material_boundaries;
            material_sizes = obj.Materials.sizes_vec;

            % We will generate a mesh taking the material grid as a basis
            % In the smallest material cell we will create #refinment nodes
            % In the other cells we will keep the ratio #refinment/length
            % of the smallest cell
            [min_length, ~] = min(material_sizes);
            nRegions = length(materials_vec);

            % Calculate how many nodes for every material cell beforehand
            % (in order to know the length of the nodes_vec 
            % and material_region_indices arrays before defining them)
            % and the corresponding delta_region
            nNodes_material_cell = zeros(nRegions, 1);

            for i = 1:nRegions
                length_i = material_sizes(i);
                ref_i = ceil(refinment * length_i / min_length);
                nNodes_material_cell(i) = ref_i;
                obj.delta_region(i) = (material_boundaries(i+1) - material_boundaries(i)) / (ref_i + 1);
            end

            % Total nodes in the whole system
            obj.nNodes = sum(nNodes_material_cell);

            % Pre-allocate exactly the right size of the arrays
            obj.nodes_vec = zeros(obj.nNodes, 1);
            obj.material_region_indices = zeros(obj.nNodes, 1);
            
            % For every material region, for every node we want in it,
            % we compute the coordinates of the node and store the index of
            % which material region it is
            ind = 1; 
            for i = 1:nRegions
                left_boundary = material_boundaries(i);
                for j = 1:nNodes_material_cell(i)
                    obj.nodes_vec(ind) = left_boundary + obj.delta_region(i) * j;
                    obj.material_region_indices(ind) = i; % Stores which region this node belongs to
                    ind = ind + 1;
                end
            end
        end

        function displayMesh(obj)
            fprintf('=================\n    MESH DATA    \n=================\n\n')
            fprintf('Number of nodes: %d\n\n', obj.nNodes)
            fprintf('Coordinates of the nodes:\n')
            disp(obj.nodes_vec)
            fprintf('Material region indices of the nodes:\n')
            disp(obj.material_region_indices)
        end
    end
end
