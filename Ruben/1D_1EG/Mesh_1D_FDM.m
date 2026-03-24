classdef Mesh_1D_FDM
    %CLEAN_MESH Summary of this class goes here
    %   Detailed explanation goes here
    properties
        % Input attributes
        region_lengths (:,1) double       % Length of each physical region
        cells_per_region (:,1) double     % Number of discrete cells in each region
        
        % Computed attributes
        nCells (1,1) double               % Total number of cells in the slab
        cell_centers (:,1) double         % Cartesian coordinates of the cell centers
        cell_region_idx (:,1) double      % The physical region index each cell belongs to
        delta_x (:,1) double              % The step size (dx) for every single cell
        region_boundaries (:,1) double    % Cartesian coordinates of the region interfaces
    end
    
    methods
        function obj = Mesh_1D_FDM(region_lengths, cells_per_region)
            %Mesh_1D_FDM Constructor of the class
            
            % Check if the input leads to a valid mesh
            if length(region_lengths) ~= length(cells_per_region)
                error("Invalid input: Arguments region_lengths and cells_per_region must have the same dimensions.")
            end
            
            % Store inputs as column vectors
            obj.region_lengths = region_lengths;
            obj.cells_per_region = cells_per_region;
            
            % Compute region boundaries and total cells
            obj.region_boundaries = [0; cumsum(obj.region_lengths)];
            obj.nCells = sum(obj.cells_per_region);
            
            % Pre-allocate arrays
            obj.cell_centers = zeros(obj.nCells, 1);
            obj.cell_region_idx = zeros(obj.nCells, 1);
            obj.delta_x = zeros(obj.nCells, 1); 
            
            % Populate the cells
            ind = 1; 
            for i = 1:length(obj.region_lengths)
                % Step size for a given region
                dx = obj.region_lengths(i) / obj.cells_per_region(i);
                left_boundary = obj.region_boundaries(i);
                
                for j = 1:obj.cells_per_region(i)
                    obj.delta_x(ind) = dx;
                    % Place the calculation point in the middle of the cell
                    obj.cell_centers(ind) = left_boundary + dx * (j - 0.5);
                    obj.cell_region_idx(ind) = i; % Tag cell with its physical region
                    ind = ind + 1;
                end
            end
        end
        
        function displayMesh(obj)
            %displayMesh Display on screen mesh information
            
            fprintf('==========================\n');
            fprintf('        MESH DATA         \n');
            fprintf('==========================\n\n');
            
            fprintf('Total Slab Length: %.2f cm\n', obj.region_boundaries(end));
            fprintf('Total Cells:       %d\n', obj.nCells);
            fprintf('Number of Regions: %d\n\n', length(obj.region_lengths));
            
            % Display Table of Regions
            fprintf('--- Region Breakdown ---\n');
            Region = (1:length(obj.region_lengths))';
            Length = obj.region_lengths;
            Cells = obj.cells_per_region;
            
            T = table(Region, Length, Cells);
            disp(T);
        end
    end
end

