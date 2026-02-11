% This class implements a mesh for a squared-shaped reactor for order 2 FDM

classdef Mesh2D
    properties
        L (1,1) double {mustBePositive} % Length of the side
        M (1,1) integer {mustBeInteger} % Number of nodes on each side
        cell_size (:,:) double % Matrix containing the length of the side of every square cell

        X (:,:)
        Y (:,:)
    end
    
    methods
        function m = Mesh2D(L, M, cell_size)
            m.L = L;
            m.M = M;
            
            if sum(cell_size,"all") ~= L^2 || all(size(cell_size)) ~= M
                error("Invalid size: cell size vector must have length L = %d for each row/column and its elements must sum M = %d", L, M);
            else
                m.cell_size = cell_size;
            end
            
            x_nodes = linspace(0, L, M+1);
            y_nodes = linspace(0, L, M+1);

            [m.X,m.Y] = meshgrid(x_nodes,y_nodes);

        end
    end
end
