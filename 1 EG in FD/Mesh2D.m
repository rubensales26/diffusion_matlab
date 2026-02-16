% This class implements a mesh for a squared-shaped reactor for FDM
% The points in the mesh are centered inside of the cells
% using a 5 points stencil

classdef Mesh2D
    properties
        L (1,1) double % Length of the side
        nx (1,1) % Number of nodes on side x
        ny (1,1) % Number of nodes on side y
        hx (1,1) double % Mesh step size x
        hy (1,1) double % Mesh step size y
        nNodes (1,1) double

        X (:,:)
        Y (:,:)
    end
    
    methods
        function obj = Mesh2D(L, nCellsX, nCellsY)
            obj.L = L;
            obj.nx = nCellsX;
            obj.ny = nCellsY;
            obj.hx = L/obj.nx;
            obj.hy = L/obj.ny;
            obj.nNodes = obj.nx * obj.ny;
            
            x_c = (obj.hx/2) : obj.hx : (L - obj.hx/2);
            y_c = (obj.hy/2) : obj.hy : (L - obj.hy/2);

            [obj.X, obj.Y] = meshgrid(x_c, y_c);
        end

        function k = getIndexInVector(obj, i, j)
            % Get the lexicographical order of point (i,j)
            k = i + (j-1)*obj.ny;
        end

        function [x,y] = getPlaneCoordinates(obj, i, j)
            x = obj.X(i,j);
            y = obj.Y(i,j);
        end
    end
end
