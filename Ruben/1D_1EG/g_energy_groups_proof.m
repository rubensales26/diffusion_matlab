clearvars;
close all;

%%

% DATA
region_materials = [1;2;1];
region_lengths = [10; 20; 10]; % in cm
cells_per_region = [1; 2; 1];
D_lib = [1.446; 0.776]; % in cm
sigma_a_lib = [0.0077; 0.0244]; % in 1/cm
nu_sigma_f_lib = [0; 0.0260]; % Mat 1 has no fission; 1/cm

% Material grid and mesh definition
mesh = Mesh_1D_FDM(region_lengths, cells_per_region);
mesh.displayMesh;

materials = Materials(region_materials, D_lib, sigma_a_lib, nu_sigma_f_lib);
materials.displayMaterials

solver = Solver_1D_FDM(mesh, materials);
solver = solver.assembleMatrices;
disp(solver.A)
disp()
%solver.displayProblem;
%solver.plotPhi;