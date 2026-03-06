clearvars;
close all;

%%

% DATA
region_materials = [1;2;1];
region_lengths = [10; 20; 10]; % in cm
nodes_per_region = [1; 2; 1];
D_lib = [1.446; 0.776]; % in cm
sigma_a_lib = [0.0077; 0.0244]; % in 1/cm
nu_sigma_f_lib = [0; 0.0260]; % Mat 1 has no fission; 1/cm

% Material grid and mesh definition
mesh = Mesh_1D_FDM(region_lengths, region_materials, nodes_per_region);
mesh.displayMesh;

materials = Materials_1D_1eg_FDM(region_materials, D_lib, sigma_a_lib, nu_sigma_f_lib);
materials.displayMaterials

problem = Problem_1D_1eg_FDM(mesh, materials);
problem = problem.assembleMatrices;
problem = problem.solveEigenvalues(1);
problem.displayProblem;
problem.plotPhi;
%%

c_mesh = clean_mesh(region_lengths, nodes_per_region);
c_materials = clean_materials(region_materials, D_lib, sigma_a_lib, nu_sigma_f_lib)

c_problem = clean_solver(c_mesh, c_materials);

c_problem = c_problem.assembleMatrices.solveEigenvalues(1);
c_problem.displayProblem;
c_problem.plotPhi;
