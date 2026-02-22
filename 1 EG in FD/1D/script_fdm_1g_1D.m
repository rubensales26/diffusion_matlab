clearvars;
close all; 
%%

% DATA
materials_vec = [1; 2; 1];
sizes_vec = [1; 2; 1];
D = [1.446; 0.776];
sigma_a = [0.0077; 0.0244];
nu_sigma_f = [0; 0.0260]; % Mat 1 has no fission

% Material grid and mesh definition
Materials = Materials1gD1(materials_vec, sizes_vec, D, sigma_a, nu_sigma_f);
Mesh = Mesh1DFD(Materials, 1);
Mesh.displayMesh;

% Solver initialization
Problem = DiffusionProblem_1D_1g_FDM(Mesh,Materials);
Problem = Problem.assembleMatrices();
Problem = Problem.solveEigenvalues(3);

Problem.displayProblem;