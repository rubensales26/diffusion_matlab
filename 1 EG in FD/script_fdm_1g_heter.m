% DESCRIPTION:
% Main script for solving the difusion equation for 1 energy group with
% heterogeneus materials.
% It runs on FDM.
% Configuration "Reflector-Fuel-Reflector":
clearvars;
close all; 
%%
L = 75;
N_nodes_side = 25;


% Fisical properties (Mat 1: Reflector, Mat 2: Fuel)
D = [1.446, 0.776];
sigma_a = [0.0077, 0.0244];
nu_sigma_f = [0, 0.0260]; % Mat 1 has no fission

% Mesh definition and material distribution
grid_materials = [1, 1, 1; 1, 2, 1; 1, 1, 1];
cell_size = 25;

mesh = Mesh2D(L, N_nodes_side, N_nodes_side);

materials = Materials2D1g(mesh, grid_materials, cell_size, D, sigma_a, nu_sigma_f);

[ind, D, sigma_a, nu_sigma_f] = materials.GetProperties(25,25);

problem = DifusionProblem_2D_1g_FD_O2(mesh, materials);

problem = problem.Solve();

fprintf('k_eff = %.5f\n', problem.k_eff);

% Visualize the Flux
figure('Name', 'Reactor Flux Results');

% Plot 1: The Flux Map
subplot(1, 2, 1);
imagesc([0 L], [0 L], problem.Phi);
set(gca, 'YDir', 'normal'); % Flip Y so (0,0) is bottom-left
colormap('jet');
colorbar;
title(['Neutron Flux (\phi), k_{eff} = ' num2str(problem.k_eff, '%.5f')]);
xlabel('X (cm)'); ylabel('Y (cm)');
axis square;

% Plot 2: A slice through the center (y = L/2)
subplot(1, 2, 2);
mid_index = round(mesh.ny / 2);
x_axis = linspace(0, L, mesh.nx);
plot(x_axis, problem.Phi(mid_index, :), 'LineWidth', 2);
grid on;
title('Flux Cross-Section (y = 37.5 cm)');
xlabel('X (cm)'); ylabel('Flux');
xlim([0 L]);