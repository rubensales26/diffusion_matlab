# README

Diffusion equation in 1D with 1 energy group solver using the Finite Difference Method.

## Structure

- **DiffusionProblem_1D_1g_FDM**: A class that handles the problem data and obtaining its solution.
- **Mesh1DFD**: A class that handles the mesh (discretization) of the 1D reactor (line segment).
- **Materials1gD1**: A class that handles the material properties of the regions in the reactor.
- **script_fdm_1g_1D.m**: Script with data that solves an instance of the problem.

---

# Classes

## Materials1gD1

## Mesh1DFD

The constructor takes as input a *Materials* object and a *refinment* parameter. The latter represents the number of nodes in the smallest material region of the reactor. In order to avoid possible problems such as a node falling in one of the material boundaries, we will generate a "pseudo-uniform" mesh, considering each material region and dividing it such that we have the following equality between the ratios:

$$nodes_i/length_i = refinment/minlength, \forall i material region$$.

Thus, to get the number of nodes for the material region i we just compute (and pass it through the ceil function just in case);

$$nodes_i = ceil(refinment length_i/minlength)$$.

## DiffusionProblem_1D_1g_FDM

Instead of using `for` loops to iterate through each node and populate the tridiagonal system element by element, we perform simultaneous operations on vectors mapped to the nodal structure.

#### 1. Nodal Property Mapping
First, the algorithm extracts the material properties and the step size for every individual node. Using the `region_idx` vector, it maps the cross-sections and diffusion coefficients from each material region down to the nodes:
$$D_i, \Sigma_{a,i}, \nu\Sigma_{f,i}, \Delta x_i \quad \text{for } i = 1, \dots, nNodes$$

#### 2. Coupling Coefficients
To define the leakage of neutrons between adjacent nodes $i$ and $i+1$, we calculate the coupling coefficient $d_i$. Based on the finite difference approximation with harmonic averaging for the interface between different materials and varying mesh sizes, the formula is:

$$d_i = \frac{2 D_i D_{i+1}}{\Delta x_{i+1} D_i + \Delta x_i D_{i+1}}$$

In the code, this is vectorized by creating two shifted arrays of length $nNodes-1$: a "current" array (indices $1$ to $nNodes-1$) and a "next" array (indices $2$ to $nNodes$). MATLAB performs the element-wise arithmetic (`.*`, `./`) on these arrays instantly to produce `d_vec`.



#### 3. Boundary Leakage Coefficients
Assuming a Zero Flux (vacuum) boundary condition at the left and right edges of the reactor, the leakage coefficients at the boundaries are computed as:

$$L_{left} = \frac{2 D_1}{\Delta x_1}$$

$$L_{right} = \frac{2 D_N}{\Delta x_N}$$

#### 4. The Loss/Transport Matrix (A)
The matrix $\mathbf{A}$ represents the balance of neutrons lost via absorption and transport. We build it as a sparse tridiagonal matrix.

**Diagonal Elements ($A_{i,i}$):** These represent the total removal of neutrons from node $i$. For an internal node, it is the sum of the absorption in the node volume and the leakage to its left and right neighbors:
$$A_{i,i} = \Delta x_i \Sigma_{a,i} + d_{i-1} + d_i$$

For the boundary nodes, the outward coupling is replaced by the boundary leakage $L$:
$$A_{1,1} = \Delta x_1 \Sigma_{a,1} + d_1 + L_{left}$$

$$A_{N,N} = \Delta x_N \Sigma_{a,N} + d_{N-1} + L_{right}$$

In MATLAB, shifting the `d_vec` array up and down with zeros allows us to add the left and right couplings to the entire diagonal vector in one operation.

**Off-Diagonal Elements ($A_{i,i+1}$ and $A_{i+1,i}$):**
These represent the neutrons entering node $i$ from its neighbors:
$$A_{i, i+1} = A_{i+1, i} = -d_i$$



Using `spdiags`, these vectors are efficiently mapped into the main diagonal (`0`) and the upper/lower diagonals (`1` and `-1`).

#### 5. The Production Matrix (Q)
The matrix $\mathbf{Q}$ represents the source of new neutrons generated from fission. Because fission neutrons produced in node $i$ are born in node $i$, this is purely a diagonal matrix where each element is the production rate scaled by the nodal volume:

$$Q_{i,i} = \Delta x_i \nu\Sigma_{f,i}$$

This is easily vectorized as `diag_Q = Delta_vec .* NuSigF_vec` and placed into a sparse diagonal matrix using `spdiags`.