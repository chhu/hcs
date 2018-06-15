
#include "includes.hpp"
#include "solver.hpp"



int main(int argc, char **argv) {

	H2 h2;

	const int max_level = 9;
	SparseScalarField2 x(max_level);

	cout << "Solving Poisson equation Nabla^2(x) = 1 with zero Dirichlet boundary condition in 1D.\n";

	Solver<data_t, SparseScalarField2> solver;
	Matrix<data_t, SparseScalarField2> M;

	// The Laplacian stencil, a matrix-free implementation. Works at all levels.
	M.setStencil([](coord_t coord, ScalarField2& x)->ScalarField2::coeff_map_t {
		level_t l = x.hcs.GetLevel(coord);
		coord_t level = (coord_t)1 << x.hcs.GetLevel(coord);
		data_t dist = 1. / level;  // Neighbor distance, assuming all directions with equal scale and scale = 1
		data_t vol = dist * dist;
		ScalarField2::coeff_map_t coeffs;
		for (int ne_idx = 0; ne_idx < x.hcs.parts; ne_idx++) {
			coord_t ne_coord = x.hcs.getNeighbor(coord, ne_idx);
			data_t coeff  = (1. ) / (vol);
			if (x.hcs.IsBoundary(ne_coord))
				coeff *= 2;  // If we hit a boundary, its only half the distance
			ScalarField2::coeff_map_t coeffs_ne;
			x.getCoeffs(ne_coord, coeffs_ne, false);
			coeffs[coord] -= coeff;
			for (const auto &e : coeffs_ne) {
				coord_t level_ne = x.hcs.GetLevel(e.first);
				coeffs[e.first] += coeff * e.second;
			}
		}
		return coeffs;
	});

	x = 0;
	SparseScalarField2 b = x;
	b = 1;
	/*
	b = 0;
    b[h2.createFromUnscaled(9, {256, 256})] = 1;//(512 *512);
    b[h2.createFromUnscaled(9, {255, 255})] = 1;//(512 *512);
    b[h2.createFromUnscaled(9, {256, 255})] = 1;//(512 *512);
    b[h2.createFromUnscaled(9, {255, 256})] = 1;//(512 *512);
	 */
	// Turns the right-side boundary value to 1

	x.boundary[0] = [](ScalarField2 *self, coord_t c)->data_t {
		//return self->get(self->hcs.removeBoundary(c)); // Neumann, derivative == 0
		return 0;	// Dirichlet to 1
	};
	x.boundary_propagate[0] = false; // propagation must be disabled for Dirichlet != 0 or solver will fail. For Neumann, propagation must stay enabled.

	cout << "Finished. Iterations: " << solver.solve(M, x, b, 10000, 1e-5, 1e-20) << "\n";
	x.write("test6.raw");
}
