#include "includes.hpp"
#include "solver.hpp"

int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif

	H1 h1;
	H2 h2;
	H3 h3;

	const int max_level = 5;

	ScalarField1 x;
	x.createEntireLevel(max_level);
	cout << "Solving Poisson equation Nabla^2(x) = 1 with zero Dirichlet boundary condition in 1D.\n";

	Solver<data_t, ScalarField1> solver;
	Matrix<data_t, ScalarField1> M;

	// The Laplacian stencil, a matrix-free implementation. Works at all levels.
	M.setStencil([](coord_t coord, ScalarField1& x)->ScalarField1::coeff_map_t {
		level_t l = x.hcs.GetLevel(coord);
		coord_t level = (coord_t)1 << x.hcs.GetLevel(coord);
		data_t dist = 1. / level;  // Neighbor distance, assuming all directions with equal scale and scale = 1
		data_t vol = dist;
		ScalarField1::coeff_map_t coeffs;
		for (int ne_idx = 0; ne_idx < x.hcs.parts; ne_idx++) {
			coord_t ne_coord = x.hcs.getNeighbor(coord, ne_idx);
			data_t coeff  = (1. ) / (dist);
			if (x.hcs.IsBoundary(ne_coord))
				coeff *= 2;  // If we hit a boundary, its only half the distance
			ScalarField1::coeff_map_t coeffs_ne;
			x.getCoeffs(ne_coord, coeffs_ne, false);
			coeffs[coord] -= coeff / vol;
			for (const auto &e : coeffs_ne) {
				coord_t level_ne = x.hcs.GetLevel(e.first);
				coeffs[e.first] += (coeff / vol) * e.second;
			}
		}
		return coeffs;
	});

	x = 0;
	coord_t hole = h1.createFromPosition(4,{0.49});
	coord_t hole2 = h1.createFromPosition(3,{0.24});
	x.coarse(hole);
	x.coarse(hole2);
	x.refineTo(h1.createFromPosition(9,{0.51}));
	ScalarField1 b = x;
	b = 1;


	for (auto e : b) if (!b.isTop(e.first)) e.second = 0;

	data_t tot = 0;
	for (auto e : b) if (b.isTop(e.first)) {
		level_t l = h1.GetLevel(e.first);
		tot += e.second * 1./((coord_t)1 << l);
		//e.second *= ();
	}
	cout << tot << endl;
	//b[h1.createFromPosition(max_level, {0.88})] = 1000; // creates a spike at the highest level in b
	//b.propagate();  // always a good idea. preserves integral of b at lower levels


	// Turns the right-side boundary value to 1

	x.boundary[0] = [](ScalarField1 *self, coord_t c)->data_t {
		//return self->get(self->hcs.removeBoundary(c)); // Neumann, derivative == 0
		return 1;	// Dirichlet to 1
	};
	x.boundary_propagate[0] = false; // propagation must be disabled for Dirichlet != 0 or solver will fail. For Neumann, propagation must stay enabled.

//	M.mul(b, x);

	cout << "Finished. Iterations: " << solver.solve(M, x, b, 10000, 1e-12, 1e-20) << "\n";

	ScalarField1 res = b;
	M.mul(x, res);
	res = b - res;
	//x=res;

	map<data_t, data_t> x_val_map;
	for (auto e : x) {
		if (!x.isTop(e.first))
			continue;
		data_t x_pos = h1.getPosition(e.first)[0];
		x_val_map[x_pos] = e.second;
		//cout << h1.toString(e.first) << " V: " << e.second << endl;
	}

	ofstream out("test6.txt");
	for (auto e : x_val_map)
		out << e.first << " " << e.second << endl;
	out.close();
	/*
	for (int o_level = 0; o_level <= max_level; o_level++) {
		ofstream out("test6_level_" + to_string(o_level) + ".txt");
		out << "0 " << x.get(h1.getNeighbor(0, 1)) << endl; // Left Boundary value
		for (data_t ux = 0; ux < pow(2, o_level); ux++) {
			coord_t c = h1.createFromUnscaled(o_level, {ux});
			out << h1.getPosition(c)[0] << " " << x.get(c) << endl;
		}
		out << "1 " << x.get(h1.getNeighbor(0, 0)) << endl; // Right Boundary value
		out.close();
	}
*/
}
