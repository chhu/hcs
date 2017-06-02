#include "includes.hpp"
#include "solver.hpp"



ScalarField2 criteria(ScalarField2 &f) {
	ScalarField2 result;
	VectorField2 grad_f;
	grad_f.takeStructure(f);
	result.takeStructure(f);
	grad<2>(f, grad_f);
	div<2>(grad_f, result);
	result.convert<data_t>(result, [](coord_t c, ScalarField2 &f)->data_t{
		return pow(f.get(c) * 1./(512.*512.), 2);//1/pow(1U << f.hcs.GetLevel(c), 2);
	});
	result.propagate();
	return result;
}

void refinement(ScalarField2 &f, ScalarField2 &criteria, data_t sensitivity, level_t lowest_level, level_t highest_level, coord_t start = 0) {
	H2 &h = f.hcs;
	level_t current = h.GetLevel(start);
	data_t critical = criteria.get(start, true) * sensitivity;
//	cout << h.toString(start) << " CRIT: " << critical << endl;
	bool keep = critical >= 1 || current < lowest_level;
	if (keep) {
		if (!f.exists(start))
			f.refineTo(start);
		if (current == highest_level)
			return;
		for (uint16_t i = 0; i < h.parts; i++) {
			coord_t next = h.IncreaseLevel(start, i);
			refinement(f, criteria, sensitivity, lowest_level, highest_level, next);
		}
	} else {
		f.coarse(start);
	}

}

int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif



	// Dimension setup
	typedef H2 HCS;
	typedef ScalarField2 ScalarField;
	typedef VectorField2 VectorField;
	typedef Vec2 Vec;

	// Resolution setup
	const int max_level = 9;	// Max level for "C" field
	const int min_level = 4;	// Min level for "C" field
	const int vel_level = 4;	// Level for velocity field


	HCS hcs;

	ScalarField c;
	VectorField v;

	//cout << Vec2(hcs.getDirectionNormal(0)) << endl;return 0;

	c.createEntireLevel(max_level);
	v.createEntireLevel(vel_level);


	// Set velocity to divergence free rotating center
	// Generate simple divergence free vector field {Y, -X}
	for (auto e : v) {
		HCS::pos_t pos = hcs.getPosition(e.first);
		pos[0] -= 0.5; // center around 0,0
		pos[1] -= 0.5;
		e.second.x = pos[1] * 2;
		e.second.y = -pos[0] * 2;
	}
	v.propagate();


	// Set C to a sphere and then coarsen
	Vec center_sph({0.5, 0.75}); 	// sphere center
	data_t r_sph = 0.2;				// ... and radius
	for (auto e : c) {
		Vec pos(hcs.getPosition(e.first));
		e.second = (pos - center_sph).length() > r_sph ? 0 : 1;
	}
	c.propagate();

/*
	int n_refinements = 10;
	// Randomly coarse into a coordinate.
	// Create random structure between level 5 and 10 coords.
	for (int i = 0; i < n_refinements; i++) {
		double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		double _l = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		level_t l = double(_l * (max_level - min_level) + min_level);
		coord_t c2 = hcs.createFromPosition(l, {x, y});
		c.coarse(c2);
	}
*/
	// make v stress-free at all boundaries
	for (auto &boundary : v.boundary)
		boundary = [](VectorField *self, coord_t c)->Vec {
			return self->get(self->hcs.removeBoundary(c)); // Neumann BC, derivative == 0
		};

	data_t sense = 256;

	write_pgm("c_init.pgm", c, max_level);


	ScalarField2 gm = criteria(c);
	//write_pgm("crit_init.pgm", gm, max_level);
	gm.propagate(true);

	refinement(c, gm, sense, min_level, max_level);
	//write_pgm_level("ar.pgm", c);
	c.propagate();
	//write_pgm("c1.pgm", c, max_level);



	Matrix<data_t, ScalarField> M;
	Solver<data_t, ScalarField> solver;

	data_t time = 0;
	data_t time_step = 0.001;

	// Implicit first-order upwind finite-volume stencil, no diffusion
	M.setStencil([&v, &time_step](coord_t coord, ScalarField &x)->ScalarField::coeff_map_t {
		level_t l = x.hcs.GetLevel(coord);
		data_t dist = 1. / ((coord_t)1 << x.hcs.GetLevel(coord));	// distance to neighbors
		data_t face_area = pow(dist, x.hcs.GetDimensions() - 1);	// all face area around our box are equal
		data_t vol = pow(dist, x.hcs.GetDimensions());

		data_t result = 0;
		Vec vel = v.get(coord);	// interpolated velocity for coord
		ScalarField::coeff_map_t coeffs;
		for (int neighbor_direction = 0; neighbor_direction < x.hcs.parts; neighbor_direction++) {
			coord_t ne_coord = x.hcs.getNeighbor(coord, neighbor_direction);
			Vec vel_neighbor = v.get(ne_coord);
			Vec face_vel = (vel + vel_neighbor) * 0.5;	// average to get face-velocity
			Vec face_normal(x.hcs.getDirectionNormal(neighbor_direction));
			data_t flux = (face_vel * face_normal) * face_area / vol;

			ScalarField::coeff_map_t coeffs_ne;
			x.getCoeffs(ne_coord, coeffs_ne, false);

			if (flux > 0)
				coeffs[coord] += flux;
//			coeffs[coord] += flux * 0.5;
			for (const auto &e : coeffs_ne) {
				if (flux < 0)
					coeffs[e.first] += flux * e.second;
			}

		}

		coeffs[coord] += 1. / time_step;
		return coeffs;
	});

	// Benchmark mat-mul
	ScalarField a = c;
	auto t1 = high_resolution_clock::now();
	for (int i = 0; i < 10; i++) {
		M.mul(c, a);
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "10 mul " << " took " << duration << "ms.\n";

	// Run 1000 time steps
	for (int step = 1; step < 1000; step++) {
		auto t1 = high_resolution_clock::now();
		M.clearCache();
		ScalarField rhs = c * 1. / time_step;  // right-hand-side or b-vector. Solve M * c = rhs
		int its = solver.solve(M, c, rhs, 1000, 1e-8, 1e-14);  // ... with BiCGStab
		c.propagate();
		write_pgm("c_" + to_string(step) + ".pgm", c, max_level);
		ScalarField2 gm = criteria(c);
		write_pgm("p_" + to_string(step) + ".pgm", gm, max_level);
		gm.propagate(true);
		refinement(c, gm, sense, min_level, max_level);
		write_pgm_level("l_" + to_string(step) + ".pgm", c);
		c.propagate();

		auto t2 = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(t2-t1).count();
		cout << "Time-step " << step << " took " << its << " iterations and " << duration << "ms. Elements: " << c.nElementsTop() << endl;
	}

}
