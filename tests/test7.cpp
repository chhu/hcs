#include "includes.hpp"
#include "solver.hpp"

int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif



	// Dimension setup
	typedef H2 HCS;
    typedef DenseScalarField2 ScalarField;
    typedef DenseVectorField2 VectorField;
    typedef ScalarField2 ScalarFieldBase;
    typedef VectorField2 VectorFieldBase;
	typedef Vec2 Vec;

	// Resolution setup
	const int c_level = 8;	// Max level for "C" field
	const int vel_level = 4;	// Level for velocity field


	HCS hcs;

	ScalarField c(c_level);
	VectorField v(vel_level);

	// Staging: Set initial conditions for c & v
	if (HCS::GetDimensions() == 1) {
		// Generate simple v = 1
		v = Vec(1);
		for (auto e : c) {
			Vec pos_(hcs.getPosition(e.first));
			data_t pos = pos_[0];
			e.second = pos > 0.2 && pos < 0.3  ? 1 : 0;
		}
	}

	//coord_t coarsed_coord = hcs.createFromPosition(max_level-4,{0.75,0.5});

	if (HCS::GetDimensions() == 2) {
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
		//c.coarse(coarsed_coord);
		c.propagate();
//		c[coarsed_coord] = 0;

	}

	// make v stress-free at all boundaries
	for (auto &boundary : v.boundary)
		boundary = [](VectorFieldBase *self, coord_t c)->Vec {
			return self->get(self->hcs.removeBoundary(c)); // Neumann BC, derivative == 0
		};

	write_pgm("c_init.pgm", c, c_level);
//return 0;
	Matrix<data_t, ScalarField> M;
	Solver<data_t, ScalarField> solver;

	data_t time = 0;
	data_t time_step = 0.01;

	// Implicit first-order upwind finite-volume stencil, no diffusion
	M.setStencil([&v, &time_step](coord_t coord, ScalarFieldBase &x)->ScalarFieldBase::coeff_map_t {
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
	for (int i = 0; i < 100; i++) {
		M.mul(c, a);
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "100 mul " << " took " << duration << "ms.\n";

	// Run 1000 time steps
	for (int step = 1; step < 1000; step++) {
		auto t1 = high_resolution_clock::now();

		ScalarField rhs = c * 1. / time_step;  // right-hand-side or b-vector. Solve M * c = rhs
		int its = solver.solve(M, c, rhs, 1000, 1e-8, 1e-14);  // ... with BiCGStab
		c.propagate();
		write_pgm("c_" + to_string(step) + ".pgm", c, c_level);

		auto t2 = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(t2-t1).count();
		cout << "Time-step " << step << " took " << its << " iterations and " << duration << "ms.\n";
	}

}
