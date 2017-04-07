#include "includes.hpp"
#include "solver.hpp"

int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif

	H2 h2;

	const int max_level = 9;	// Max level for "C" field
	const int min_level = 4;	// Min level for "C" field

	const int vel_level = 4;	// Level for velocity field


	ScalarField2 c;
	ScalarField2 c_comp;
	VectorField2 v;

	c.createEntireLevel(max_level);
	c_comp.createEntireLevel(max_level);
	v.createEntireLevel(vel_level);

	// Set velocity to divergence free rotating center
	// Generate simple divergence free vector field {Y, -X}
	for (auto e : v) {
		H2::pos_t pos = h2.getPosition(e.first);
		pos[0] -= 0.5; // center around 0,0
		pos[1] -= 0.5;
		e.second.x = pos[1] * 2;
		e.second.y = -pos[0] * 2;
	}

	// Set C to a sphere and then coarsen
	Vec2 center_sph({0.5, 0.75});
	data_t r_sph = 0.2;
	for (auto e : c) {
		Vec2 pos(h2.getPosition(e.first));
		e.second = (pos - center_sph).length() > r_sph ? 0 : 1;
	}

	c_comp = 1. - c;  // complement

	// our velocity field has inflow-outflow from our box, for complement we have to "feed" 1 at the boundary
	c_comp.boundary[0] = c_comp.boundary[1] = c_comp.boundary[2] = c_comp.boundary[3] =
		[](ScalarField2 *self, coord_t c)->data_t {
			return 1;	// Dirichlet to 1
		};
	for (auto &bp : c_comp.boundary_propagate)	// Disable all boundary propagation
		bp = false;

	// make v stress-free at boundary
	v.boundary[0] = v.boundary[1] = v.boundary[2] = v.boundary[3] =
		[](VectorField2 *self, coord_t c)->Vec2 {
			return self->get(self->hcs.removeBoundary(c)); // Neumann BC, derivative == 0
		};


	write_pgm("c_init.pgm", c, max_level);
	write_pgm("cc_init.pgm", c_comp, max_level);

	Matrix<data_t, ScalarField2> M;
	Solver<data_t, ScalarField2> solver;

	data_t time = 0;
	data_t time_step = 0.01;

	// Implicit first-order upwind finite-volume stencil, no diffusion
	M.mul_stencil = [&v, &time_step](coord_t coord, data_t val, ScalarField2 &x)->data_t {
		level_t l = x.hcs.GetLevel(coord);
		data_t dist = 1. / ((coord_t)1 << x.hcs.GetLevel(coord));
		data_t face_area = pow(dist, x.hcs.GetDimensions() - 1);
		data_t vol = pow(dist, x.hcs.GetDimensions());

		data_t result = 0;
		Vec2 vel = v.get(coord);	// interpolated velocity for coord
		ScalarField2::coeff_map_t coeffs;
		for (int neighbor_direction = 0; neighbor_direction < x.hcs.parts; neighbor_direction++) {
			coord_t ne_coord = x.hcs.getNeighbor(coord, neighbor_direction);
			Vec2 vel_neighbor = v.get(ne_coord);
			Vec2 face_vel = (vel + vel_neighbor) * 0.5;	// average to get face-velocity
			Vec2 face_normal(x.hcs.getDirectionNormal(neighbor_direction));
			data_t flux = (face_vel * face_normal) * face_area / vol;

			ScalarField1::coeff_map_t coeffs_ne;
			x.getCoeffs(ne_coord, coeffs_ne, false);
			coeffs[coord] += flux > 0 ? flux : 0;
			for (const auto &e : coeffs_ne) {
				// if (e.first == coord) ... flux invert?
				coeffs[e.first] += flux < 0 ? flux * e.second : 0;
			}
		}
		coeffs[coord] += 1. / time_step;
		for (auto e : coeffs)
			result += x.get(e.first) * e.second;
		return result;
	};

	for (int step = 1; step < 1000; step++) {
		ScalarField2 rhs = c * 1. / time_step;  // or b-vector
		ScalarField2 rhs_comp = c_comp * 1. / time_step;  // or b-vector
		solver.solve(M, c, rhs, 1000, 1e-8, 1e-10);  // BiCGStab
		solver.solve(M, c_comp, rhs_comp, 1000, 1e-8, 1e-10);  // BiCGStab
		write_pgm("c_" + to_string(step) + ".pgm", c, max_level);
		write_pgm("cc_" + to_string(step) + ".pgm", c_comp, max_level);

		ScalarField2 c_diff = 1. - c_comp;
		c_diff -= c;
		for (auto e : c_diff)
			e.second = abs(e.second);

		write_pgm("cc_diff_" + to_string(step) + ".pgm", c_diff, max_level);

	}



}
