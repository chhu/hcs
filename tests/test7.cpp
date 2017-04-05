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
	VectorField2 v;

	c.createEntireLevel(max_level);
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

	write_pgm("c_init.pgm", c, max_level);

	Matrix<data_t, ScalarField2> M;

	// Implicit first-order upwind finite-volume stencil, no diffusion
	M.mul_stencil = [&v](coord_t coord, data_t val, ScalarField2 &x)->data_t {
		level_t l = x.hcs.GetLevel(coord);
		data_t dist = 1. / ((coord_t)1 << x.hcs.GetLevel(coord));
		data_t face_area = pow(dist, x.hcs.GetDimensions() - 1);
		data_t vol = pow(dist, x.hcs.GetDimensions());

		data_t result = 0;
		Vec2 vel = v.get(coord);	// interpolated velocity for coord
		ScalarField2::coeff_map_t coeffs;
		for (int ne_idx = 0; ne_idx < x.hcs.parts; ne_idx++) {
			coord_t ne_coord = x.hcs.getNeighbor(coord, ne_idx);
			Vec2 vel_ne = v.get(ne_coord);
			Vec2 face_vel = (vel + vel_ne) * 0.5;
			Vec2 face_normal(x.hcs.getDirectionNormal(ne_idx));
			data_t flux = (face_vel * face_normal) * face_area / vol;

			ScalarField1::coeff_map_t coeffs_ne;
			x.getCoeffs(ne_coord, coeffs_ne, false);
			coeffs[coord] -= flux > 0 ? flux : 0;
			for (const auto &e : coeffs_ne) {
				// if (e.first == coord) ... flux invert?
				coeffs[e.first] += flux < 0 ? flux * e.second : 0;
			}
		}
		for (auto e : coeffs)
			result += x.get(e.first) * e.second;
		return result;
	};



}
