#include "includes.hpp"




int main(int argc, char **argv) {

	H2 h2;  // The 2D coordinate system

	ScalarField2 x2;
	VectorField2 v2;

	v2.createEntireLevel(9);  // 2D box with 256x256 vectors
	x2.createEntireLevel(9);
	// Generate simple divergence free vector field {Y, -X}
	for (auto e : v2) { 
		H2::pos_t pos = h2.getPosition(e.first);
		pos[0] -= 0.5; // center around 0,0
		pos[1] -= 0.5;
		e.second.x = pos[1] * 200;
		e.second.y = -pos[0] * 200;
	}
	// Optional play with refinement / coarsening
	//v2.refineTo(h2.createFromList({0,3,1,1,1,1,1,1,1}));
	v2.propagate(); // averages all lower-level coords
	//v2.coarse(h2.createFromList({1,2,1}));
	//v2.propagate();
	cout << v2.getHighestLevel() << endl;
	
	auto t1 = high_resolution_clock::now();
	for (int i = 0; i < 100; i++) {
		//v2[h2.createFromUnscaled(8, {127, 127})] = Vec2({100,i}); // create some divergence in the center
		
		// Use convert to calculate explicit divergence
		//x2.takeStructure<Vec2>(v2);

		x2.convert<Vec2>(v2, [&h2, &v2](coord_t c, Vec2& t2)->data_t {
			// Divergence calculation with finite-volume stencil

			level_t l = h2.GetLevel(c);
			data_t dist = 2 * (h2.scales[0] / data_t(1U << l)); // neighbor distance at that level, assuming all scales equal
			data_t vol = pow(dist, 2);  // dist is also face area and therefore dist*dist is the volume of the cell
			data_t divergence = 0;
			for (int n_idx = 0; n_idx < h2.parts; n_idx++) {  // Traverse all 2^D neighbors
				coord_t c_ne = h2.getNeighbor(c, n_idx);
				if (h2.IsBoundary(c_ne))
					continue;	// assume zero flux at BC
				Vec2 ne_val = v2.getDirect(c_ne); // convert() ensures that c exists in v2.
				Vec2 normal = Vec2(v2.hcs.getDirectionNormal(n_idx));
				Vec2 interp = (t2  + ne_val) * 0.5; // The vector at the face
				divergence += (interp * normal) * (dist / vol);
			}
			return divergence;
		});
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "100 Divergence calculations with 256^2 vec took " << duration << "ms. \n";

	// Write divergence as grey-scale image
	x2.propagate();
	write_pgm("test4.pgm", x2, 9);

	// turn vector field into scalar field containing only x (or u) component
	x2.convert<Vec2>(v2, [](coord_t c, Vec2& t2)->data_t {
		return t2.x;
	});
	write_pgm("test4x.pgm", x2, 8);

	x2.convert<Vec2>(v2, [](coord_t c, Vec2& t2)->data_t {
		return t2.y;
	});
	write_pgm("test4y.pgm", x2, 8);

}
