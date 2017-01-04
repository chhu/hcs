#include "includes.hpp"




int main(int argc, char **argv) {


	// Test suite
	H2 h2;


	ScalarField2 x2('x', &h2);
	VectorField2 v2('v', &h2);

	v2.createEntireLevel(8);

	cout << v2.getHighestLevel() << endl;
	for (auto e : v2) {
		H2::pos_t pos = h2.getPosition(e.first);
		pos[0] -= 0.5;
		pos[1] -= 0.5;
		e.second.x = pos[1] * 200;
		e.second.y = -pos[0] * 200;
	}
	//v2.refineTo(h2.createFromList({0,3,1,1,1,1,1,1,1}));
	v2.propagate();
//	v2.coarse(h2.createFromList({1,2,1}));
//	v2.propagate();
	cout << v2.getHighestLevel() << endl;
	auto t1 = high_resolution_clock::now();
	for (int i = 0; i < 100; i++) {
		v2[h2.createFromUnscaled(8, {127, 127})] = Vec2({100,i});
		x2.convert<Vec2>(v2, [&h2, &v2](coord_t c, Vec2& t2)->data_t {
			// Explicit divergence calculation with finite-volume stencil
			level_t l = h2.GetLevel(c);
			data_t dist = 2 * (h2.scales[0] / pow(2, l)); // neighbor distance at that level, assuming all scales equal
			data_t vol = dist * dist;  // dist is also face area and therefore dist*dist is the volume of the cell
			data_t divergence = 0;
			for (int n_idx = 0; n_idx < h2.parts; n_idx++) {  // Traverse all neighbors
				coord_t c_ne = h2.getNeighbor(c, n_idx);
				if (h2.IsBoundary(c_ne))
					return 0;	// assume zero flux at BC
				Vec2 ne = v2.get(c_ne);
				Vec2 normal;
				switch (n_idx) { // there is probably a better way...
				case 0: normal = {1 , 0 };break;
				case 1: normal = {-1, 0 };break;
				case 2: normal = {0 , 1 };break;
				case 3: normal = {0 ,-1 };break;
				}
				Vec2 interp = t2 * 0.5 + ne * 0.5;
				divergence += (interp * normal) * dist / vol;
			}
			return divergence;
		});
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "100 Divergence calculations with 256^2 vec took " << duration << "ms. \n";


	x2.propagate();
	write_pgm("test4.pgm", x2, 8);

	x2.convert<Vec2>(v2, [](coord_t c, Vec2& t2)->data_t {
		return t2.x;
	});
	write_pgm("test4x.pgm", x2, 8);

	x2.convert<Vec2>(v2, [](coord_t c, Vec2& t2)->data_t {
		return t2.y;
	});
	write_pgm("test4y.pgm", x2, 8);

}
