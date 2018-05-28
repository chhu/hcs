#include "includes.hpp"




int main(int argc, char **argv) {

	H2 h2;  // The 2D coordinate system


	level_t max_level = 9;
	DenseScalarField2 x2(max_level);
	DenseVectorField2 v2(max_level);

	// Generate simple divergence free vector field {Y, -X}
	for (auto e : v2) { 
		H2::pos_t pos = h2.getPosition(e.first);
		pos[0] -= 0.5; // center around 0,0
		pos[1] -= 0.5;
		e.second.x() = pos[1] * 200;
		e.second.y() = -pos[0] * 200;
	}
	// Optional play with refinement / coarsening
	//v2.refineTo(h2.createFromList({0,3,1,1,1,1,1,1,1}));
	v2.propagate(); // averages all lower-level coords
	//v2.coarse(h2.createFromList({1,2,3}));

	cout << v2.getHighestLevel() << endl;
	
	auto t1 = high_resolution_clock::now();
	for (int i = 0; i < 100; i++) {
		div<2>(v2, x2);
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "100 Divergence calculations with " << (1U << v2.getHighestLevel()) << "^2 vec took " << duration << "ms. \n";

	x2.propagate();

	write_pgm_level("test4l.pgm", x2);

	// Write divergence as grey-scale image
	write_pgm("test4.pgm", x2, max_level);

	// turn vector field into scalar field containing only x (or u) component
	x2.convert<Vec2>(v2, [](coord_t c, VectorField2 &v)->data_t {
		return v.get(c).x();
	});
	x2.propagate();	// convert() operates on TLCs only
	write_pgm("test4x.pgm", x2, max_level);

	x2.convert<Vec2>(v2, [](coord_t c, VectorField2 &v)->data_t {
		return v[c].y();
	});
	x2.propagate();
	write_pgm("test4y.pgm", x2, max_level);

}
