#include "includes.hpp"

int main(int argc, char **argv) {

	int n_refinements = 100000;

	H3 h3;
	H2 h2;
	H1 h1;

//	cout << __builtin_clzll(1) << " " << __builtin_clzll(0) << endl;
//	cout << __count_leading_zeros(1) << " " << __count_leading_zeros(0) << endl;
//	return 0;
//	cout << h2.toString(h2.getNeighbor(1, 0)) << endl;
//	cout << h2.toString(h2.getNeighbor(1, 1)) << endl;
//	cout << h2.toString(h2.getNeighbor(1, 2)) << endl;
//	cout << h3.toString(h3.getNeighbor(1, 5)) << endl;
//	return 0;
	ScalarField3 f;
	ScalarField2 u;
	/*
	ScalarField1 x;
	x.createEntireLevel(2);
	for (auto e = x.begin(); e != x.end(); ++e)
		cout << h1.toString((*e).first) << endl;

	cout << endl;
//	x.coarse(h1.createFromList({1}));
	x.refineTo(h1.createFromList({1,1,0,1}));
	for (auto e = x.begin(true); e != x.end(); ++e)
		cout << h1.toString((*e).first) << endl;

	*/
	auto t1 = high_resolution_clock::now();

	// Create random structure between level 5 and 10 coords.
	for (int i = 0; i < n_refinements; i++) {
		double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		double z = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		double _l = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		level_t l = double(_l * 7 + 3);
//		H3::pos_t pos3 = {x,y,z};
//		H2::pos_t pos2 = {x,y};
		coord_t c3 = h3.createFromPosition(l, {x, y, z});
		coord_t c2 = h2.createFromPosition(l, {x, y});
		f.refineTo(c3);
		//if (i > 1000) continue;
		u.refineTo(c2);
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();

	cout << "Created " << u.nElements() << " elements in " << duration << "ms.\n";
	t1 = high_resolution_clock::now();
	// Now iterate and set to a certain function

	for (auto e : f) {
		H3::pos_t pos = h3.getPosition(e.first);
		e.second = pos[1];
	}
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();

	cout << "Filled TLC in " << duration << "ms.\n";
	t1 = high_resolution_clock::now();
	f.propagate();
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Propagate() in " << duration << "ms.\n";

	for (auto e : u) {
		H2::pos_t pos2 = h2.getPosition(e.first);
		e.second = pos2[1];
	}
	u.propagate();

	// BC 0 = BC 1 = reflect (X+ X-
	// BC 3 = 0 (default) (Y+)
	// BC 2 = 1 (Y-)
	u.boundary[2] = [](ScalarField2Base *self, coord_t cc)->data_t { return 1;};
	u.boundary[0] = u.boundary[1] = [](ScalarField2Base *self, coord_t cc)->data_t { coord_t c = self->hcs.removeBoundary(cc); return self->get(c);};

	write_pgm("test31.pgm", u, 10);
	write_pgm_level("test32.pgm", u);

}
