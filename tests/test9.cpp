#include "includes.hpp"
#include "solver.hpp"

// Dimension setup


void poisson(ScalarField1 &x, ScalarField1 &b) {
	b.propagate();
	level_t highest = x.getHighestLevel();
	H1 &h =x.hcs;
	data_t div = x.nElementsTop();

	ScalarField1 fake;	// easy way to obtain correct coeffs.

	// zero special since there is no "parent"
	data_t l_val = x.get(h.getNeighbor(0, 1), true);
	data_t r_val = x.get(h.getNeighbor(0, 0), true);
	x[0] = 0.5 * l_val + 0.5 * r_val + b[0]/4;
	for (level_t l = 1; l <= highest; l++) {

		for (auto it = x.begin(false, l); it != x.end(); ++it) {
			cout << x.hcs.toString((*it).first);
			coord_t coord = (*it).first;
			data_t &value = (*it).second;
			//value = 0;
			data_t b_delta = 0;// / (l + 4);// / (1 );
			ScalarField1::coeff_map_t coeffs;
			fake.getCoeffs(coord, coeffs, false);
			for (auto e : coeffs) {
				value += x.get(e.first, true) * e.second;
				b_delta += b.get(e.first, true) * e.second;
				//cout << "[ " << h.toString(e.first) << " :: " << e.second << " ]";
			}
			//b[coord] = b_delta;// + b.get(coord, true) / (l + 2);
			cout << " VAL: " << value <<  " B DELTA: " << b_delta << endl;
		}
		fake.clear();
		fake.createEntireLevel(l);
	}
}

int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif


	H1 h;

	// Resolution setup
	const int max_level = 3;	// Max level for "C" field



	ScalarField1 x, b;
	b.createEntireLevel(max_level);
	x.createEntireLevel(max_level);

	//x.boundary[0] = [](ScalarField1 *self, coord_t cc)->data_t { return 1;};

	x = 0;
	b = 0;
	b[h.createFromUnscaled(max_level, {4})] = 1;
	b.boundary[0] = b.boundary[1] = [](ScalarField1 *self, coord_t cc)->data_t { return 0;};
	b.propagate();

	poisson(x,b);
	poisson(x,b);

	write_txt("test9.txt", x);
	write_txt("test9b.txt", b,false);



}
