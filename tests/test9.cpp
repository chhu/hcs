#include "includes.hpp"
#include "solver.hpp"

// Dimension setup

void poisson_err(ScalarField1 &err, ScalarField1 &x, ScalarField1 &b) {
	err = x;
	for (auto e : x) {
		coord_t c = e.first;
		data_t &value = e.second;
		data_t dist = x.hcs.getDistance(c, 0);
		data_t sum = 0;
		for (int ne = 0; ne < x.hcs.parts; ne++) {
			coord_t ne_coord = x.hcs.getNeighbor(c, ne);
			data_t ne_value = x.get(ne_coord, true);
			data_t ne_dist = x.hcs.IsBoundary(ne_coord) ? 0.5 * dist : dist;
			sum += (ne_value - value) / (ne_dist * ne_dist);
		}
		err[c] = b[c] - sum;
	}
}

void poisson(ScalarField1 &x, ScalarField1 &b) {
	level_t highest = x.getHighestLevel();
//b/=pow(2,highest);
	//b.subHarmonics();
	//b.propagate();
	x.propagate();
	H1 &h =x.hcs;
	data_t div = x.nElementsTop();
	//b = b*b;
	ScalarField1 fake;	// easy way to obtain correct coeffs.
//x=b;
	// zero special since there is no "parent"
	data_t l_val = x.get(h.getNeighbor(0, 1), true);
	data_t r_val = x.get(h.getNeighbor(0, 0), true);
	//x[0] = 0.5 * l_val + 0.5 * r_val + b[0];
	x[0] = b[0] * 0.25;
	for (level_t l = 1; l <= highest; l++) {
		data_t level_dist = 0.5/ (1U << l);
		for (auto it = x.begin(false, l); it != x.end(); ++it) {
			//cout << x.hcs.toString((*it).first);
			coord_t coord = (*it).first;
			data_t &value = (*it).second;
			value = 0;
			//value *= level_dist;
			//data_t b_delta = b[coord];// - bb[x.hcs.ReduceLevel(coord)] ;// / (l + 4);// / (1 );
			ScalarField1::coeff_map_t coeffs;
			fake.getCoeffs(coord, coeffs, false);
			for (auto e : coeffs) {
				value += x.get(e.first, true) * e.second ;// + b[e.first];//.get(e.first, true) * e.second;
				//b_delta += b.get(e.first, true) * e.second;
				//cout << "[ " << h.toString(e.first) << " :: " << e.second << " ]";
			}

			//value += b_delta;
//			b[coord] = value;// + b.get(coord, true) / (l + 2);
			//cout << " VAL: " << value <<  " B ORIG: " << b[coord] << " B DELTA: " << b_delta /l<< endl;
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
	const int max_level = 9;	// Max level for "C" field



	ScalarField1 x, b, err;
	b.createEntireLevel(max_level);
	err.createEntireLevel(max_level);
	x.createEntireLevel(max_level);

	//x.boundary[0] = [](ScalarField1 *self, coord_t cc)->data_t { return 1;};

	x = 0;
	b = 0;
	//b[0] = 1/8.;
	//b[h.createFromList({0})] = 1 / 8.;
	//b[h.createFromList({1})] = 1 / 8.;
	//b[h.createFromList({0,1,1,1,1,1,1})] = 1 / 8.;
	//b[h.createFromList({1,0,0,0,0,0,0})] = 1 / 8.;
//	b[h.createFromList({1,0})] = 1 / 8.;
//	b[h.createFromList({1,1})] = 1 / 8.;
	//b[h.createFromList({0,0})] = 1 / 8.;
//	b[h.createFromList({0,1})] = 1 / 8.;
//	b[h.createFromList({0,1,1})] = 1 / (32.*32.);
	/*
	b[h.createFromList({0,1,1,1})] = 1 / 64.;
	b[h.createFromList({0,1,1,1,1})] = 1 / 128.;
	b[h.createFromList({0,1,1,1,1,1})] = 1 / 256.;
	b[h.createFromList({0,1,1,1,1,1,1})] = 1 / 512.;
	b[h.createFromList({0,1,1,1,1,1,1,1})] = 1 / 1024.;
	b[h.createFromList({0,1,1,1,1,1,1,1,1})] = 1 / 2048.;
/*
	b = 0;
	b[0] = 1/2048.;//1/4.;
	b[h.createFromList({0})] = 1/2048.;//1 / 8.;
	b[h.createFromList({0,1})] = 1/2048.;//1 / 16.;
	b[h.createFromList({0,1,1})] = 1/2048.;//1 / 32.;
	b[h.createFromList({0,1,1,1})] = 1/2048.;//1 / 64.;
	b[h.createFromList({0,1,1,1,1})] = 1/2048.;//1 / 128.;
	b[h.createFromList({0,1,1,1,1,1})] = 1/2048.;//1 / 256.;
	b[h.createFromList({0,1,1,1,1,1,1})] = 1/2048.;//1 / 512.;
	b[h.createFromList({0,1,1,1,1,1,1,1})] = 1/2048.;//1 / 1024.;
	b[h.createFromList({0,1,1,1,1,1,1,1,1})] = 1 / 2048.;
	//b[h.createFromUnscaled(max_level, {4})] = 1;*/
	b.boundary[0] = b.boundary[1] = [](ScalarField1Base *self, coord_t cc)->data_t {
	//	return self->get(self->hcs.removeBoundary(cc));
		return 0;
	};

	// B = f(x) = x**3
	// Analytic solution for Dirichlet = 0: x / 20 - x**5/20
	//for (auto e : b)
	//	e.second = pow(h.getPosition(e.first)[0], 3);
	//b.propagate();
	//b[h.createFromPosition(max_level, {0.5})] = 1;
	b = 1;
	for (int i = 0; i < 1; i++)
		poisson(x,b);
	x.propagate();
	b.propagate();

	poisson_err(err,x,b);
	//b = 1;

	//b.subHarmonics();

	write_txt("test9e.txt", err,false);
	write_txt("test9.txt", x);
	write_txt("test9b.txt", b,false);



}
