#include "includes.hpp"


int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif

	H1 h1;
	H2 h2;
	H3 h3;

	const int max_level = 4;

	ScalarField1 x;
	x.createEntireLevel(max_level);


	x = 0;

	coord_t cpos = h1.createFromPosition(3,{0.42});
	x.coarse(cpos);
/*
	map<coord_t, data_t> result;
	x.getCoeffs(cpos+1, result, false);
	for (auto it : result)
		cout << h1.toString(it.first) << " " << it.second << endl;
//exit(2);
*/
	ScalarField1 b = x;
	b = 1.;

	for (auto it = x.begin(true); it != x.end(); ++it)
		cout << h1.toString((*it).first) << " VAL: " << (*it).second << endl;
	/*
	x.boundary[0] = [](ScalarField1 *self, coord_t c)->data_t {
		//return self->get(self->hcs.removeBoundary(c)); // Neumann, derivative == 0
		return 1;	// Dirichlet to 1
	};
	*/


	auto stencil = [&x](coord_t coord)->map<coord_t, data_t>  {
		map<coord_t, data_t> result;
		coord_t level = (coord_t)1 << x.hcs.GetLevel(coord);
		data_t dist = 1. / level;  // Neighbor distance, assuming all directions with equal scale and scale = 1
		data_t vol = dist;
		//data_t row_result = 0;
		for (int ne_idx = 0; ne_idx < x.hcs.parts; ne_idx++) {
			coord_t ne_coord = x.hcs.getNeighbor(coord, ne_idx);
			map<coord_t, data_t> current;
			x.getCoeffs(ne_coord, current, false);

			data_t coeff  = (1. ) / (dist * vol);

			if (x.hcs.IsBoundary(ne_coord)) {
				coeff *= 2;  // If we hit a boundary, its only half the distance
			}

			result[coord] += coeff;
			for (const auto &e : current) {
				result[e.first] -= coeff * e.second;
			}
		}
		return result;
	};

	ScalarField1::coeff_map_t cf;
	map<coord_t, int> coord_to_row;
	int row = 0;
	int bc_row = x.nElementsTop();
	for (auto e = x.begin(true); e != x.end(); ++e) {
		coord_t c = (*e).first;
		coord_to_row[c] = row++;
		for (int ne_idx = 0; ne_idx < x.hcs.parts; ne_idx++) {
			coord_t ne_coord = x.hcs.getNeighbor(c, ne_idx);
			if (x.hcs.IsBoundary(ne_coord)) {
				coord_to_row[ne_coord] = bc_row++;
			}
		}
	}
	int n_all = bc_row;
	double mat[18][18] {0};

	row = 0;
	for (auto e = x.begin(true); e != x.end(); ++e) {
		cout << h1.toString((*e).first) << " >>> ";
		cf = stencil((*e).first);
		for (auto it = cf.begin(); it != cf.end(); ++it) {
			coord_t coeff_coord = it->first;
			data_t coeff = it->second;
			auto cit = coord_to_row.find(coeff_coord);
			if (cit == coord_to_row.end())
				cout << " U!" << h1.toString(coeff_coord) << " " << coeff << " |";
			else {
				int col = coord_to_row[coeff_coord];
				printf("(%d, %d) %g ", col, row, coeff);
				mat[col][row] = coeff;
			}
		}
		printf("\n");
		row++;
	}
	while (bc_row-- != row)
		mat[bc_row][bc_row] = 1.;
//	printf("(%d, %d) %g\n ", bc_row, bc_row, 1.);

	printf("\n");
	for (int j = 0; j < 18; j++) {
		for (int i = 0; i < 18; i++)
			printf("%g ", mat[i][j]);
		printf("\n");
	}
	printf("\n");



}
