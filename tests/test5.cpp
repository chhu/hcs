#include "includes.hpp"

#include <bitset>

int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif

	cout << "Interactive test for HCS and Field class\n";

	H1 h1;
	H2 h2;
	H3 h3;

	coord_t c;

	cout << h1.GetLevel(h1.CreateMaxLevel(8)) << endl;

	c = h1.createFromUnscaled(8, {127});
	cout << "1D-center level-8: " << h1.toString(c) << endl;

	c = h2.createFromUnscaled(8, {127,127});
	cout << "2D-center level-8: " <<  h2.toString(c) << endl;

	c = h3.createFromUnscaled(8, {127,127,127});
	cout << "3D-center level-8: " <<  h3.toString(c) << endl;

	cout << "3D Neighbor directions from HCS.getNeighborDirection():\n";
	for (uint8_t i = 0; i < 6; i++) {	// 6 neighbors for 3 dimensions 2 * D
		cout << "Neighbor direction " << (int)i << " points toward: " << Vec3(h3.getDirectionNormal(i)) << endl;
	}

	cout << "2D Coefficient Test\n";
	cout << "Operating on a complete 2D level 4 scalar field, X and Y scaled as unit-cube.\nPrints interpolation coeffs for a provided coord.\n";
	cout << "Value of the field is set to its X position.\n";

	SparseScalarField2 x;
	x.createEntireLevel(4);
	x.bracket_behavior = ScalarField2::BR_INTERP;
	for (auto e : x) {
		e.second = h2.getPosition(e.first)[0];//h2.getPosition(e.first)[0];	// Set the value of the field to the X Cartesian-coord
	}
	x.propagate();
	while (true) {
		level_t level;
		H2::pos_t uc;
		bool use_non_top;
		cout << "Level: "; cin >> level;
		cout << "X    : "; cin >> uc[0];
		cout << "Y    : "; cin >> uc[1];
		//cout << "Unscaled Z: "; cin >> uc[2];
	//	cout << "Use non-top: "; cin >> use_non_top;
		coord_t c = h2.createFromPosition(level, uc);
		cout << "Closest H2 coord: " << h2.toString(c) << endl;
		//ScalarField2::coeff_map_t coeffs;
		//x.coeff_down_count = x.coeff_up_count = 0;
		//x.getCoeffs(c, coeffs, use_non_top);
		//auto coeffs = h2.getCoeffCoords(c);
		auto coeffs = h2.getCoeffs(c);
		int count =  0;
		data_t sum = 0;
		for (auto coeff : coeffs) {
			cout << "Coeff " << count++ << ": " << h2.toString(coeff.first) << " Weight: " << coeff.second << endl;
			sum += coeff.second;
//			cout << "Coeff " << count++ << ": " << h2.toString(coeff) << endl;

		}
		count =  0;
		/*

		cout << "Non-Existent recursive calls: " << x.coeff_down_count <<
				" Top-averaging recursive calls: " << x.coeff_up_count <<
				" Total (should be 1): " << sum << endl <<
				" Value: " << x.get(c) << endl <<
				"Neighbor coords and values:\n";
		for (int ne_idx = 0; ne_idx < x.hcs.parts; ne_idx++) {
			coord_t ne_coord = h2.getNeighbor(c, ne_idx);
			cout << h2.toString(ne_coord) << " " << x.get(ne_coord) << endl;
		}
		*/

	}

}
