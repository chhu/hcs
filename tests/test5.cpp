#include "includes.hpp";



int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#endif


	// Test suite
	H1 h1;
	H2 h2;
	H3 h3;

	coord_t c = h2.createFromUnscaled(9, {127,256});
	cout << h2.toString(c) << endl;

	c = h3.createFromUnscaled(8, {127,127,127});
	cout << h3.toString(c) << endl;

	c = h1.createFromUnscaled(8, {127});
	cout << h1.toString(c) << endl;

	ScalarField2 x('x', &h2);
	x.createEntireLevel(8);
	cout << "Coefficient Test\n";
	cout << "Operating on a complete 2D level 8 scalar field, X and Y scaled as unit-cube.\nPrints interpolation coeffs for a provided coord.\n";
	while (true) {
		level_t level;
		H2::pos_t uc;
		bool use_non_top;
		cout << "Level: "; cin >> level;
		cout << "X    : "; cin >> uc[0];
		cout << "Y    : "; cin >> uc[1];
		//cout << "Unscaled Z: "; cin >> uc[2];
		cout << "Use non-top: "; cin >> use_non_top;
		coord_t c = h2.createFromPosition(level, uc);
		cout << "Closest H2 coord: " << h2.toString(c) << endl;
		ScalarField2::coeff_map_t coeffs;
		x.coeff_down_count = x.coeff_up_count = 0;
		x.getCoeffs(c, coeffs, use_non_top);
		int count =  0;
		data_t sum = 0;
		for (auto coeff : coeffs) {
			cout << "Coeff " << count++ << ": " << h2.toString(coeff.first) << " Weight: " << coeff.second << endl;
			sum += coeff.second;
		}
		cout << "Non-Existent recursive calls: " << x.coeff_down_count <<
				" Top-averaging recursive calls: " << x.coeff_up_count <<
				" Total (should be 1): " << sum << endl;

	}

}
