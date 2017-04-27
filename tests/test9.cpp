#include "includes.hpp"
#include "solver.hpp"


template <typename DTYPE, typename HCSTYPE>
class PField: public Field<DTYPE, HCSTYPE> {
	void poisson(Field<DTYPE, HCSTYPE> &b) {

	}
};


int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif



	// Dimension setup
	typedef H1 HCS;
	typedef PField<data_t, H1> ScalarField;

	// Resolution setup
	const int max_level = 9;	// Max level for "C" field

	HCS hcs;

	ScalarField c;
	c.createEntireLevel(4);
	for (auto e : c) e.second = 1;
	write_txt("test9.txt", c);



}
