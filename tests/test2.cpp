#include "includes.hpp"

int main(int argc, char **argv) {


	// Test suite
	H3 h3;

	VectorField3 v1('v', &h3);
	VectorField3 v2('v', &h3);
	ScalarField3 x('u', &h3);
	cout << sizeof(Vec4) << endl;

	v1.createEntireLevel(8);
	v2.createEntireLevel(8);

	cout << sizeof(Vec3) << endl;
	auto t1 = high_resolution_clock::now();
	uint64_t count = 0;
	for (auto e : v1) {
		Vec3 vec({count + 1, count , -(data_t)count});
		vec.normalize();
		e.second = vec;
		count++;
	}
	for (auto e : v2) {
		Vec3 vec({count + 1, count , -(data_t)count});
		vec.normalize();
		e.second = vec;
		count++;
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();

	cout << "Setting vectors of 2 complete level-8 fields took " << duration << "ms.\n";

	t1 = high_resolution_clock::now();
	x.merge<Vec3>(v1, v2, [](coord_t c, Vec3 v1v, Vec3 v2v)->data_t {
		return v1v * v2v;
	});

	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Merging dot product into new ScalarField of level 8 took " << duration << "ms.\n";

	t1 = high_resolution_clock::now();
	x.convert<Vec3>(v1,[](coord_t c, Vec3 v1v)->data_t {
		return v1v.length();
	});

	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Converting vector field into new ScalarField of level 8 took " << duration << "ms.\n";

	t1 = high_resolution_clock::now();
	ScalarField3 y = x;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Copy of ScalarField of level 8 took " << duration << "ms.\n";

	t1 = high_resolution_clock::now();
	y *= x;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Multiply *= ScalarField of level 8 took " << duration << "ms.\n";

	ScalarField3 l7s('x', &h3);
	l7s.createEntireLevel(7);
	t1 = high_resolution_clock::now();
	y *= l7s;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Multiply *= ScalarField of level 8 with level 7 took " << duration << "ms.\n";

	ScalarField3 l6s('x', &h3);
	l6s.createEntireLevel(6);
	t1 = high_resolution_clock::now();
	y *= l6s;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Multiply *= ScalarField of level 8 with level 6 took " << duration << "ms.\n";
}
