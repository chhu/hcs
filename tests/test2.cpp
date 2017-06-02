#include "includes.hpp"

int main(int argc, char **argv) {

	// Test suite
	H3 h3;

	VectorField3 v1;
	VectorField3 v2;
	ScalarField3 x;
	cout << "3D - level 8 test, fully populated 256x256x256 box\n" << endl;

	v1.createEntireLevel(8);
	v2.createEntireLevel(8);
	x.takeStructure(v1);

	int v3size = sizeof(Vec3);
	cout << "Single vector3 size: " << v3size << "bytes, nTop = " << v1.nElementsTop() << endl;
	auto t1 = high_resolution_clock::now();
	Vec3 vec({1,2,3});
	v1 = vec;
	v2 = vec;
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();

	double total_bytes_written = v3size * v1.nElements();
	duration /= 2;
	cout << "Setting vector field to a constant took " << duration << "ms.\n";
	cout << "Throughput: " << (double)(total_bytes_written / duration * 1000) / 1024 / 1024 << " MByte/s\n";

	t1 = high_resolution_clock::now();
	v1.propagate();
	v2.propagate();
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Propagating vector field of level 8 took " << duration /2 << "ms.\n";
	t1 = high_resolution_clock::now();
	x.merge<Vec3>(v1, v2, [](coord_t c, Vec3 v1v, Vec3 v2v)->data_t {
		return v1v * v2v;
	});

	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Merging dot product into new ScalarField of level 8 took " << duration << "ms.\n";

	t1 = high_resolution_clock::now();
	x.convert<Vec3>(v1,[](coord_t c, VectorField3 &s)->data_t {
		return s.get(c).length();
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
	x = y;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Copy of ScalarField of level 8 into equal level took " << duration << "ms.\n";

	t1 = high_resolution_clock::now();
	y *= x;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Multiply *= ScalarField of level 8 took " << duration << "ms.\n";

	ScalarField3 l7s;
	l7s.createEntireLevel(7);
	t1 = high_resolution_clock::now();
	y *= l7s;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Multiply *= ScalarField of level 8 with level 7 took " << duration << "ms.\n";

	ScalarField3 l6s;
	l6s.createEntireLevel(4);
	t1 = high_resolution_clock::now();
	y *= l6s;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Multiply *= ScalarField of level 8 with level 6 took " << duration << "ms.\n";
}
