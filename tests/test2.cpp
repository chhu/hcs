#include "includes.hpp"
#include <bitset>
int main(int argc, char **argv) {

//	H3 h;
	ScalarField3 x, o,m;

	x.createEntireLevel(8);
	x = 111;
    o.createEntireLevel(8);
    m.createEntireLevel(8);
	o = 89;
	m = x + o;
	cout << m[1] << endl;
//	return 0;
	size_t count = 0;
	auto t1 = high_resolution_clock::now();
	for (int i = 0; i < 1; i++) {
		count += x.nElementsTop() / 10;
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "10 counts took " << duration << "ms. "<< count << "\n";

/*
	coord_t tot = 0;
	coord_t start = 0;
	for (int i = 0; i <= 9; i++) {
		coord_t min_l = h.CreateMinLevel(i);
		coord_t max_l = h.CreateMaxLevel(i);
		bitset<32> min_(min_l);
		bitset<32> max_(max_l);
		coord_t diff = max_l - min_l + 1;
		cout << "LEVEL " << i << endl;
		cout << min_ << endl << max_ << endl << min_l << " " << max_l << "  " << diff <<  endl;
		bitset<32> start_(start);
		cout << start_ << "   " << start << " -> " << start + diff << endl;
		cout << "C2I " << h.coord2index(min_l);
		cout << "---" << endl;
		tot += diff;
		start += diff;

		coord_t tot = 0;
	}
	cout << tot << endl;
//	return 0;
*/
	// Test suite
	x.clear();
	H3 h3;

	VectorField3 v1;
	VectorField3 v2;
	//ScalarField3 x;
	cout << "3D - level 8 test, fully populated 256x256x256 box\n" << endl;

	v1.createEntireLevel(8);
	v2.createEntireLevel(8);
	x.takeStructure(v1);

	int v3size = sizeof(Vec3);
	cout << "Single vector3 size: " << v3size << "bytes, nTop = " << v1.nElementsTop() << endl;
	t1 = high_resolution_clock::now();
	Vec3 vec({1,2,3});
	v1 = vec;
	v2 = vec;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();

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
	x.convert<Vec3>(v1,[](coord_t c, VectorField3Base &s)->data_t {
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

/*
	ScalarField3 l6s;
	l6s.createEntireLevel(6);
	t1 = high_resolution_clock::now();
	y *= l6s;
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Multiply *= ScalarField of level 8 with level 6 took " << duration << "ms.\n";
*/
	// Setup raw-upscale test
	HCS<6> h;
	const int source_level = 3;

	// Source Field
	coord_t ls_min = h.CreateMinLevel(source_level);
	coord_t ls_max = h.CreateMaxLevel(source_level);

	vector<data_t> lsraw(ls_max - ls_min + 1, 1.);

	// Target Field
	coord_t lt_min = h.CreateMinLevel(source_level + 1);
	coord_t lt_max = h.CreateMaxLevel(source_level + 1);

	vector<data_t> ltraw(lt_max - lt_min + 1, 0.);

	t1 = high_resolution_clock::now();
	for (coord_t c = lt_min; c <= lt_max; c++) {
		auto coeffs = h.getCoeffs(c);
		size_t target_idx = c - lt_min;
		for (auto coeff : coeffs) { // coeffs _must_ only contain L8 coords because we queried a L9 coord
			if (h.IsBoundary(coeff.first)) { // assume all boundaries have 1
				ltraw[target_idx] += coeff.second;
			} else {
				size_t source_idx = coeff.first - ls_min;
				ltraw[target_idx] += lsraw[source_idx] * coeff.second;
			}
		}
	}
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Raw upscale of " << h.GetDimensions() << "D L" << h.GetLevel(ls_min) << "->L" << h.GetLevel(lt_min) << " took " << duration << "ms. Target Elements : " << ltraw.size() << " Interpolations/ms : " << data_t(ltraw.size()) / (duration) << endl;

	bool err = false;
	for (coord_t c = lt_min; c <= lt_max; c++) {
		size_t target_idx = c - lt_min;
		if (abs(1. - ltraw[target_idx]) > 1e-15) {
			cout << "Interpolation Error: " << h.toString(c) << " ERR: " << abs(1.-ltraw[target_idx]) << endl;
			break;
		}
	}
}
