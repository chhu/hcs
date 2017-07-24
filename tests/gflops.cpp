#include "includes.hpp"

int main(int argc, char **argv) {

	if (argc != 2) {
		cout << "Use: ./gflops <MByte to test>\n";
		exit(1);
	}
	size_t mem_size = stoi(string(argv[1]), nullptr) * 1024 * 1024;
	size_t elements = mem_size / sizeof(double) / 2;

	vector<double> a(elements, 0.);
	vector<double> b(elements, 0.);

	for (size_t i = 0; i < elements; i++) {
		a[i] = double(elements) / (i + 1.);
		b[i] = 1. - a[i];
	}
	double result = 0;
	cout << "Dot product with 2 " << elements << " double arrays\n";
	auto t1 = high_resolution_clock::now();

	for (size_t i = 0; i < elements; i += 8)
		result += a[i + 0] * b[i + 0]	// Vectorizing hint
			    + a[i + 1] * b[i + 1]
				+ a[i + 2] * b[i + 2]
				+ a[i + 3] * b[i + 3]
				+ a[i + 4] * b[i + 4]
				+ a[i + 5] * b[i + 5]
				+ a[i + 6] * b[i + 6]
				+ a[i + 7] * b[i + 7];

	for (size_t i = 0; i < elements; i += 8)
		result += a[i + 0] * b[i + 0]	// Vectorizing hint
			    + a[i + 1] * b[i + 1]
				+ a[i + 2] * b[i + 2]
				+ a[i + 3] * b[i + 3]
				+ a[i + 4] * b[i + 4]
				+ a[i + 5] * b[i + 5]
				+ a[i + 6] * b[i + 6]
				+ a[i + 7] * b[i + 7];

	for (size_t i = 0; i < elements; i += 8)
		result += a[i + 0] * b[i + 0]	// Vectorizing hint
			    + a[i + 1] * b[i + 1]
				+ a[i + 2] * b[i + 2]
				+ a[i + 3] * b[i + 3]
				+ a[i + 4] * b[i + 4]
				+ a[i + 5] * b[i + 5]
				+ a[i + 6] * b[i + 6]
				+ a[i + 7] * b[i + 7];

	for (size_t i = 0; i < elements; i += 8)
		result += a[i + 0] * b[i + 0]	// Vectorizing hint
			    + a[i + 1] * b[i + 1]
				+ a[i + 2] * b[i + 2]
				+ a[i + 3] * b[i + 3]
				+ a[i + 4] * b[i + 4]
				+ a[i + 5] * b[i + 5]
				+ a[i + 6] * b[i + 6]
				+ a[i + 7] * b[i + 7];

	auto t2 = high_resolution_clock::now();

	auto duration_i = duration_cast<microseconds>(t2-t1).count();
	double duration = double(duration_i) / 1e6 / 4.;  // 4 times loop

	cout << "Took " << duration << "s. Or "
		 << (1e-9) * ((elements) / (double)duration) << " GFLOP/s (FMA is 1 FLOP).Memory bandwidth: "
		 << double(elements) * 2. * double(sizeof(double)) / double(duration) / 1024. / 1024. / 1024. << "GByte/s Result: " << result << endl;

}
