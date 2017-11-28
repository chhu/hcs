#include <stdio.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <chrono>

#ifdef USE_MPI
#include <mpi.h>
#endif


using namespace std;
using namespace chrono;

double dot(const vector<double> &a, const vector<double> &b) {
//  return std::inner_product(a.begin(), a.end(), b.begin(), 0);

	double result = 0;
	for (size_t i = 0; i < a.size(); i += 8)
		result += a[i + 0] * b[i + 0]	// Vectorizing hint
			    + a[i + 1] * b[i + 1]
				+ a[i + 2] * b[i + 2]
				+ a[i + 3] * b[i + 3]
				+ a[i + 4] * b[i + 4]
				+ a[i + 5] * b[i + 5]
				+ a[i + 6] * b[i + 6]
				+ a[i + 7] * b[i + 7];
	return result;
}

int main(int argc, char **argv) {

	size_t mem_size = 0;
	unsigned loops = 0;
	int rank = 0, cpu_count = 1;
#ifdef USE_MPI
	MPI_Init( &argc, &argv );
	MPI_Comm communicator = MPI_COMM_WORLD;

	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &cpu_count);
	if (rank == 0) printf("MPI Enabled, node %d of %d\nOutput is [min of single rank, max of single rank, sum]\nThis benchmark does not account for MPI overhead. A large dot product is performed that quickly evades CPU cache.\n", rank, cpu_count);
#endif

	if (argc != 3) {
		if (rank == 0) cout << "Use: ./gflops <MByte to test> <N loops>\nDefault: 100MB / rank, 20 loops.\n";
		mem_size = 100 * (size_t)1024 * 1024;
		loops = 20;
		
	} else {
		 mem_size = (size_t)stoi(string(argv[1]), nullptr) * (size_t)1024 * 1024;
		loops = stoi(string(argv[2]), nullptr);
	
	}
    

	size_t elements = mem_size / sizeof(double) / 2;
	vector<double> best_gflop;
	vector<double> a(elements, 0.);
	vector<double> b(elements, 0.);


	for (size_t i = 0; i < elements; i++) {
		a[i] = double(elements) / (i + 1.);
		b[i] = 1. - a[i];
	}
	if (rank == 0) cout << "Dot product with 2 " << elements * cpu_count << " double arrays\n";
	
#ifdef USE_MPI
	double single_ref_gflop = 1;
	if (rank == 0) {
	auto t1 = high_resolution_clock::now();
	volatile double result = 0;
	result += dot(a,b);
	result += dot(a,b);
	result += dot(a,b);
	result += dot(a,b);


	auto t2 = high_resolution_clock::now();

	auto duration_i = duration_cast<microseconds>(t2-t1).count();
	double duration = double(duration_i) / 1e6 / 4.;  // 4 times loop
	
	single_ref_gflop = (1e-9) * ((elements) / (double)duration);
	printf("Single reference for speedup: %5.3g\n", single_ref_gflop);
	}
	MPI_Barrier(communicator);
#endif

	while (loops--) {
	
	volatile double result = 0; // Avoid optimization
	
	auto t1 = high_resolution_clock::now();

	result += dot(a,b);
	result += dot(a,b);
	result += dot(a,b);
	result += dot(a,b);


	auto t2 = high_resolution_clock::now();

	auto duration_i = duration_cast<microseconds>(t2-t1).count();
	double duration = double(duration_i) / 1e6 / 4.;  // 4 times loop
	
	double gflops = (1e-9) * ((elements) / (double)duration);
	
	double mem_speed = double(elements) * 2. * double(sizeof(double)) / double(duration) / 1024. / 1024. / 1024.; // GB/s
#ifdef USE_MPI
	double send_data[2] = {gflops, mem_speed};
	double receive_data_sum[2];
	double receive_data_min[2];
	double receive_data_max[2];
	MPI_Allreduce(send_data, receive_data_sum, 2, MPI_DOUBLE, MPI_SUM, communicator);
	MPI_Allreduce(send_data, receive_data_min, 2, MPI_DOUBLE, MPI_MIN, communicator);
	MPI_Allreduce(send_data, receive_data_max, 2, MPI_DOUBLE, MPI_MAX, communicator);
	best_gflop.push_back(receive_data_sum[0]);
	if (rank == 0)	printf("Run %u: [%5.3g, %5.3g, %5.3g] GFLOP/s (FMA is 1 FLOP), [%5.3g, %5.3g, %5.3g] GByte/s. Speedup: %.4g\n", loops, receive_data_min[0],  receive_data_max[0], receive_data_sum[0],  receive_data_min[1], receive_data_max[1],  receive_data_sum[1], receive_data_sum[0] / single_ref_gflop); 
#else
	best_gflop.push_back(gflops);
	if (rank == 0)	printf("Run %u: %5.3g GFLOP/s (FMA is 1 FLOP), %5.3g GByte/s.\n", loops, gflops, mem_speed); 
#endif
	}
	if (rank == 0) printf("Peak performance: %g GFLOP/s.\n", *max_element(best_gflop.begin(), best_gflop.end()));
#ifdef USE_MPI
	MPI_Finalize();
#endif
}

