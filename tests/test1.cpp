#include "includes.hpp"

// TEST1: HCS neighbor+GetPosition speed

int main(int argc, char **argv) {

#ifdef __BMI2__
	printf("Using BMI2 instructions. (AMD screwed up BMI2 in Zen and Zen+ cores! Deactivate manually.) \n");
#endif

	H2 h2;
	H3 h3;
	coord_t c = h3.createFromList({0,0,0,0,0,0,0,0,0,0}); // Level 10 corner
	coord_t origin = c;
	coord_t current;
	uint64_t count = 0;

	auto t1 = high_resolution_clock::now();
	for (int i = 0; i < 100000; i++) {
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor(c, 0); // X+
			count++;
		}
		c = current;  // one back
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor(c, 2); // Y+
			count++;
		}
		c = current;
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor(c, 4); // Z+
			count++;
		}
		c = current;
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor(c, 3); // Y-
			count++;
		}
		c = current;  // one back
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor(c, 5); // Z-
			count++;
		}
		c = current;
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor(c, 1); // X-
			count++;
		}
		c = current;
		assert(c == origin);	// We should be back where we started
	}
	auto t2 = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Traversing " << count << " neighbors with getNeighbor() took " << duration << "ms. Or " << (count / (double)duration) << " lookups per ms.\n";

	count = 0;
	t1 = high_resolution_clock::now();
	for (int i = 0; i < 100000; i++) {
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor2(c, 0); // X+
			count++;
		}
		c = current;  // one back
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor2(c, 2); // Y+
			count++;
		}
		c = current;
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor2(c, 4); // Z+
			count++;
		}
		c = current;
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor2(c, 3); // Y-
			count++;
		}
		c = current;  // one back
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor2(c, 5); // Z-
			count++;
		}
		c = current;
		while (!h3.IsBoundary(c)) {
			current = c;
			c = h3.getNeighbor2(c, 1); // X-
			count++;
		}
		c = current;
		assert(c == origin);	// We should be back where we started
	}
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Traversing " << count << " neighbors with getNeighbor2() took " << duration << "ms. Or " << (count / (double)duration) << " lookups per ms.\n";


	c = h3.CreateMinLevel(9);
	coord_t c_end = h3.CreateMaxLevel(9);


	cout << h3.toString(c) << endl;
	cout << h3.toString(c_end) << endl;

	t1 = high_resolution_clock::now();

	count = 0;
	data_t sum = 0;

	for (; c <= c_end; c++) {
		H3::pos_t pos = h3.getPosition(c);
		sum += pos[0] + pos[1] + pos[2]; // prevnts opt-out of loop
		count++;
	}
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Level-9 pos lookup " << count << " times took " << duration << "ms. Or " << (count / (double)duration) << " lookups per ms. #" << int(sum) << endl;

	t1 = high_resolution_clock::now();
	c = h3.CreateMinLevel(8);
	c_end = h3.CreateMaxLevel(8);
	count = 0;
	sum = 0;

	for (; c <= c_end; c++) {
		auto coeffs = h3.getCoeffs(c);
		for (auto c : coeffs)
			sum += c.second;
		count++;
	}
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout.precision(20);
	cout << "Level-8 coeff lookup " << count << " times took " << duration << "ms. Or " << int(count / (double)duration) << " lookups per ms. #" << (sum) << endl;

}
