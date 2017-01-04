#include "includes.hpp"

// TEST1: HCS neighbor+GetPosition speed

int main(int argc, char **argv) {

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
	cout << "Traversing " << count << " neighbors took " << duration << "ms. Or " << (count / (double)duration) << " lookups per ms.\n";

	t1 = high_resolution_clock::now();
	c =     h3.createFromList({0,0,0,0,0,0,0,0,0}); // Level 10 corner
	coord_t c_end = h3.createFromList({7,7,7,7,7,7,7,7,7}); // Level 10 corner

	count = 0;

	for (; c <= c_end; c++) {
		H3::pos_t pos = h3.getPosition(c);
		 assert(pos[0] < 1);
		 count++;
	}
	t2 = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(t2-t1).count();
	cout << "Level-9 pos lookup " << count << " times took " << duration << "ms. Or " << (count / (double)duration) << " lookups per ms.\n";

}
