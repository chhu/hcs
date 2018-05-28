/*
 * test_helper.hpp
 *
 *  Created on: Dec 28, 2016
 *      Author: snk
 */

#pragma once

// STL includes

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <set>
#include <cassert>
#include <list>
#include <algorithm>
#include <valarray>
#include <string>
#include <map>
#include <functional>
#include <numeric>
#include <cstring>
#include <limits>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <array>
#include <memory>
#include <chrono>
#include <functional>
#include <bitset>


// Own includes
#include "hcs.hpp"
#include "tensor.hpp"
#include "field.hpp"
#include "sparsefield.hpp"
#include "densefield.hpp"
#include "numerics.hpp"

using namespace std;
using namespace chrono;


// Write a 2D scalar field to a PGM. Scales automatically.
// Level determines resolution.
void write_pgm(string filename, ScalarField2 &field, level_t level) {
	auto& hcs = field.hcs;

	int width = 1U << level;
	int height= 1U << level;

	valarray<data_t> buffer(width * height);

	coord_t level_start = hcs.CreateMinLevel(level);
	field.bracket_behavior = ScalarField2::BR_INTERP;
	for (coord_t c = level_start; c < level_start + width * height; c++) {
		H2::unscaled_t pos = hcs.getUnscaled(c);
		int x = pos[0];
		int y = pos[1];
		buffer[x + width * y] = field.get(c, true);
	}
	data_t f_min = buffer.min();
	data_t f_max = buffer.max();
	cout << "Field min: " <<f_min << " max: " << f_max << " sum: " << buffer.sum() << endl;

	// Byte-scale
	if (f_min == f_max)
		f_min = f_max - 1;
	else {
		buffer = (buffer - f_min) / (f_max - f_min);
		buffer *= 255;
	}

	ofstream f(filename);
	f << "P2" << endl << "# MIN: " << f_min << " MAX: " << f_max << endl << width << " " << height << endl << "255" << endl;
	for (int v = 0; v < width * height; v++) {
		f << (int)buffer[v] << " ";
		if (v % 20 == 0)
			f << endl;
	}
	f.close();
}

// Writes a PGM with the resolution of the highest level found in field, marking the level of each TL coord.
void write_pgm_level(string filename, ScalarField2 &field) {
	level_t highest = field.getHighestLevel();
	coord_t c_lo = field.hcs.CreateMinLevel(highest);
	coord_t c_hi = field.hcs.CreateMaxLevel(highest);
	DenseScalarField2 level_field;
	level_field.createEntireLevel(highest);
	for (coord_t c = c_lo; c < c_hi; c++) {
		coord_t cc = c;
		while (!field.exists(cc))
			cc = field.hcs.ReduceLevel(cc);
		level_field[c] = field.hcs.GetLevel(cc);
	}
	write_pgm(filename, level_field, highest);
}

void write_txt(string filename, ScalarField1 &field, bool top_only = true) {
	ofstream f(filename);
	auto& hcs = field.hcs;
	map<double, double> x_val;	// to output sorted...
	for (auto it = field.begin(top_only); it != field.end(); ++it) {
		double x = hcs.getPosition((*it).first)[0];
		x_val[x] = (*it).second;
	}
	for (auto e : x_val)
		f << e.first << " " << e.second << endl;
	f.close();
}



void read_raw2(string filename, ScalarField2 &field, level_t level) {
    ifstream f(filename);
    unsigned n = 1U << level;
    cout << "Reading " << filename << " with " << n << "² elements.\n";
    field.clear();
    field.createEntireLevel(level);
    for (int y = 0; y < n; y++)
        for (int x = 0; x < n; x++) {
            coord_t c = field.hcs.createFromUnscaled(level, {x,y});
            f.read((char *)&field[c], sizeof(data_t));
        }
    f.close();
}

