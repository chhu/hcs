

/* The H-Coordinate System
 * =======================
 *
 * This include / template library provides an extended, multi-dimensional recursive coordinate system, based on
 * the Z order space-filling curve or Morton-codes. It extends the Z curve as it captures the recursive nature
 * of the H fractal. Each recursive refinement is called level.
 * The coordinate type can be configured here and is intentionally not a template parameter. The bits available in
 * the coordinate type in combination with the dimensions you want to use determines how many levels you can go.
 *
 * The default coordinate type is unsigned 64bit, allowing a 3D recursion depth of 19 iterations (levels),
 * resulting in a closest distance of 1/2^19 for a scale of one (1x1x1 box). There is a level-marker bit that
 * determines the level of the coordinate, the lesser-significant bits are Morton codes for that level. This allows
 * almost linear storage between levels, so coord can be used as an array index.
 * The most significant bit marks a boundary coord, with the following N bits marking which boundary was hit, while
 * the remaining bits still describe a valid interior coordinate that "hit" or requested the boundary.
 *
 * Because the H fractal best illustrates this, its called the H-coordinate-system
 *
 *      6|        |7
 *       |        |
 *       |--------|
 *       |   /    |
 *   2| 4|  #  |3 |5
 *    |    /   |
 *    |--------|
 *    |        |
 *   0|        |1
 *
 *    Z+
 *    |   Y+
 *    |  /
 *    | /
 *    |------> X+
 *
 *    A single coordinate for a level (3-bit, 0-7 values for 3D) is as follows: LSB is X+/-, next bit Y+/-, next bit Z+/-
 *	  Example:
 *	    Level-1 coordinate:                                                  X        Y      Z
 *	                    L1
 *	     0b0000 .. 001 001            Would result in position (center) + (+scale, -scale, -scale) / 2
 *	                   ZYX
 *	                 ^ LEVEL-MARKER-BIT
 *
 *	    Level-2 coordinate:
 *	                    L1  L2
 *	     0b0000 .. 001 011 110        Would result in position (center) + (+scale, -scale, +scale) / 2 + (-scale, +scale, +scale) / 4
 *	                   ZYX ZYX
 *	                 ^ LEVEL-MARKER-BIT
 *
 *	 Remark: The LSBs (Least Significant Bit) are always the ones for the highest level (from level part).
 *	 This guarantees shortest linear distance in memory for neighbors and is compatible with Morton-codes.
 *
 *	 What does "unscaled" mean:
 *	 Methods like getPosition() or createFromPosition() use the scalings provided in center + scales to
 *	 compute the H coordinate for a certain level in floating-point units.
 *	 Unscaled means they operate on integers representing the whole level. A level-8 coordinate in 2D has
 *	 2^8 x 2^8 possible locations, so an unscaled level-8 coordinate would be in _unscaled_ Cartesian space from
 *	 X= 0 -> 255 and Y=0 -> 255, while a level-9 has 2^9, so X= 0 -> 511 and Y= 0 -> 511. An unscaled coordinate
 *	 has a completely different "true" location than the same coordinate at a different level!
 *
 *	 The HCS template library is independent of most other libs, exception is the toString() method.
 *	 Intel/AMD CPUs greatly profit from BMI2 instruction set.
 */

#pragma once

#ifdef __BMI2__
#include <immintrin.h>

#endif

namespace hcs {

// Configured to store coordinate in 64 bit
// The uint type used limits the maximal level to store, also depending on how many dimensions you use.
typedef uint64_t coord_t;

#define HCS_COORD_BITS sizeof(coord_t) * 8

// Configure in this routine the fastest way to count leading zeros of your coord_t.
// This depends on compiler and CPU. The reference implementation is terribly slow.
// Check out https://en.wikipedia.org/wiki/Find_first_set
inline  __attribute__((always_inline)) int __count_leading_zeros(const coord_t c) {
	//return c == 0 ? 64 :__builtin_clzll(c);
	return __builtin_clzll(c);
	//return __lzcnt64(c);
};


/*
 * Boundary coordinates
 * ====================
 * Generally, the "Special" bit (highest bit, or HCS_COORD_BITS-1) marks a boundary coordinate.
 * <dimension> bits are reserved for boundary direction from bits HCS_LEVEL_BIT-2 .. HCS_LEVEL_BIT-2 - D,
 * telling you what boundary is meant. The embedded coordinate is the one inside the domain that "hit"
 * the boundary through neighbor search.
 */

// A type used to carry separated level information, should be able to count to HCS_COORD_BITS,
// so all unsigned integer types will do.
typedef uint16_t level_t;

// Data precision
typedef double data_t;

using namespace std;

// The H-coordinate system (HCS) class just stores the scaling and position parameters and calculates some other useful stuff.
// The HCS class does not store data!

template<size_t dimensions = 3>
class HCS {
public:

	HCS() {
		part_mask = (1U << dimensions) - 1;
		parts = 1U << dimensions;
		// Initializes origin + scale
		for (int d = 0; d < dimensions; d++)
			center[d] = scales[d] = 0.5;	// makes 1x1x... box from 0-> 1
		max_level = (HCS_COORD_BITS - 2 - dimensions) / dimensions;
		boundary_mask = ~(coord_t)0 << (HCS_COORD_BITS - dimensions - 1);
		for (int d = 0; d < dimensions * 2; d += 2) {
			coord_t single = ~(1U << (d / 2)) & part_mask;
			successor_mask[d] = single;
			for (level_t l = 1; l < max_level; l++)
				successor_mask[d] |= single << (dimensions * l);
			successor_mask[d + 1] = ~successor_mask[d];
			bmi_mask[d >> 1] = successor_mask[d+1] & ~boundary_mask;

		}
		level_bounds[0].first = level_bounds[0].second = 1;

		for (int l = 1; l < level_bounds.size(); l++) {
			level_bounds[l].first = CreateMinLevel(l);
			level_bounds[l].second = CreateMaxLevel(l);
		}
	}

	// A type to store a Cartesian position
	typedef array<data_t, dimensions> pos_t;
	typedef array<uint32_t, dimensions> unscaled_t; // and the raw morton->Cartesian coords

	pos_t center;
	pos_t scales;

	coord_t part_mask;  // A bit mask that covers a single level
	level_t parts;		// How many directions has a single level (nodes on the H)

	coord_t boundary_mask; // masks all bits representing boundary information
	level_t max_level;	// The highest recursion depth (level) of this setup

	array<coord_t, dimensions * 2> successor_mask;
	array<coord_t, dimensions> bmi_mask;
	array<pair<coord_t, coord_t>, 64> level_bounds; // min-max pair of coords for each level

	// Test for most-significant bit
	static bool IsBoundary(coord_t coord) {
		return __count_leading_zeros(coord) == 0;
		//return bool(coord & ((coord_t)1 << (HCS_COORD_BITS - 1)));
	}

	// returns true if coord is int the dead-zone.
	static bool IsValid(coord_t coord) {
		if (coord == 0)
			return false;
		uint32_t first_bit = HCS_COORD_BITS - __count_leading_zeros(coord) - 1;
	// 2D:	1 3 5 7 9
	// 3D:  1 4 7 10 13 -1: 0 3 6 9 12
		return first_bit % dimensions == 0;
	}

	// ..because derived templates cannot access HCS template params anymore
	static inline level_t GetDimensions() {
		return dimensions;
	}

	// If coord is marked as boundary, retrieve which boundary. (0 = X+, 1 = X-, 2= Y+ ...)
	static coord_t GetBoundaryDirection(coord_t coord) {
		return (coord << 1) >> (HCS_COORD_BITS - dimensions);
	}

	// If a coord is marked as boundary, retrieve the originating coord
	coord_t removeBoundary(coord_t coord) {
		return (coord << (dimensions + 1)) >> (dimensions + 1);
	}

	// Return neighbor for a certain direction. 0=X+, 1=X-, 2=Y+, 3=Y-,...
	// This uses overflow arithmetic, aka the successor formula
	// https://en.wikipedia.org/wiki/Moser%E2%80%93de_Bruijn_sequence
	coord_t getNeighbor(coord_t coord, uint8_t direction) {
		coord_t s_mask = successor_mask[direction];
		coord_t result = 0;
		if (direction & 1) { // negative direction
			result = (coord & s_mask) - 1;
			result &= s_mask;
			result |= (~s_mask) & coord; // bits that have nothing to do with direction stay unchanged
		} else {
			result = (coord | s_mask) + 1;
			result &= ~s_mask;
			result |= s_mask & coord; // bits that have nothing to do with direction stay unchanged
		}
		// Did we hit the boundary?
		if (__count_leading_zeros(coord) == __count_leading_zeros(result))
			return result;

		// Mark result as boundary
		result = coord;
		result |= (coord_t)1 << (HCS_COORD_BITS - 1);
		result |= (coord_t)direction << (HCS_COORD_BITS - 1 - dimensions);

		return result;
	}

	// Alternative approach using unscaled Cartesian coords. On systems with BMI2 instructions,
	// this is on-par with original.
	coord_t getNeighbor2(coord_t coord, uint8_t direction) {
#ifdef __BMI2__
		level_t bit_pos = GetLevelBitPosition(coord);
		coord_t mask = _bzhi_u64(bmi_mask[direction >> 1], bit_pos);
		uint32_t unscaled =  _pext_u64(coord, mask);
		unscaled += direction & 1 ? -1 : 1;

		if (unscaled >= 1U << (bit_pos / dimensions)) {
			coord |= (coord_t)1 << (HCS_COORD_BITS - 1);
			coord |= (coord_t)direction << (HCS_COORD_BITS - 2);
		} else {
			coord &= ~mask; // clear bits for dim while leaving level bits untouched
			coord |=  _pdep_u64(unscaled, mask);
		}
		return coord;
#else

		uint32_t unscaled = getSingleUnscaled(coord, direction >> 1);
		level_t l = GetLevel(coord);
		uint32_t max_coord = (1U << l);
		unscaled += direction & 1 ? -1 : 1;
		if (unscaled >= max_coord) {
			coord |= (coord_t)1 << (HCS_COORD_BITS - 1);
			coord |= (coord_t)direction << (HCS_COORD_BITS - 2);
		} else {
			setSingleUnscaled(coord, l, direction >> 1, unscaled);
		}
		return coord;
#endif
	}

	// Returns a normal vector for the provided direction
	pos_t getDirectionNormal(uint8_t direction) {
		pos_t result {}; // all zero
		result [direction >> 1] = direction & 1 ? -1. : 1.;
		return result;
	}

	data_t getDistance(coord_t coord, uint8_t direction) {
		return (2 * scales[direction]) / data_t(1U << GetLevel(coord));
	}

	// Returns the iteration level of this coordinate. Higher level coordinates carry more information.
	// CAREFUL: The Special bit is not cleared here for performance reasons!
	static inline level_t GetLevel(coord_t coord) {
		return GetLevelBitPosition(coord) / dimensions;
	}

	static inline level_t GetLevelBitPosition(coord_t coord) {
		return (HCS_COORD_BITS - 1 - __count_leading_zeros(coord));
	}

private:
	// Turns the coordinate into an inconsistent state by removing the level-marker bit
	// Returns the position of the level marker.
	static level_t RemoveLevel(coord_t &coord) {
		level_t bit_pos = GetLevelBitPosition(coord);
#ifdef __BMI2__
		coord = _bzhi_u64(coord, bit_pos);
#else
		coord = coord & ((1U << bit_pos ) - 1);
#endif
		return bit_pos;
	}

	// Remove bits for level
	static coord_t RemoveLevel(coord_t coord, level_t level) {
#ifdef __BMI2__
		coord = _bzhi_u64(coord, level * dimensions);
#else
		coord = coord & ((1U << level * dimensions) - 1);
#endif
		return coord;
	}

	// Set the level-marker bit of a (raw) coordinate
	static void SetLevel(coord_t &coord, level_t level) {
		coord_t level_bit = 1U << (level * dimensions);
		coord |= level_bit;
	}

public:

	// Return the closest coordinate at the next lower level.
	static coord_t ReduceLevel(coord_t coord) {
		if (IsBoundary(coord) || coord <= 1)
			return coord;
		coord >>= dimensions;  // remove highest level, marker bit wanders too
		return coord;
	}

	// Increases the level of coord, setting the highest level to a new sub-coordinate (must be between 0 and part_mask)
	static coord_t IncreaseLevel(coord_t coord, uint8_t new_level_coord) {
		if (IsBoundary(coord))
			return coord;
		//assert(new_level_coord < (1U << dimensions));
		return (coord << dimensions) + new_level_coord;
	}


	// Extract the single-level-coord, not checking validity!!
	// The order is reversed here, level 0 is the highest level!
	uint16_t extract(coord_t coord, level_t level) {
		return (coord >> (dimensions * level)) & part_mask;
	}

	pos_t getPosition(coord_t coord) {
		pos_t result = center;
		getPosition(coord, result);
		return result;
	}

	void getPosition(coord_t coord, pos_t &result) {
		if (IsBoundary(coord) || coord <= 1)
			return;
		unscaled_t unscaled = getUnscaled(coord);
		level_t level = GetLevel(coord);
		data_t scale_divisor = 1./ data_t(1 << level);
		for (uint8_t dim = 0; dim < dimensions; dim++)
			result[dim] = scales[dim] * ((data_t)unscaled[dim] * scale_divisor * 2 + scale_divisor) - center[dim] + scales[dim];
	}

	// Returns the coordinate closest to the provided Cart. coordinates for specific level
	coord_t createFromPosition(level_t level, pos_t pos) {
		unscaled_t unscaled;
		data_t scale_divisor = data_t(1 << level);
		for (uint8_t dim = 0; dim < dimensions; dim++)
			unscaled[dim] = floor(((pos[dim] - center[dim]) / (scales[dim] * 2) + 0.5) * scale_divisor);
		return createFromUnscaled(level, unscaled);
	}

	// create a coord from a sub-coordinate list (a sub-coord is between 0 and (2^dimension)-1)
	// Example: coord_t my_lower_left_third_level_coord = createFromList({0,0,0});
	coord_t createFromList(initializer_list<uint8_t> sub_coords) {
		coord_t result = 1;
		for (auto sub : sub_coords)
			result = IncreaseLevel(result, sub);
		return result;
	}


	array<coord_t, 1 << dimensions> getCoeffCoords(coord_t coord) {
		array<coord_t, 1 << dimensions> result;
		uint16_t high_part = extract(coord, 0); 	//
		coord_t origin = ReduceLevel(coord);

		for (uint32_t i = 0; i < (1 << dimensions); i++) {
			coord_t current = origin;
			bool bc_hit = false;
			for (uint32_t j = 0; j < dimensions; j++) {
				if (((i >> j) & 1) == 0)
					continue;

				coord_t prev_current = current;
				current = getNeighbor(current, 2 * j + (((high_part >> j) & 1) ? 0 : 1));
				if (IsBoundary(current)) {
					result[i] = current;
					current = prev_current;
					bc_hit = true;
				}
			}
			if (bc_hit)
				continue;
			result[i] = current;
		}
		return result;
	}

	map<coord_t, data_t> getCoeffs(coord_t coord) {
		map<coord_t, data_t> result;
		// Spawn a rectangle of lower-level coords around missing coord
		// A (hyper)cubical interpolation (2D bi-linear, 3D tri-linear,...) is the best choice,
		// simplexes (triangle, tetrahedron, ...) are not unique in orthogonal spaced coordinates.
		// The neighborhood-search that returns 2^D coordinates that cover our coord is
		// surprisingly straight forward, and the interpolation factors follow the same schema.
		// The originating coord is the one from reducing coord. It is always our closest corner,
		// should therefore get the highest interpolation factor.
		// From there the high_part of coord determines the first D search directions.
		// Hypercube search pattern that surrounds coord from a lower level:
		// 2D: Requires 4 coords (box). The first is _aways_ the level-reduced version of
		//	   coord itself, others are determined by the reduced direction (high_part) of coord.
		// high_part = 0b11 -> X+ Y+ (X+)Y+ <<- SAME ->> (Y+)X+  = 3 neighbors
		// 			   0b00 -> X- Y- (X-)Y- <<- SAME ->> (Y-)X-
		// 			   0b01 -> X+ Y- (X+)Y- <<- SAME ->> (Y-)X+
		// Coords 3D : Box with 8 corners, one is known.
		// 0b101 -> X+ Y- Z+ (X+)Y- (X+)Z+ (Y-)Z+ ((X+)Y-)Z+
		//        The order is not important. Many combinations lead to the same coord.
		//		  This combination follows a bit-order from ordinary counting!
		//		  Three bits for three dimensions, the order is not important, the
		//		  neighborhood direction from high_part is!
		//        X+ Y- Z+
		//        0  0  0     (nothing, the origin point)
		//        0  0  1     Z+
		//        0  1  0     Y-
		//		  0  1  1     Y- -> Z+
		//	      1  0  0     X+
		//		  1  0  1     X+ -> Z+
		//        1  1  0     X+ -> Y-  (the neighbor of X+ in Y- direction)
		//		  1  1  1     X+ -> Y- -> Z+ (the one on the opposite site)
		// Weights 3D:
		//     0 = 0.75, 1 = 0.25
		//
		//        0  0  0  =  0.75³         = 0.4219
		//        0  0  1  =  0.25  * 0.75² = 0.1406
		//        0  1  0  =  0.25  * 0.75² = 0.1406
		//        0  1  1  =  0.25² * 0.75  = 0.0469
		//		  1  0  0  =  0.25  * 0.75² = 0.1406
		//        1  0  1  =  0.25² * 0.75  = 0.0469
		//        1  1  0  =  0.25² * 0.75  = 0.0469
		//        1  1  1  =  0.25³         = 0.0156
		//						TOTAL	    = 1 :)
		// This principle is universal for all dimensions!
		uint16_t high_part = extract(coord, 0); 	//
		coord_t origin = ReduceLevel(coord);

		array<bool, 1 << dimensions> boundary_quench;
		for (uint8_t j = 0; j < dimensions; j++) {
			bool plus = ((high_part >> j) & 1);
			boundary_quench[j] = IsBoundary(getNeighbor(origin, 2 * j + (plus ? 0 : 1)));
		}

		for (uint32_t i = 0; i < (1 << dimensions); i++) {
			coord_t current = origin;
			array<coord_t, dimensions> bc_collector;
			uint32_t bc_collector_count = 0;

			data_t weight = 1;

			for (uint32_t j = 0; j < dimensions; j++)
				if (boundary_quench[j])
					weight *= 0.5;
				else
					weight *= ((i >> j) & 1) ? 0.25 : 0.75;

			for (uint32_t j = 0; j < dimensions; j++) {
				if (((i >> j) & 1) == 0)
					continue;

				coord_t prev_current = current;
				current = getNeighbor(current, 2 * j + (((high_part >> j) & 1) ? 0 : 1));
				if (IsBoundary(current)) {
					bc_collector[bc_collector_count++] = current;
					current = prev_current;
				}
			}
			if (bc_collector_count > 0) {
				for (uint32_t k = 0; k < bc_collector_count; k++)
					result[bc_collector[k]] += weight / bc_collector_count;
				continue;
			}
			result[current] += weight;
		}

		return result;
	}

	// Create "largest" linear coord for provided level
	static coord_t CreateMaxLevel(level_t level) {
		return ((coord_t)1U << (level * dimensions + 1)) - 1;
	}

	// Create "smallest" linear coord for provided level
	static coord_t CreateMinLevel(level_t level) {
		return ((coord_t)1U << (level * dimensions));
	}

	// Inspired by https://github.com/Forceflow/libmorton/
	coord_t createFromUnscaled(level_t level, unscaled_t cart_coord) {
		coord_t result = 0;
#ifdef __BMI2__
		for (uint8_t dim = 0; dim < dimensions; dim++)
			result |=  _pdep_u64(cart_coord[dim], RemoveLevel(bmi_mask[dim], level));
#else
		level_t l = level + 1;
		while (l--) {
			coord_t shift = (coord_t)1 << l;
			for (uint8_t dim = 0; dim < dimensions; dim++)
				result |= (cart_coord[dim] & shift) << ((dimensions - 1) * l + dim);
		}
#endif
		SetLevel(result, level);
		return result;
	}

	// Alters a single unscaled Cartesian component
	void setSingleUnscaled(coord_t &result, level_t level, uint8_t dim, uint32_t unscaled_coord) {
		coord_t mask = RemoveLevel(bmi_mask[dim], level);
		result &= ~mask; // clear bits for dim while leaving level bits untouched

#ifdef __BMI2__
		result |=  _pdep_u64(unscaled_coord, mask);
#else
		level++;
		while (level--) {
			coord_t shift = (coord_t)1 << level;
			result |= (unscaled_coord & shift) << ((dimensions - 1) * level + dim);
		}
#endif
	}

	// Inspired by https://github.com/Forceflow/libmorton/
	uint32_t getSingleUnscaled(coord_t c, uint8_t dim) {
		level_t level = RemoveLevel(c) + dimensions;
#ifdef __BMI2__
		return  _pext_u64(c, bmi_mask[dim]);
#else
		uint32_t result = 0;
		coord_t one = 1;
		while (level -= dimensions) {
			uint8_t shift_selector = level;
			uint8_t shiftback = (dimensions - 1) * level / dimensions;
			result |= (c & (one << (shift_selector + dim))) >> (shiftback + dim);
		}
		return result;
#endif
	}

	// Inspired by https://github.com/Forceflow/libmorton/
	unscaled_t getUnscaled(coord_t c) {
		unscaled_t result {};
		level_t level = RemoveLevel(c) + dimensions;
#ifdef __BMI2__
		for (uint8_t dim = 0; dim < dimensions; dim++)
			result[dim] =  _pext_u64(c, bmi_mask[dim]);
#else
		while (level -= dimensions) {
			uint8_t shift_selector = level - dimensions;
			uint8_t shiftback = (dimensions - 1) * (level - dimensions) / dimensions;
			for (uint8_t dim = 0; dim < dimensions; dim++)
				result[dim] |= (c & (1U << (shift_selector + dim))) >> (shiftback + dim);
		}
#endif
		return result;
	}

	// Increment / Decrement coord in-place. Returns true if a level-transition happened.
	// !! IMPORTANT: The inc/dec routines only work correctly if a valid coord is passed. REASON: isValid() too expensive.
	bool inc(coord_t &coord) {
		level_t l = GetLevel(coord);
		coord++;
		if (coord > level_bounds[l].second) {
			coord = level_bounds[l + 1].first;
			return true;
		}
		return false;
	}

	// Does not decrease below 1
	bool dec(coord_t &coord) {
		if (coord <= 1)
			return true;
		level_t l = GetLevel(coord);
		coord--;
		if (coord < level_bounds[l].first) {
			coord = level_bounds[l - 1].second;
			return true;
		}
		return false;
	}

	// Increment / Decrement coord in-place by the amount of parts. Returns true if a level-transition happened.
	// !! IMPORTANT: The inc/dec routines only work correctly if a valid coord is passed. REASON: isValid() too expensive.
	bool incParts(coord_t &coord) {
		level_t l = GetLevel(coord);
		coord += parts;
		if (coord > level_bounds[l].second) {
			coord = level_bounds[l + 1].first;
			return true;
		}
		return false;
	}

	// Does not decrease below 1
	bool decParts(coord_t &coord) {
		if (coord <= 1)
			return true;
		level_t l = GetLevel(coord);
		coord -= parts;
		if (coord < level_bounds[l].first) {
			coord = level_bounds[l - 1].second;
			return true;
		}
		return false;
	}

	// Turns c into a linear-index that includes the level.
	// HCS coords are almost linear but have gaps between the levels.
	size_t coord2index(coord_t c) {
		level_t l = RemoveLevel(c);
#ifdef __BMI2__
		coord_t mask = _bzhi_u64(bmi_mask[0], l);
#else
		coord_t mask = bmi_mask[0] & ((1U << l) - 1);
#endif
		return c + mask;
	}

	// Turns a linear index back into a HCS coord. SLOW!
	// Should be avoided.
	coord_t index2coord(size_t index) {
		for (level_t l = 1; l < max_level; l++) {
			size_t start_idx = coord2index(level_bounds[l].first);
			size_t end_idx = coord2index(level_bounds[l].second);
			if (index >= start_idx && index <= end_idx) {
				index -= start_idx;
				coord_t result = index;
				SetLevel(result, l);
				return result;
			}
		}
		cout << "R: " << index << endl;
		throw range_error("index2coord oob ");
		//return 0;
	}



	string toString(coord_t coord) {
		stringstream result;
		if (coord == 0)
			return string("(SPECIAL/ZERO)");
		if (coord == 1)
			return string("(CENTER)");

		if (IsBoundary(coord)) {
			coord_t origin = removeBoundary(coord);
			if (IsValid(origin))
				result << "(BOUNDARY: " << (GetBoundaryDirection(coord)) << " ORIGIN : " << toString(origin) << ")";
			else
				result << "(BOUNDARY / INVALID ORIGIN)";
			return result.str();
		}

		if (!IsValid(coord)) {
			return string("(INVALID)");
		}
		int level = GetLevel(coord);
		result << "(" << level << ") [";
		for (int i = 1; i <= level; i++)
			result << extract(coord, level - i) << (i < level ? ", " : "]");
		pos_t pos = getPosition(coord);
		result << " (";
		for (uint8_t d = 0; d < dimensions; d++)
			result << pos[d] << (d == dimensions-1 ? "": ", ");
		result << ")";
		return result.str();
	}

};

};






