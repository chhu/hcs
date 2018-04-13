/*
 * hcs-config.h
 *
 *  Created on: Feb 2, 2018
 *      Author: C Huettig
 *
 *  Configuration for the HCS class. This file will only be included by hcs.hpp.
 *  Theoretically this would be possible via templates too, but these kind of definitions likely
 *  do not change within a project.
 */

#ifndef HCS_CONFIG_INC_
#define HCS_CONFIG_INC_

// Configured to store a H coordinate in 64 bit. Depending on max used level and dimension, this can be 32 bit.
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


// A type used to carry separated level information, should be able to count to HCS_COORD_BITS,
// so all unsigned integer types will do.
typedef uint16_t level_t;

// Data precision
typedef double data_t;



#endif /* HCS_CONFIG_INC_ */
