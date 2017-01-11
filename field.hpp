/*
 * field.hpp
 *
 * Sparse storage in HCS
 *
 * This is a sparse storage class for the H coordinate system.
 * - dedicated refinement / coarsening
 * - only complete "H"s exist
 * - lower-level coords always exist, but top-level get marked as such. (Top-Level-Coordinate TLC)
 * - iterator class that allows fast iteration over all top-level or all existing coords or
 * 	 existing coords of a specific level
 * - bi-linear interpolation of non-existing coords, providing coefficients for TLC
 * - Arbitrary data type that needs to support some basic arithmetic
 * - Field supplies basic arithmetic operators
 * - A bracket operator for coordinates is implemented, with adjustable behavior for non-existing coords.
 * - The performance of exists() relies on STL's map::lower_bound O(log) complexity
 * - The center coordinate (0) always exists
 * - boundary conditions can be implemented as lambdas
 *
 */

using namespace std;


template <typename DTYPE, typename HCSTYPE>
class Field {

	// Forward declaration of sub-classes
private:
	class Bucket;	// The storage class (private)
 public:
	class iterator; // The public iterator class

	Field(char symbol, HCSTYPE hcs_) : symbol(symbol), _current(NULL), hcs(hcs_), bracket_behavior(BR_THROW) {
		coeff_up_count = coeff_down_count = 0;
		// Create single-value center bucket, the only coordinate that always exists. [0]
		data[0] = new Bucket(0, 0);
		data[0]->setTop(0, true);
		for (auto &bf : boundary)
			bf = nullptr;
	}

	Field(char symbol) : Field(symbol, HCSTYPE()) {}
	Field(HCSTYPE hcs) : Field('x', hcs) {}
	Field() : Field('x', HCSTYPE()) {}


	// The copy constructor, to make quick copies of the field and its structure
	// Field<??> a = b; Or Field<??> a(b);
	Field(const Field<DTYPE, HCSTYPE> &f) {
		cout << "FCOPY\n"; // debug hint, this is an expensive op and can happen when you least expect it
		this->symbol = f.symbol;
		this->hcs = f.hcs;
		this->bracket_behavior = f.bracket_behavior;
		this->data = f.data;
		this->boundary = f.boundary;
		coeff_up_count = coeff_down_count = 0;
		// The buckets are pointers, so in order to not get a reference to the values, we need to copy separately.
		for (auto & bucket : data) {
			bucket.second = new Bucket(*(bucket.second)); // Calls implicit copy constructor of Bucket.
		}
		_current = NULL;
	}

	// Because we "new"d Buckets, we need to release them.
	~Field() {
		for (auto e : data)
			delete e.second;
	}

	 // Any other type of Field is a friend.
	 template <typename U, typename V>
	 friend class Field;

	// The H-coordinate system to operate on
	HCSTYPE	hcs;

	// The boundary functions
	array<function<DTYPE(coord_t origin)>, 64> boundary; // max 32

	/*
	 * Open properties, without getter / setter, feel free to modify them at any time.
	 */
	char		symbol; // Single character symbol, like 'T' to distinguish

	// If a value is accessed via [], and if that value does not exist:
	//   BR_THROW: throws range_error, slow if it happens often.
	//   BR_REFINE: brings requested coord into existence via refineToCoord(), might be very slow
	//   BR_INTERP: useful if you only read from the coord. The intermediate is filled with the interpolated value (via get()).
	//				writing to the returned reference just sets the intermediate.
	//   BR_NOTHING: Just return a reference to the intermediate. Fastest version.
	//				Set intermediate to a value that marks non-existence and check return...
	//				writing to the returned reference just sets the intermediate.
	enum { BR_THROW, BR_REFINE, BR_INTERP, BR_NOTHING } bracket_behavior;

	//  Used as reference for the [] operator if coord does not exist, see above
	DTYPE		intermediate;

	//  The type to store a list of coords and their coefficients.
	//  It is a map instead of a vector because of unique coord elimination.
	typedef map<coord_t, data_t> coeff_map_t;


 private:
	// The actual data and useful typedefs.
	typedef typename map<coord_t, Bucket*, greater<coord_t> >::iterator map_iter_t;
	typedef typename map<coord_t, Bucket*>::iterator map_iter_rev_t;
	typedef map<coord_t, Bucket*, greater<coord_t> > map_t;

	map_t 		data; 				// re-arrange key sort so we can use lower_bound().


	// The last successful bucket of a exists() query.
	// Saves a lot of calls to map.lower_bound() which is expensive
	Bucket*	_current;


 public:

	// Returns the number of available elements for this field
	size_t nElements() {
		size_t sum = 0;
		for (auto const & kv : data)
			sum += kv.second->data.size();
		return sum;
	}

	// Returns the number of top-level elements for this field
	size_t nElementsTop() {
		size_t sum = 0;
		for (auto const & kv : data)
			sum += count(kv.second->top.begin(), kv.second->top.end(), true);
		return sum;
	}

	// Read-write access to existing coords. For (probably) non-existing, use get() and retrieve interpolated values.
	// the value set to bracket_behavior applies.
    DTYPE& operator[](coord_t coord) {
    	if (!this->exists(coord)) {
    		switch (bracket_behavior) {
    		case BR_THROW:
    			throw range_error("[]: Coord does not exist");
    		case BR_INTERP:
    			intermediate = get(coord);
    			return intermediate;
    		case BR_REFINE:
    			refineTo(coord);
    	    	return this->_current->get(coord);
    		case BR_NOTHING:
    			return intermediate;
    		}
    	}
    	return this->_current->get(coord);
    }

	// Do we have a value for this coord? And if yes, make sure it is in _current
    // A bucket's end coord is its last existing coord
	bool exists(coord_t coord) 	{
		if (this->_current != NULL && coord >= _current->start && coord <= _current->end) {
			return true;
		}
		if (hcs.IsBoundary(coord))
			return false;
		map_iter_t result = data.lower_bound(coord);
		if (result == data.end() || result->second->end < coord)
			return false;
		this->_current = result->second;
		return true;
	}

	// Does not query coefficients, throws if coord does not exist
	DTYPE& getDirect(coord_t coord) {
    	if (!this->exists(coord))
			throw range_error("[]: Coord does not exist");
    	return this->_current->get(coord);
	}

	// Returns value for coord, if not present, interpolates.
	// if it is not TLC, return value anyway. To retrieve proper values from non-TLC
	// call propagate() first
	DTYPE get(coord_t coord) {
		if (exists(coord)) {
			return _current->get(coord);
		}
		DTYPE result = 0;
		get(coord, result, true);
		return result;

		coeff_map_t coeffs;
		getCoeffs(coord, coeffs, true);
		for (auto coeff : coeffs) {
			coord_t c = coeff.first;
			if (hcs.IsBoundary(c)) {
				uint8_t boundary_index = hcs.GetBoundaryDirection(c);
				if (boundary[boundary_index] != nullptr) {
					result += boundary[boundary_index](c) * coeff.second; // Ask the provided boundary callback
				}	// If no BC provided, assume zero, do nothing
			} else {
				result += getDirect(c) * coeff.second;
			}
		}
		return result;
	}

	map<coord_t, DTYPE> upscale_cache;
	void get(coord_t coord, DTYPE& result, bool use_non_top = true) {
		if (hcs.IsBoundary(coord)) {
			uint8_t boundary_index = hcs.GetBoundaryDirection(coord);
			if (boundary[boundary_index] != nullptr)
				result = boundary[boundary_index](coord);
			else
				result = 0;
			return;

		}
		if (exists(coord)) {
			if (use_non_top || isTop(coord)) {
				result += _current->get(coord);
				return;
			} else {
				for (uint16_t direction = 0; direction < hcs.parts; direction++) {
					coeff_up_count++;
					DTYPE partial = 0;
					//getCoeffs(hcs.IncreaseLevel(coord, direction), partial, use_non_top, recursion + 1);
					get(hcs.IncreaseLevel(coord, direction), partial, use_non_top);
					partial /= (data_t)hcs.parts;
					result += partial;
				}
			}
		} else {

			uint16_t high_part = hcs.extract(coord, 0);
			coord_t origin = hcs.ReduceLevel(coord);

			array<bool, 64> boundary_quench;

			for (uint8_t j = 0; j < hcs.GetDimensions(); j++) {
				bool plus = ((high_part >> j) & 1);
				boundary_quench[j] = hcs.IsBoundary(hcs.getNeighbor(origin, 2 * j + (plus ? 0 : 1)));
			}


			for (uint8_t i = 0; i <= hcs.part_mask; i++) {
				coord_t current = origin;

				data_t weight = 1;

				for (uint8_t j = 0; j < hcs.GetDimensions(); j++)
					weight *= boundary_quench[j] ? 0.5 : (((i >> j) & 1) ? 0.25 : 0.75);

				set<coord_t> boundary_shares;

				for (uint8_t j = 0; j < hcs.GetDimensions(); j++) {
					if (((i >> j) & 1) == 0)
						continue;

					coord_t prev_current = current;
					current = hcs.getNeighbor(current, 2 * j + (((high_part >> j) & 1) ? 0 : 1));
					if (hcs.IsBoundary(current)) {
						boundary_shares.insert(current);
						current = prev_current;
					}
				}


				// get coeffs for current

				if (!boundary_shares.empty()) {
					for (auto b_coord : boundary_shares) {
						uint8_t boundary_index = hcs.GetBoundaryDirection(b_coord);
						if (boundary[boundary_index] != nullptr)
							result += boundary[boundary_index](b_coord) * (weight / (data_t)boundary_shares.size()); // Ask the provided boundary callback
					}
					continue;
				}
				if (exists(current)) {
					result += _current->get(current) * weight;
				} else {
					auto it = upscale_cache.find(current);
					if (it != upscale_cache.end()) {
						result += it->second * weight;
						continue;
					}
					DTYPE partial = 0;
					coeff_down_count++;
					get(current, partial, use_non_top);
					upscale_cache[current] = partial;
					result += partial * weight;
				}
			}
		}

	}

	// Do coordinates exist in a higher level?
	bool isTop(coord_t coord) {
		if (!exists(coord))
			throw range_error("isTop coord does not exist!");
		return this->_current->isTop(coord);
	}

	// Average all non-top coords from top-level
	void propagate() {
		// Use the <greater> sorting from our data map that will deliver top-level coords first
		// We don't store the averaged values immediately to avoid excessive lower_bound() lookups
		vector<tuple<coord_t, DTYPE> > lower_level_cache;
		level_t highest = getHighestLevel();

		for (auto& entry : data) {
			level_t current_level = hcs.GetLevel(entry.first);
			if (current_level < highest) {
				// We have reached the next lower level. Write cache.
				for (auto& centry : lower_level_cache)
					getDirect(std::get<0>(centry)) = std::get<1>(centry);
				lower_level_cache.clear();
				highest = current_level;
			}
			if (highest == 0) // 0-Bucket has only one coord that should have been filled.
				break;
			// A Bucket must have a multiple of hcs.parts
			Bucket *b = entry.second;
			for (coord_t c = b->start; c <= b->end; c += hcs.parts) {
				DTYPE total = 0;
				for (int j = 0; j < hcs.parts; j++)
					total += b->get(c + j);
				total /= hcs.parts;
				lower_level_cache.push_back(make_tuple(hcs.ReduceLevel(c), total));
			}
		}

	}


	size_t coeff_up_count, coeff_down_count;
	// Return interpolation coeffs and their associated >existing< coords.
	// The first value of the pair is the coefficient, always >0 and <=1.
	// If the coord exists the returning vector will be of size 1 and first=1., second=coord.
	// use_non_top = true uses coefficients from existing, but not top-level coordinates. Use propagate() first
	// to set non-top values to their averaged versions from top-level values.
	// never use recursion parameter, its purely internal to protect the stack.
	void getCoeffs(const coord_t coord, coeff_map_t &coeffs, bool use_non_top = true, int recursion = 0) {
		if (hcs.IsBoundary(coord)) {
			coeffs[coord] = 1.;
			return;
		}
		if (recursion > hcs.max_level) {
			cout << "RECURSION LIMIT REACHED (" << hcs.max_level << ") coord: " << hcs.toString(coord) << endl;
			printBucketInfo();
			exit(1);
		}
		if (exists(coord)) {
			if (isTop(coord) || use_non_top) {
				coeffs[coord] = 1.;
				return;
			} else {
				for (uint16_t direction = 0; direction < hcs.parts; direction++) {
					coeff_up_count++;
					coeff_map_t partial;
					getCoeffs(hcs.IncreaseLevel(coord, direction), partial, use_non_top, recursion + 1);
					for (auto &coeff : partial)
						coeff.second /= hcs.parts;
						//std::get<1>(coeff) /= hcs.parts;
					coeffs.insert(partial.begin(), partial.end());
				}
			}
		} else {

			//typename HCSTYPE::neighbor_t ne;
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
			//     0 = 0.75, 1=0.25
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
			uint16_t high_part = hcs.extract(coord, 0); 	//
			coord_t origin = hcs.ReduceLevel(coord);
			//hcs.getNeighbors(origin, ne);
			array<bool, 64> boundary_quench;


			for (uint8_t j = 0; j < hcs.GetDimensions(); j++) {
				bool plus = ((high_part >> j) & 1);
				boundary_quench[j] = hcs.IsBoundary(hcs.getNeighbor(origin, 2 * j + (plus ? 0 : 1)));
			}


			array<tuple<coord_t, data_t, int>, 64> collection;
			for (uint8_t i = 0; i <= hcs.part_mask; i++) {
				coord_t current = origin;
				std::get<2>(collection[i]) = 0;
				int boundaries_involved = 0;
				vector<coord_t> bc_collector;

				data_t weight = 1;

				for (uint8_t j = 0; j < hcs.GetDimensions(); j++)
					weight *= boundary_quench[j] ? 0.5 : (((i >> j) & 1) ? 0.25 : 0.75);


				for (uint8_t j = 0; j < hcs.GetDimensions(); j++) {
					if (((i >> j) & 1) == 0)
						continue;

					coord_t prev_current = current;
					current = hcs.getNeighbor(current, 2 * j + (((high_part >> j) & 1) ? 0 : 1));
					if (hcs.IsBoundary(current)) {
						bc_collector.push_back(current);
						current = prev_current;
					}
				}

				if (!bc_collector.empty()) {
					for (auto bcc : bc_collector)
						coeffs[bcc] += weight / bc_collector.size();
					continue;
				}

				if (exists(current)) {
					coeffs[current] += weight;
				} else {
					coeff_map_t partial;
					coeff_down_count++;
					getCoeffs(current, partial, use_non_top, recursion + 1);
					for (auto &coeff : partial)
						coeffs[coeff.first] += coeff.second * weight;
				}
			}
		}
	}

	// .. and all levels below.
	// This routine DELETES everything in the field and is meant as an initializer.
	// If there are elements present, it throws.
	void createEntireLevel(level_t level) {
		if (data.size() > 1)
			throw range_error("Not empty!");
		data[0]->setTop(0, false);
		for (level_t l = 1; l <= level; l++) {
			coord_t level_start = 0;
			coord_t level_end = ((coord_t)1 << (l * hcs.GetDimensions())) - 1;
			hcs.SetLevel(level_start, l);
			hcs.SetLevel(level_end, l);
			Bucket* bucket = new Bucket(level_start, level_end);
			data[level_start] = bucket;
			fill(bucket->top.begin(), bucket->top.end(), l == level);
		}
	}

	// refine one level up from _existing_ coordinate
	// creates 2^d new coordinates.
	void refineFrom(coord_t coord) {
		if (!exists(coord))
			throw range_error("refineFrom() Trying to refine from a coord that does not exist!");
		Bucket* coord_bucket = _current;  // bucket that holds coord
		coord_t lower_corner = hcs.IncreaseLevel(coord, 0);
		coord_t upper_corner = hcs.IncreaseLevel(coord, hcs.part_mask);
		if (exists(lower_corner))
			return;	// ? nothing to do...
		Bucket* bucket = new Bucket(lower_corner, upper_corner);
		data[lower_corner] = bucket;
		fill(bucket->top.begin(), bucket->top.end(), true);	// Mark as top
		fill(bucket->data.begin(), bucket->data.end(), coord_bucket->get(coord));	// Set values from orig coord
		// Now the original coord is not top anymore...
		coord_bucket->setTop(coord, false);
		_current = bucket;	// grant immediate access to new coords
	}

	// refine up until coord exists
	void refineTo(coord_t coord) {
		coord_t existing = coord;

		// Traverse down until a coord exists (worst-case is 0 or center)
		int i = 0;
		while (!exists(existing)) {
			existing = hcs.ReduceLevel(existing);
			i++;
		}
		// Refine upwards
		while (i-- > 0) {
			refineFrom(existing);
			existing = hcs.IncreaseLevel(existing, hcs.extract(coord, i));
		}
	}

	// Remove all coords on higher level above coord
	void coarse(coord_t coord) {
		if (!exists(coord))
			return;

		if (_current->isTop(coord))
			return; // Nothing on top

		Bucket* coord_bucket = _current;

		coord_t first_up = hcs.IncreaseLevel(coord, 0);

		if (!exists(first_up))
			throw range_error("Inconsistency detected!");

		bool part_range_top = true;
		for (int i = 0; i < hcs.parts; i++)
			part_range_top &= _current->isTop(first_up + i);
		if (part_range_top) {
			removeCoords(first_up);
			coord_bucket->setTop(coord, true);
		} else {
			for (int8_t i = hcs.part_mask; i >= 0; i--) // reverse coarsening leads to tail-shrinks instead of expensive head-shrinks.
				coarse(first_up + i);
			coarse(coord);
		}
	}

	// look for buckets that can be combined.
	// Returns amount of combined buckets.
	size_t optimize() {
		size_t result = 0;
		Bucket *last_bucket = data.find((coord_t)0)->second;
		assert(last_bucket != NULL);
		for (auto b = data.begin(); b != data.end();) {
			Bucket *bucket = b->second;
			if (last_bucket->start == bucket->end+1) {
				result++;
				bucket->end = last_bucket->end;
				bucket->data.insert(bucket->data.end(), last_bucket->data.begin(), last_bucket->data.end());
				bucket->top.insert(bucket->top.end(), last_bucket->top.begin(), last_bucket->top.end());
				assert(bucket->data.size() == (bucket->end - bucket->start + 1));
				b = data.erase(--b);
				delete last_bucket;
			} else {
				last_bucket = bucket;
				++b;
			}
		}
		_current = NULL;
		return result;
	}

	// Return highest stored coord-level
	level_t getHighestLevel() {
		return hcs.GetLevel(data.begin()->second->start); // map's sort order is "greater", so highest-level bucket is first.
	}

	// Iterator methods & class
    iterator begin(bool top_only = false, int only_level = -1) {
        return iterator(this, top_only, only_level);
    }

    iterator end() {	// Just dummy, the begin iterator determines termination
    	return iterator(this);
    }

    // Arithmetic Ops, preserving structure of current refinement.
    // Exampe: a * b keeps sparse structure of a and multiplies with (possible) interpolates from b
    // while b * a keeps sparse structure of b. A generic merge() can specify merged structure and arbitrary ops.

    // Assignment operator requires equal structure, dirty-check with data.size()
	Field &operator=(const Field& f){
		cout << "XCOPY\n";
		assert(data.size() == f.data.size());
		auto iter_this = data.begin();
		auto iter_f = f.data.begin();
		while (iter_this != data.end()) {
			Bucket *b_this = iter_this->second;
			Bucket *b_f = iter_f->second;
			b_this->data = b_f->data;
			++iter_this;
			++iter_f;
		}
		_current = NULL;
		return *this;
	};

    Field<DTYPE, HCSTYPE>& operator=(DTYPE f) {
    	for (auto e : *this)
    		e.second = f;
    	return *this;
    }
    Field<DTYPE, HCSTYPE> operator+(Field<DTYPE, HCSTYPE> &f)  {
    	Field<DTYPE, HCSTYPE> result = *this;
    	result += f;
		return result;
	} // add

    Field<DTYPE, HCSTYPE>& operator+=(Field<DTYPE, HCSTYPE>& f) {
    	f.upscale_cache.clear();
    	for (auto e : *this)
    		e.second += f.get(e.first);
    	f.upscale_cache.clear();
    	return *this;
	} // add

    Field<DTYPE, HCSTYPE> operator+(DTYPE f) {
    	Field<DTYPE, HCSTYPE> result = *this;
    	result += f;
		return result;
	} // add

    Field<DTYPE, HCSTYPE>& operator+=(DTYPE f) {
    	for (auto e : *this)
    		e.second += f;
    	return *this;
	} // add

    Field<DTYPE, HCSTYPE> operator-(Field<DTYPE, HCSTYPE> &f)  {
    	Field<DTYPE, HCSTYPE> result = *this;
    	result -= f;
		return result;
	} // add

    Field<DTYPE, HCSTYPE>& operator-=(Field<DTYPE, HCSTYPE>& f) {
    	f.upscale_cache.clear();
    	for (auto e : *this)
    		e.second -= f.get(e.first);
    	f.upscale_cache.clear();
    	return *this;
	} // add

    Field<DTYPE, HCSTYPE> operator*(Field<DTYPE, HCSTYPE> &f)  {
    	Field<DTYPE, HCSTYPE> result = *this;
    	result *= f;
		return result;
	} // add

    Field<DTYPE, HCSTYPE>& operator*=(Field<DTYPE, HCSTYPE>& f) {
    	f.upscale_cache.clear();
    	for (auto e : *this) {
    		e.second *= f.get(e.first);
    	}
    	f.upscale_cache.clear();
    	return *this;
	} // add

    Field<DTYPE, HCSTYPE> operator*(DTYPE f) {
    	Field<DTYPE, HCSTYPE> result = *this;
    	result *= f;
		return result;
	} // add

    Field<DTYPE, HCSTYPE>& operator*=(DTYPE f) {
    	for (auto e : *this) {
    		e.second *= f;
    	}
    	return *this;
	} // add

    const Field<DTYPE, HCSTYPE> operator/(Field<DTYPE, HCSTYPE> &f) const {
    	Field<DTYPE, HCSTYPE> result = *this;
    	result /= f;
		return result;
	} // add

    Field<DTYPE, HCSTYPE>& operator/=(Field<DTYPE, HCSTYPE>& f) {
    	f.upscale_cache.clear();
    	for (auto e : *this)
    		e.second /= f.get(e.first);
    	f.upscale_cache.clear();

    	return *this;
	} // add

    // Converts a Field with another DTYPE according to convert function.
    // The convert function must have a single argument of the foreign DTYPE2 and return DTYPE.
    // Empties "this" first.
    // This example turns a "vector" field into a scalar field marking the
    // length of each vector:
	//  ScalarField2 v2x('x', &h2);
	//  v2x.convert<Tensor1<data_t, 2> >(v2, [](Tensor1<data_t, 2> t2)->data_t {return t2.length();});
	template <typename DTYPE2>
	void convert(const Field<DTYPE2, HCSTYPE> &f, function<DTYPE(coord_t, DTYPE2&)> convert) {
		clear();
		for (auto e : f.data) {
			auto *b = e.second;
			Bucket *bn = new Bucket(b->start, b->end);
			bn->top = b->top;
			for (coord_t c = bn->start; c <= bn->end; c++)
				bn->get(c) = convert(c, b->get(c));
			data[b->start] = bn;
		}
	}

	// Merge 2 fields with possible foreign data type into "this".
	// Arbitrary operations possible through the converter function.
	// clear() is called first so everything in "this" will be deleted.
	// merger function must have 2 arguments of foreign DTYPE2& and return DTYPE.
	// The resulting structure will be the one of f1!
	template <typename DTYPE2>
	void merge(const Field<DTYPE2, HCSTYPE> &f1, Field<DTYPE2, HCSTYPE> &f2, function<DTYPE(coord_t, DTYPE2&, DTYPE2)> merger) {
		clear();
		for (auto e : f1.data) {
			auto *b = e.second;
			Bucket *bn = new Bucket(b->start, b->end);
			bn->top = b->top;
			for (coord_t c = bn->start; c <= bn->end; c++)
				bn->get(c) = merger(c, b->get(c), f2.get(c));
			data[b->start] = bn;
		}
	}

	// Empties all data
	void clear() {
		for (auto e : data)
			delete e.second;
		data.clear();
		data[0] = new Bucket(0, 0);
		data[0]->setTop(0, true);
	}

	// C++ goodies, with this operator you can iterate over all existing coords in a field
	class iterator {
	  public:
	    iterator(Field<DTYPE, HCSTYPE>* field, bool top_only = true, int only_level = -1) : global_index(0), field(field),
		bucket(NULL), bucket_index(0), only_level(only_level), top_only(top_only) {
	    	map_iter = field->data.begin();
	    	if (only_level >= 0) {
	    		while (map_iter != field->data.end()) {
	    			++map_iter;
		    		coord_t start = map_iter->first;
		    		if (field->hcs.GetLevel(start) == only_level)
		    			break;
	    		}
	    	}
	    	at_end = !(map_iter != field->data.end());
	    	if (!at_end) {
	    		bucket = map_iter->second;

	    		// Skip eventual non-tops
	    		if (top_only)
	    			while (!at_end && !bucket->top[bucket_index])
	    				increment();

	    	}
	    }

	    // these three methods form the basis of an iterator for use with
	    // a range-based for loop
	    bool operator!= (const iterator& other) const {
	        return !at_end;
	    }

	    // this method must be defined after the definition of IntVector
	    // since it needs to use it
	    //    pair<coord_t, array<data_t, components> > operator* () const {
	    //       return make_pair<coord_t, array<data_t, components> >(0, [1]);
	    // }
	    pair<coord_t, DTYPE&> operator* () const {
	       if (at_end)
	    	   throw range_error("Iterator reached end and was queried for value!");
	       return pair<coord_t, DTYPE&>(bucket->start + bucket_index, bucket->data[bucket_index]);	// This should not happen... Other containers return garbage
	    }

	    iterator& operator++ () {
	    	if (top_only) {
	    		do {
	    			increment();
	    		} while (!at_end && !bucket->top[bucket_index]);
	    	} else
	    		increment();
	    	global_index++;
	        return *this;
	    }


	  private:
	    void increment() {
	    	bucket_index++;
	    	if (bucket_index >= bucket->data.size()) {
	    		++map_iter;
	    		if (!(map_iter != field->data.end())) {
	    			this->at_end = true;
	    			return;
	    		}
	    		bucket = map_iter->second;
	    		bucket_index = 0;
	    		if (only_level >= 0 && field->hcs.GetLevel(bucket->start) < only_level) {
	    			at_end = true;
	    			return;
	    		}
	    	}
	    }

	    bool at_end, top_only;
	    map_iter_t map_iter;
	    size_t global_index, bucket_index;
	    int only_level;
	    Bucket *bucket;
	    Field<DTYPE, HCSTYPE> *field;
	};

    /*
     * IO routines
     *
     * format:
     * 		"HCF" 	- 3 char magic
     * 		"0"   	- 1 char version
     * 				- coordinate info:
     * 		<uint8>	- HCS_COORD_BITS; info about coord format (must be multiple of 8)
     * 		<uint8> - HCS_LEVEL_BIT
     * 		<uint8> - dimensions
     * 	D x	<double>- center position Cartesian position of coord 0.
     * 	D x <double>- scales
     * 				- FIELD INFO:
     * 		<char>  - DTYPE identifier; "f" = float "i" = integer "u" = unsigned int
     * 		<uint8> - DTYPE bits per component
     * 		<uint8> - components (1 = scalar, 3= vector (or RGB) etc.
     * 		<uint64>- N coord : value pairs
     * 		Following a combination of N
     * 		    coord             value
     * 		<HCS_COORD_BITS> + <DTYPE bits> * components
     */
    // Expects binary stream
    bool save(ofstream &stream) {

    	return true;
    }

    bool load(ofstream &stream) {
    	return true;
    }

	void printBucketInfo() {

		for (auto b : data) {
			size_t n_top = count(b.second->top.begin(), b.second->top.end(), true);
			cout << "Bucket: N = " << b.second->data.size() << " N_top = " << n_top << " Start: " << hcs.toString(b.second->start) << " End: " << hcs.toString(b.second->end) <<endl;
		}
	}
  private:

	// unconditionally remove without checking hierarchy, start until start + part_mask get thrown away.
	void removeCoords(coord_t start) {
		start = start & (~hcs.part_mask); // make sure first sub-coord is zero
		coord_t end = start | hcs.part_mask;

		map_iter_t result = data.lower_bound(start);

		if (result == data.end() || result->second->end < start)
			return;
		Bucket *b = result->second;
		coord_t b_start = b->start;
		coord_t b_end = b->end;
		if (start < b_start || end > b_end)
			throw range_error("Inconsistency while coarsing");

		if (b_start == start && b_end == end) {	// Bucket matches range
			// Throw the whole bucket away
			delete b;
			data.erase(result);
			_current = NULL;
		} else if (b_end == end && b_start < start) { // Bucket needs shrinking (tail clip)
			b->end = start - 1;
			b->data.resize(b->end - b->start + 1);
			b->top.resize(b->end - b->start + 1);
		} else if (b_end > end && b_start == start) { // Bucket needs shrinking (head clip)
			_current = new Bucket(end + 1, b_end, b);
			delete b;
			data.erase(result);
			data[end+1] = _current;

		} else {	// coords are within a bucket, need to split...
			_current = new Bucket(end + 1, b_end, b);
			data[end+1] = _current;

			b->end = start - 1;
			b->data.resize(b->end - b->start + 1);
			b->top.resize(b->end - b->start + 1);
		}
	}

	// A simple storage container that exposes an STL vector and associates it with a coord-range.
	// data.size() == top.size() == (end-start+1)

	class Bucket {
		Bucket(coord_t _start, coord_t _end) : start(_start), end(_end) {
			this->data.resize(_end-_start + 1);
			this->top.resize(_end-_start + 1);
		}

		// Assumes that from_bucket's range is bigger than own start-end and copies values into own vectors
		Bucket(coord_t _start, coord_t _end, Bucket* from_bucket) : Bucket(_start, _end) {
			for (coord_t i = start; i <= end; i++) {
				size_t bindex = from_bucket->index(i); // foreign index for coord
				size_t oindex = index(i);  // own index for coord
				this->data[oindex] = from_bucket->data[bindex];
				this->top[oindex] = from_bucket->top[bindex];
			}
		}

		template<typename U, typename V>
		friend class Field;	// Field can access private

		coord_t			start, end;
		vector<DTYPE> 	data;
		vector<char>	top;

		size_t index(coord_t coord) {
			assert(coord >= start && coord <= end);
			return coord - this->start;
		}

		bool isTop(coord_t coord) {
			return top[coord - start];
		}

		DTYPE& get(coord_t coord) {
			size_t idx = index(coord);
			//assert(idx < data.size());
			return data[idx];
		}

		void setTop(coord_t coord, bool is_top) {
			top[coord - start] = is_top;
		}
	};

};


