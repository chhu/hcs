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
using namespace hcs;

template <typename DTYPE, typename HCSTYPE>
class Field {

 public:
	class iterator; // The public iterator class
	class dual_iterator; // The public iterator class

	Field(char symbol, HCSTYPE hcs_) : symbol(symbol), hcs(hcs_), bracket_behavior(BR_THROW) {
		coeff_up_count = coeff_down_count = 0;
		clear();
		// Create single-value center bucket, the only coordinate that always exists. [0]

		for (auto &bf : boundary)
			bf = nullptr;
		for (bool &bf_prop : boundary_propagate)
			bf_prop = true;
	}

	Field(char symbol) : Field(symbol, HCSTYPE()) {}
	Field(HCSTYPE hcs) : Field('x', hcs) {}
	Field() : Field('x', HCSTYPE()) {}


	// The copy constructor, to make quick copies of the field and its structure
	// Field<??> a = b; Or Field<??> a(b);
	Field(const Field<DTYPE, HCSTYPE> &f) {
		//cout << "FCOPY\n"; // debug hint, this is an expensive op and can happen when you least expect it
		this->symbol = f.symbol;
		this->hcs = f.hcs;
		this->bracket_behavior = f.bracket_behavior;
		this->data = f.data;
		this->tree = f.tree;
		boundary_propagate = f.boundary_propagate;
		for (int i = 0; i < 64; i++)
			this->boundary[i] = boundary_propagate[i] ? f.boundary[i] : nullptr;
		coeff_up_count = coeff_down_count = 0;
	}

	~Field() {}

	// Any other type of Field is a friend.
	template <typename U, typename V>
	friend class Field;

	// The H-coordinate system to operate on
	HCSTYPE	hcs;

	// The boundary functions
	array<function<DTYPE(Field<DTYPE, HCSTYPE> *self, coord_t origin)>, 64> boundary; // max 32 dimensions
	array<bool, 64> boundary_propagate;												  // if this field is copied, is the boundary function copied too?

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

	// The actual data is stored linear to coord for efficiency. Data storage is _not_ sparse!
	vector<DTYPE> data;

	// The tree has the same size as index2coord(data.size()). It is therefore in the wasteful
	// HCS coord space. tree[coord] reveals if coord exists or not.
	vector<bool> tree;

 public:

	// Returns the number of available elements for this field
	size_t nElements() {
		//return count(tree.begin(), tree.end(), true);
		size_t sum = 0;
		for (const auto  e : tree)
		//for (size_t i = 0; i < tree.size(); i++)
			sum += e;
		return sum;
	}

	// Returns the number of top-level elements for this field
	size_t nElementsTop() {
		size_t sum = 0;

		for (auto it = begin(true); it != end(); ++it)
			sum++;
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
    	    	return data[hcs.coord2index(coord)];
    		case BR_NOTHING:
    			return intermediate;
    		}
    	}
    	return data[hcs.coord2index(coord)];
    }

	// Do we have a value for this coord? And if yes, make sure it is in _current
    // A bucket's end coord is its last existing coord
	bool exists(coord_t coord) {
		//	if (hcs.IsBoundary(coord) || tree.size() <= coord)
		//if (tree.size() <= coord)
		//	return false;

		return (coord < tree.size()) ? tree.at(coord) : false;
	}

	// Does not query coefficients, throws if coord does not exist
	DTYPE& getDirect(coord_t coord) {
    	if (!this->exists(coord))
			throw range_error("[]: Coord does not exist");
    	return data[hcs.coord2index(coord)];
	}

	// Returns value for coord, if not present, interpolates.
	// if it is not TLC, return value anyway. To retrieve proper values from non-TLC
	// call propagate() first
	DTYPE get(coord_t coord, bool use_non_top = true) {
		DTYPE result = 0;
		get(coord, result, use_non_top);
		return result;
	}

	void get(coord_t coord, DTYPE& result, bool use_non_top = true) {
		if (hcs.IsBoundary(coord)) {
			uint8_t boundary_index = hcs.GetBoundaryDirection(coord);
			if (boundary[boundary_index] != nullptr)
				result = boundary[boundary_index](this, coord);
			else
				result = 0;
			return;

		}
		if (exists(coord)) {
			if (use_non_top || isTop(coord)) {
				result += data[hcs.coord2index(coord)];
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
							result += boundary[boundary_index](this, b_coord) * (weight / (data_t)boundary_shares.size()); // Ask the provided boundary callback
					}
					continue;
				}

				bool current_exists = exists(current);
				if (!current_exists || (current_exists && !isTop(current) && !use_non_top)) {
					// we either have a non-existent coord or an existing non-top coord that we shall not use.
					/* Caching works but needs to be invalidated every time the field changes...
					auto it = upscale_cache.find(current);
					if (it != upscale_cache.end()) {
						result += it->second * weight;
						continue;
					}
					*/
					DTYPE partial = 0;
					coeff_down_count++;
					get(current, partial, use_non_top);
					//upscale_cache[current] = partial;
					result += partial * weight;
				} else { // current_exists = true in this branch, so _current is valid.
					result += data[hcs.coord2index(current)] * weight;
				}
			}
		}
	}

	// Do coordinates exist in a higher level?
	bool isTop(coord_t coord) {
		return exists(coord) ? !exists(hcs.IncreaseLevel(coord, 0)) : false;
	}


	// Average all non-top coords from top-level
	void propagate(bool max = false) {

		uint32_t parts = hcs.parts;
		data_t inv_parts = 1. / data_t(parts);
		size_t idx = data.size() - 1;
		coord_t c = hcs.index2coord(idx);
		c -= c % parts;
		idx -= idx % parts;

		while (c > parts) {
			DTYPE sum = max ? -INFINITY : 0;
			for (size_t j = idx; j < idx + parts; j++)
				if (max)
					sum = std::max(sum, data[j]);
				 else
					sum += data[j];
			if (!max)
				sum *= inv_parts;
			data[hcs.coord2index(hcs.ReduceLevel(c))] = sum;
			hcs.decParts(c);
			idx -= parts;
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
			//printBucketInfo();
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

				bool current_exists = exists(current);
				if (!current_exists || (current_exists && !isTop(current) && !use_non_top)) {
					// we either have a non-existent coord or an existing non-top coord that we shall not use.
					coeff_map_t partial;
					coeff_down_count++;
					getCoeffs(current, partial, use_non_top, recursion + 1);
					for (auto &coeff : partial)
						coeffs[coeff.first] += coeff.second * weight;
				} else { // current_exists = true in this branch, so _current is valid.
					coeffs[current] += weight;
				}
			}
		}
	}

	// .. and all levels below.
	// This routine DELETES everything in the field and is meant as an initializer.
	// If there are elements present, it throws.
	void createEntireLevel(level_t level) {
		if (data.size() > 2)
			throw range_error("Not empty!");

		coord_t level_end_c = hcs.CreateMaxLevel(level);
		size_t level_end_idx = hcs.coord2index(level_end_c);
		data.resize(level_end_idx + 1, DTYPE(0));
		tree.resize(level_end_c + 1, false);
		for (level_t l = level; l > 0; l--) {
			coord_t level_start = hcs.CreateMinLevel(l);
			coord_t level_end = hcs.CreateMaxLevel(l);
			for (coord_t c = level_start; c <= level_end; c++) {
				tree[c] = true;
			}
		}
	}

	// refine one level up from _existing_ coordinate
	// creates 2^d new coordinates.
	void refineFrom(coord_t coord, bool interpolate_new_values = true) {
		if (!exists(coord))
			throw range_error("refineFrom() Trying to refine from a coord that does not exist!");
		if (!isTop(coord))
			return;
		coord_t lower_corner = hcs.IncreaseLevel(coord, 0);
		coord_t upper_corner = hcs.IncreaseLevel(coord, hcs.part_mask);
		size_t upper_corner_idx = hcs.coord2index(upper_corner);

		if (data.size() <= upper_corner_idx) {
			level_t l = hcs.GetLevel(upper_corner);
			//coord_t max_coord = hcs.CreateMaxLevel(l) + 2;
			//cout << "RESIZE to level: " << l << endl;
			data.resize(upper_corner_idx + 1);
			tree.resize(upper_corner + 1, false);
		}

		vector<DTYPE> interpolated(hcs.parts, data[hcs.coord2index(coord)]);
		if (interpolate_new_values)
			for (coord_t i = 0; i < hcs.parts; i++)
				interpolated[i] = get(lower_corner + i);

		for (coord_t c = lower_corner; c <= upper_corner; c++) {
			tree[c] = true;
			(*this)[c] = interpolated[c - lower_corner];
		}

	}

	// refine up until coord exists
	void refineTo(coord_t coord) {
		coord_t existing = coord;

		// Traverse down until a coord exists (worst-case is 0 or center)
		int i = 0;
		while (!exists(existing)) {	// could be faster
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

		if (isTop(coord))
			return; // Nothing on top

		coord_t lower_corner = hcs.IncreaseLevel(coord, 0);
		coord_t upper_corner = hcs.IncreaseLevel(coord, hcs.part_mask);
		for (coord_t c = lower_corner; c <= upper_corner; c++) {
			if (!isTop(c))
				coarse(c);
			tree[c] = false;
		}
	}

	// Return highest stored coord-level. Could be faster.
	level_t getHighestLevel() {
		coord_t c = tree.size() - 1;
		while (c > 1 && !tree[c])
			c--;
		return hcs.GetLevel(c);
	}


    // Assignment operator requires equal structure, dirty-check with data.size()
	// isTop is not copied because of assumption of equal structure
	Field &operator=(const Field& f){
		//cout << "XCOPY\n";
		assert(("= Operator would alter structure. if this is intended, call takeStructure(x) first!",
				data.size() == f.data.size()));

		data = f.data;

		boundary_propagate = f.boundary_propagate;
		for (int i = 0; i < 64; i++)
			this->boundary[i] = boundary_propagate[i] ? f.boundary[i] : nullptr;
		return *this;
	};

	// Assign f to all elements
    Field<DTYPE, HCSTYPE>& operator=(DTYPE f) {
    	//fill(data.begin(), data.end(), f);
    	for (auto e : *this)
    		e.second = f;
    	return *this;
    }

    // Arithmetic Ops, preserving structure of current refinement.
    // Exampe: a * b keeps sparse structure of a and multiplies with (possible) interpolates from b
    // while b * a keeps sparse structure of b. A generic merge() can specify merged structure and arbitrary ops.

    Field<DTYPE, HCSTYPE> operator-() const { Field<DTYPE, HCSTYPE> result = *this; for (auto e : result) e.second = -e.second; return result;}

    Field<DTYPE, HCSTYPE>& operator*= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second *= const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}
    Field<DTYPE, HCSTYPE>& operator/= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second /= const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}
    Field<DTYPE, HCSTYPE>& operator+= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second += const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}
    Field<DTYPE, HCSTYPE>& operator-= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second -= const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}

    Field<DTYPE, HCSTYPE>& operator*= (const DTYPE& val) { for (auto e : (*this)) e.second *= val; return *this;}
    Field<DTYPE, HCSTYPE>& operator/= (const DTYPE& val) { for (auto e : (*this)) e.second /= val; return *this;}
    Field<DTYPE, HCSTYPE>& operator+= (const DTYPE& val) { for (auto e : (*this)) e.second += val; return *this;}
    Field<DTYPE, HCSTYPE>& operator-= (const DTYPE& val) { for (auto e : (*this)) e.second -= val; return *this;}

	// Clears the field and takes the same coordinate structure as the provided field, without copying their
	// values. The provided field may have a different DTYPE. The newly created coords are initialized with zero.
	template <typename DTYPE2>
	void takeStructure(Field<DTYPE2, HCSTYPE> &f) {
		if (sameStructure(f))
			return;
		tree = f.tree;
		data.resize(hcs.coord2index(tree.size() - 1), DTYPE(0));
	}

	// Tests if the provided field has the same structure.
	// The provided field may have a different DTYPE. The newly created coords are initialized with zero.
	template <typename DTYPE2>
	bool sameStructure(Field<DTYPE2, HCSTYPE> &f) {
		if (f.data.size() != data.size())
			return false;
		auto count = std::inner_product(std::begin(tree), std::end(tree), std::begin(f.tree), 0, std::plus<bool>(), std::equal_to<bool>());
		return count == tree.size();
	}

    // Converts a Field with another DTYPE according to convert function. The structure of "this" remains.
    // The convert function must have a single argument of the foreign DTYPE2 and return DTYPE.
    // Empties "this" first.
    // Calling convention:
    //	target_field.convert< source_data_type[not Field-type!] >(source_field,
    //					[](coord_t, source_field_type) { return target_data_type;});
    // This example turns a "vector" field into a scalar field marking the
    // length of each vector:
	//  ScalarField2 vecmag;
	//  vecmag.convert<Tensor1<data_t, 2> >(v2, [](coord_t c, VectorField2 &source)->data_t {return source.get(c).length();});
	template <typename DTYPE2>
	void convert(Field<DTYPE2, HCSTYPE> &source, function<DTYPE(coord_t, Field<DTYPE2, HCSTYPE> &)> convert_fn) {

		for (auto it = begin(true); it != end(); ++it) {
			coord_t own_coord = (*it).first;
			(*it).second = convert_fn(own_coord, source);
		}

	}

	// Merge 2 fields with possible foreign data type into "this".
	// Arbitrary operations possible through the converter function.
	// The structure of "this" remains.
	// merger function must have 2 arguments of foreign DTYPE2& and return DTYPE.
	// The resulting structure will be the one of f1!
	template <typename DTYPE2>
	void merge(Field<DTYPE2, HCSTYPE> &source1, Field<DTYPE2, HCSTYPE> &source2, function<DTYPE(coord_t, DTYPE2, DTYPE2)> merge_fn) {

		for (auto it = begin(true); it != end(); ++it) {
			coord_t own_coord = (*it).first;
			(*it).second = merge_fn(own_coord, source1.get(own_coord), source2.get(own_coord));
		}
	}

	// Empties all data
	void clear() {
		data.clear();
		data.resize(2, DTYPE(0));
		tree.clear();
		tree.resize(2, false);
		tree[1] = 1;
	}

	// Iterator methods & class
	//iterator dummy = iterator(this);
    iterator begin(bool top_only = false, int only_level = -1) {
        return iterator(this, top_only, only_level);
    }

    iterator end() {	// Just dummy, the begin iterator determines termination
    	return iterator(NULL);
    }


	// C++ goodies, with this operator you can iterate over all existing coords in a field
	class iterator {
	  public:
	    iterator(Field<DTYPE, HCSTYPE>* field, bool top_only = false, int only_level = -1) : current(1), top_start(0), field(field),	only_level(only_level), top_only(top_only) {
	    	if (field == NULL)
	    		return;
	    	hcs = field->hcs;
	    	if (top_only && only_level > 0)
	    		throw range_error("Field iterator can only be top_only or only_level, not both.");
	    	at_end = !field->tree[current];

	    	if (!at_end && top_only && !field->isTop(current)) {
	    		++(*this);
	    	}
	    	if (!at_end && top_only)
	    		top_start = current;

	    	if (!at_end && only_level > 0) {
	    		current = hcs.CreateMinLevel(only_level);
	    		if (!field->exists(current))
	    			at_end = increment();
	    	}
	    }

	    // these three methods form the basis of an iterator for use with
	    // a range-based for loop
	    bool operator!= (const iterator& other) const {
	        return !at_end;
	    }

	    pair<coord_t, DTYPE&> operator* () const {
	       if (at_end)
	    	   throw range_error("Iterator reached end and was queried for value!");
//	       return pair<coord_t, DTYPE&>(current, field->data[]);	// This should not happen... Other containers return garbage
	       const size_t idx = field->hcs.coord2index(current);
	       return pair<coord_t, DTYPE&>(current, field->data[idx]);	// This should not happen... Other containers return garbage
	    }

	    pair<coord_t, DTYPE&>* operator-> () const {
	    	return &(*this);
	    }

	    iterator& operator++ () {
	    	at_end = increment();

	    	if (top_only) {
	    		if (field->isTop(current))
	    			return *this;

	    		bool level_up = at_end ? false : field->tree[current]; // if current exists and is not top, we need to inc level
	    		do {
	    			current = level_up ? hcs.IncreaseLevel(current, 0) : hcs.ReduceLevel(current);
	    		} while (!field->isTop(current));
	    		at_end = top_start == current;
	    		//while (!at_end && !field->isTop(current)) // expensive loop
		    	//	at_end = increment();
	    	} else {
		    	if (current & hcs.part_mask > 0) // current is a multiple of parts
			        return *this;
		    	while (!at_end && !field->tree[current])
		    		at_end = increment();
	    	}
	        return *this;
	    }


	  private:

	    bool at_end, top_only;
	    coord_t top_start;
	    int only_level;
	    coord_t current;
	    Field<DTYPE, HCSTYPE> *field;
	    HCSTYPE hcs;

	    // increment current to next valid coord, including level-jumps and out-of-bounds check.
	    // current may _not_ exist after call
	    // return true if at end
	    bool increment() {
	    	if (hcs.inc(current) && only_level > 0)
	    		return true;

	    	if (current >= field->tree.size()) {
	    		return true;
	    	}
	    	return false;
	    }
	};

	// Iterator methods & class
    dual_iterator begin(Field<DTYPE, HCSTYPE>* field2, bool top_only = false, int only_level = -1) {
    	return dual_iterator(this, field2, iterator(this, top_only, only_level), iterator(field2, top_only, only_level));
    }

	class dual_iterator {
	  public:
	    dual_iterator(Field<DTYPE, HCSTYPE>* field1, Field<DTYPE, HCSTYPE>* field2, iterator if1, iterator if2 ) : field1(field1), field2(field2), if1(if1), if2(if2) {
	    }

	    // these three methods form the basis of an iterator for use with
	    // a range-based for loop
	    bool operator!= (const iterator& other) const {
	        return if1 != field1->end() && if2 != field2->end();
	    }

	    tuple<coord_t, DTYPE&, DTYPE&> operator* () const {
	    	if (!(if1 != field1->end()))
	    	    throw range_error("Iterator reached end and was queried for value!");
	        auto e1 = (*if1);
	        auto e2 = (*if2);
	        if (e1.first != e2.first)
	    	    throw range_error("dual_iterator() called with inconsistent fields!");
	        return tuple<coord_t, DTYPE&, DTYPE&>(e1.first, e1.second, e2.second);
	    }

	    dual_iterator& operator++ () {
	    	++if1;
	    	++if2;
	    	return *this;
	    }

	  private:
	    Field<DTYPE, HCSTYPE> *field1;
	    Field<DTYPE, HCSTYPE> *field2;
	    iterator if1;
	    iterator if2;
	};


};

// Other non-member arithmetic ops
template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator* (const Field<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result *= rhs;
	return result;
};

template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator* (const DTYPE& val, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = rhs;
	result *= val;
	return result;
};

template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator* (const Field<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result *= val;
	return result;
};

template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator/ (const Field<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result /= rhs;
	return result;
};
template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator/ (const DTYPE& val, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = rhs;
	for (auto e : result)
		e.second = val / e.second;
	return result;
}
template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator/ (const Field<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result /= val;
	return result;
}

template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator+ (const Field<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result += rhs;
	return result;
};
template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator+ (const DTYPE& val, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = rhs;
	result += val;
	return result;
}
template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator+ (const Field<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result += val;
	return result;
}

template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator- (const Field<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result -= rhs;
	return result;
};
template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator- (const DTYPE& val, const Field<DTYPE, HCSTYPE>& rhs) {
	Field<DTYPE, HCSTYPE> result = -rhs;
	result += val;
	return result;
}
template <typename DTYPE, typename HCSTYPE> Field<DTYPE, HCSTYPE> operator- (const Field<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
	Field<DTYPE, HCSTYPE> result = lhs;
	result -= val;
	return result;
}

