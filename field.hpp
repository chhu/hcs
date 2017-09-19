#pragma once
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

    Field(HCSTYPE hcs_) : hcs(hcs_), bracket_behavior(BR_INTERP) {
        for (auto &bf : boundary)
            bf = nullptr;
        for (bool &bf_prop : boundary_propagate)
            bf_prop = true;
    }

    Field() : Field(HCSTYPE()) {}

    Field(level_t level) : Field(HCSTYPE()) {
    	createEntireLevel(level);
    }

    virtual ~Field() {}

    // Any other type of Field is a friend.
    template <typename U, typename V>
    friend class Field;

    // The H-coordinate system to operate on
    HCSTYPE	hcs;

    // The boundary functions
    array<function<DTYPE(Field<DTYPE, HCSTYPE> *self, coord_t origin)>, 64> boundary; // max 32 dimensions
    array<bool, 64> boundary_propagate;												  // if this field is copied, is the boundary function copied too?

    // If a value is accessed via [], and if that value does not exist:
    //   BR_THROW: throws range_error, slow if it happens often.
    //   BR_REFINE: brings requested coord into existence via refineToCoord(), might be very slow
    //   BR_INTERP: useful if you only read from the coord. The intermediate is filled with the interpolated value (via get()).
    //				writing to the returned reference just sets the intermediate.
    //   BR_NOTHING: Just return a reference to the intermediate. Fastest version.
    //				Set intermediate to a value that marks non-existence and check return...
    //				writing to the returned reference just sets the intermediate.
    enum { BR_THROW, BR_INTERP, BR_NOTHING, BR_REFINE } bracket_behavior;

    //  Used as reference for the [] operator if coord does not exist, see above
    DTYPE		intermediate;

    //  The type to store a list of coords and their coefficients.
    //  It is a map instead of a vector because of unique coord elimination.
    typedef map<coord_t, data_t> coeff_map_t;


public:

    // This is the (strictly forward) iterator class that needs to be implemented
    class CustomIterator {
    public:
        CustomIterator() : at_end(true), currentValPtr(NULL), currentCoord(0) {}
        virtual void increment() { cerr << "CI: INC CALLED\n"; throw bad_function_call();};
        virtual pair<coord_t, DTYPE&>* getCurrentPairPtr() { cerr << "No CustomIterator implemented\n"; throw bad_function_call();};
        virtual CustomIterator* clone() {  cerr << "CI: CLONE CALLED\n"; throw bad_function_call();};
        bool at_end;
        DTYPE* currentValPtr;
        coord_t currentCoord;
    };

    // NEVER overwrite this class.
    class Iterator {
    private:
        CustomIterator* ci;
    public:
        Iterator() : ci() {}
        Iterator(CustomIterator* ci) : ci(ci) {}
        Iterator(Iterator const& right) : ci(right.ci->clone()) {}

        ~Iterator() { delete ci; }

        Iterator& operator=(Iterator const& right) {
            delete ci;
            ci = right.ci->clone();
            return *this;
        }
        // these three methods form the basis of an iterator for use with
        // a range-based for loop
        bool operator!= (const Iterator& other) const {
            return !ci->at_end;
        }

        //pair<coord_t, DTYPE&> operator* () const { return *ci->getCurrentPairPtr();};
        pair<coord_t, DTYPE&> operator* () const { return pair<coord_t, DTYPE&>(ci->currentCoord, *ci->currentValPtr);};    // 33% faster!
        pair<coord_t, DTYPE&>* operator-> () const { return ci->getCurrentPairPtr();};

        Iterator& operator++ () { ci->increment(); return *this;};
    };

    // Iterator methods & class
    virtual Iterator begin(bool top_only = false, int only_level = -1) = 0;
    virtual Iterator end() = 0;

    // Returns the number of available elements for this field
    virtual size_t nElements() = 0;

    // Returns the number of top-level elements for this field
    virtual size_t nElementsTop() = 0;

    // Read-write access to existing coords. For (probably) non-existing, use get() and retrieve interpolated values.
    // the value set to bracket_behavior applies.
    virtual DTYPE& operator[](coord_t coord) = 0;

    // Do we have a value for this coord? And if yes, make sure it is in _current
    // A bucket's end coord is its last existing coord
    virtual bool exists(coord_t coord) = 0;

    // Does not query coefficients, throws if coord does not exist
    virtual DTYPE& getDirect(coord_t coord) = 0;

    // Returns value for coord, if not present, interpolates.
    // if it is not TLC, return value anyway. To retrieve proper values from non-TLC
    // call propagate() first
    virtual DTYPE get(coord_t coord, bool use_non_top = true) = 0;

    void get(coord_t coord, DTYPE& result, bool use_non_top = true) {
        if (hcs.IsBoundary(coord)) {
            uint8_t boundary_index = hcs.GetBoundaryDirection(coord);
            if (boundary[boundary_index] != nullptr)
                result = boundary[boundary_index](this, coord);
            else {
                result = 0;
                get(hcs.removeBoundary(coord), result, use_non_top);
            }
            return;
        }
        if (exists(coord)) {
            if (use_non_top || isTop(coord)) {
                result += getDirect(coord);
                return;
            } else {
                for (uint16_t direction = 0; direction < hcs.parts; direction++) {
                    //coeff_up_count++;
                    DTYPE partial = 0;
                    //getCoeffs(hcs.IncreaseLevel(coord, direction), partial, use_non_top, recursion + 1);
                    get(hcs.IncreaseLevel(coord, direction), partial, use_non_top);
                    partial /= (data_t)hcs.parts;
                    result += partial;
                }
            }
        } else {
        	auto coeffs = hcs.getCoeffs(coord);
        	for (auto coeff : coeffs) {
        		coord_t current = coeff.first;
        		data_t weight = coeff.second;
                bool current_exists = exists(current);
                if (!current_exists || (current_exists && !isTop(current) && !use_non_top)) {
                    // we either have a non-existent coord or an existing non-top coord that we shall not use.
                    DTYPE partial = 0;
                    get(current, partial, use_non_top);
                    result += partial * weight;
                } else { // current_exists = true in this branch, so _current is valid.
                    result += getDirect(current) * weight;
                }
        	}
        }
    }

    void correct_neumann(pair<coord_t, data_t>* coeffs) {
        int n_coeffs = 1 << hcs.GetDimensions();
        int non_neumann_bc = 0;

        for (int i = 0; i < n_coeffs; i++) {
            if (coeffs[i].second == 0)
                continue;
            if (hcs.IsBoundary(coeffs[i].first))
                if (boundary[hcs.GetBoundaryDirection(coeffs[i].first)] == nullptr) {
                    coord_t wob = hcs.removeBoundary(coeffs[i].first);
                    for (int j = 0; j < n_coeffs; j++) {
                        if (coeffs[j].first == wob && coeffs[j].second != 0)
                            coeffs[i].second = coeffs[j].second;
                    }
                    coeffs[i].first = wob;
                } else
                    non_neumann_bc++;
        }
        data_t total = 0;
        for (int i = 0; i < n_coeffs; i++)
            total += coeffs[i].second;
        if (total < 1 && non_neumann_bc > 0) {
            total = (1 - total) / data_t(non_neumann_bc);
            for (int i = 0; i < n_coeffs; i++)
                if (hcs.IsBoundary(coeffs[i].first))
                    coeffs[i].second += total;

        }

            //if (nc > 0 && dc > 0 && nc != dc) {

        //}
    }
    // Do coordinates exist in a higher level?
    virtual bool isTop(coord_t coord) = 0;	// Average all non-top coords from top-level

    // Propagates values down from top-level to lowest level by averaging them.
    // If there would be a reverse iterator, a generic algorithm would be possible here...
    virtual void propagate() = 0;

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
            exit(1);
        }
        if (exists(coord)) {
            if (isTop(coord) || use_non_top) {
                coeffs[coord] = 1.;
                return;
            } else {
                for (uint16_t direction = 0; direction < hcs.parts; direction++) {
                    coeff_map_t partial;
                    getCoeffs(hcs.IncreaseLevel(coord, direction), partial, use_non_top, recursion + 1);
                    for (auto &coeff : partial)
                        coeff.second /= hcs.parts;
                    coeffs.insert(partial.begin(), partial.end());
                }
            }
        } else { // coord does not exist
            // ask HCS for underlying coeffs
            auto sub_coeffs = hcs.getCoeffs(coord);
            for (auto sub_coeff : sub_coeffs) {
                coord_t &current = sub_coeff.first;
                data_t &weight = sub_coeff.second;
                if (weight == 0)
                    continue;
                bool current_exists = exists(current);
                if (!current_exists || (current_exists && !isTop(current) && !use_non_top)) {
                    // we either have a non-existent coord or an existing non-top coord that we shall not use.
                    coeff_map_t partial;
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
    virtual void createEntireLevel(level_t level) = 0;

    // Return highest stored coord-level
    virtual level_t getHighestLevel() = 0;

    // Assignment operator requires equal structure, dirty-check with data.size()
    // isTop is not copied because of assumption of equal structure
    //virtual Field &operator=(const Field& f) = 0;

    virtual Field<DTYPE, HCSTYPE>& operator=(DTYPE f) {
        for (auto e : (*this))
            e.second = f;

        return *this;
    }

    // Arithmetic Ops, preserving structure of current refinement.
    // Exampe: a * b keeps sparse structure of a and multiplies with (possible) interpolates from b
    // while b * a keeps sparse structure of b. A generic merge() can specify merged structure and arbitrary ops.

    //virtual Field<DTYPE, HCSTYPE> operator-() const = 0;//{ Field<DTYPE, HCSTYPE> result = *this; for (auto e : result) e.second = -e.second; return result;}

    virtual Field<DTYPE, HCSTYPE>& operator*= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second *= const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}
    virtual Field<DTYPE, HCSTYPE>& operator/= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second /= const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}
    virtual Field<DTYPE, HCSTYPE>& operator+= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second += const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}
    virtual Field<DTYPE, HCSTYPE>& operator-= (const Field<DTYPE, HCSTYPE>& rhs) { for (auto e : (*this)) e.second -= const_cast<Field<DTYPE, HCSTYPE>*>(&rhs)->get(e.first); return *this;}

    virtual Field<DTYPE, HCSTYPE>& operator*= (const DTYPE& val) { for (auto e : (*this)) e.second *= val; return *this;}
    virtual Field<DTYPE, HCSTYPE>& operator/= (const DTYPE& val) { for (auto e : (*this)) e.second /= val; return *this;}
    virtual Field<DTYPE, HCSTYPE>& operator+= (const DTYPE& val) { for (auto e : (*this)) e.second += val; return *this;}
    virtual Field<DTYPE, HCSTYPE>& operator-= (const DTYPE& val) { for (auto e : (*this)) e.second -= val; return *this;}

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
    virtual void clear() = 0;


};

// Other non-member arithmetic ops
// serve as template for overriding, result and LHS must be of inherited class, RHS can stay generic field as its iterator / copy constructpr is not needed
/*
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
 */

