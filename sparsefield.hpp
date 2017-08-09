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
 * - The center coordinate (1) always exists
 * - boundary conditions can be implemented as lambdas
 *
 */

using namespace std;
using namespace hcs;

template <typename DTYPE, typename HCSTYPE>
class SparseField : public Field<DTYPE, HCSTYPE> {

public:

    SparseField(HCSTYPE hcs_) : Field<DTYPE, HCSTYPE>(hcs_) {
        clear();
        hcs = hcs_;
    }

    SparseField() : SparseField(HCSTYPE()) {}


    // The copy constructor, to make quick copies of the field and its structure
    // Field<??> a = b; Or Field<??> a(b);
    SparseField(const SparseField<DTYPE, HCSTYPE> &f) {
        //cout << "FCOPY\n"; // debug hint, this is an expensive op and can happen when you least expect it
        this->hcs = f.hcs;
        hcs = this->hcs;
        this->bracket_behavior = f.bracket_behavior;
        this->data = f.data;
        this->tree = f.tree;
        this->boundary_propagate = f.boundary_propagate;
        for (int i = 0; i < 64; i++)
            this->boundary[i] = this->boundary_propagate[i] ? f.boundary[i] : nullptr;
    }

    ~SparseField() {}

    // Any other type of Field is a friend.
    template <typename U, typename V>
    friend class SparseField;
    // Any other type of Field is a friend.
    template <typename U, typename V>
    friend class Field;

private:

    // The actual data is stored linear to coord for efficiency. Data storage is _not_ sparse!
    vector<DTYPE> data;

    // The tree has the same size as index2coord(data.size()). It is therefore in the wasteful
    // HCS coord space. tree[coord] reveals if coord exists or not.
    vector<bool> tree;

public:
    HCSTYPE hcs;

    // C++ goodies, with this operator you can iterate over all existing coords in a field
    class SparseIterator : public Field<DTYPE, HCSTYPE>::CustomIterator {
    public:
        SparseIterator(SparseField<DTYPE, HCSTYPE>* field, bool top_only = false, int only_level = -1) : current(1), top_start(0), field(field), only_level(only_level), top_only(top_only) {
            if (field == NULL)
                return;
            hcs = field->hcs;
            if (top_only && only_level > 0)
                throw range_error("Field iterator can only be top_only or only_level, not both.");
            this->at_end = !field->tree[current];

            if (!this->at_end && top_only && !field->isTop(current)) {
                increment();
            }
            if (!this->at_end && top_only)
                top_start = current;

            if (!this->at_end && only_level > 0) {
                current = hcs.CreateMinLevel(only_level);
                if (!field->exists(current))
                    this->at_end = increment2();
            }
        }

        virtual pair<coord_t, DTYPE&> getCurrentPair() {
            if (this->at_end)
                throw range_error("Iterator reached end and was queried for value!");
            //	       return pair<coord_t, DTYPE&>(current, field->data[]);	// This should not happen... Other containers return garbage
            const size_t idx = field->hcs.coord2index(current);
            return pair<coord_t, DTYPE&>(current, field->data[idx]);	// This should not happen... Other containers return garbage
        }

        virtual void increment() {
            this->at_end = increment2();

            if (top_only) {
                if (field->isTop(current))
                    return;// *this;

                bool level_up = this->at_end ? false : field->tree[current]; // if current exists and is not top, we need to inc level
                do {
                    current = level_up ? hcs.IncreaseLevel(current, 0) : hcs.ReduceLevel(current);
                } while (!field->isTop(current));
                this->at_end = top_start == current;
                //while (!this->at_end && !field->isTop(current)) // expensive loop
                //	this->at_end = increment();
            } else {
                if (current & hcs.part_mask > 0) // current is a multiple of parts
                    return;// *this;
                while (!this->at_end && !field->tree[current])
                    this->at_end = increment2();
            }
            return;// *this;
        }

        SparseIterator* clone() {
            SparseIterator* result = new SparseIterator(field, this->top_only, only_level);
            result->current = current;
            result->top_start = top_start;
            return result;
        }

    private:
        coord_t top_start;
        SparseField<DTYPE, HCSTYPE>* field;
        coord_t current;
        int only_level;
        bool top_only;
        HCSTYPE hcs;

        // increment current to next valid coord, including level-jumps and out-of-bounds check.
        // current may _not_ exist after call
        // return true if at end
        bool increment2() {
            if (hcs.inc(current) && only_level > 0)
                return true;

            if (current >= field->tree.size()) {
                return true;
            }
            return false;
        }
    };

    // Iterator methods & class
    //iterator dummy = iterator(this);
    typename Field<DTYPE, HCSTYPE>::Iterator begin(bool top_only = false, int only_level = -1) {
        return typename Field<DTYPE, HCSTYPE>::Iterator(new SparseIterator(this, top_only, only_level));
    }

    typename Field<DTYPE, HCSTYPE>::Iterator end() {	// Just dummy, the begin iterator determines termination
        return NULL;//typename Field<DTYPE, HCSTYPE>::Iterator(new SparseIterator(NULL));//new SparseIterator(this, top_only, only_level));
    }


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
            switch (this->bracket_behavior) {
            case Field<DTYPE,HCSTYPE>::BR_THROW:
            throw range_error("[]: Coord does not exist");
            case Field<DTYPE,HCSTYPE>::BR_INTERP:
            this->intermediate = get(coord);
            return this->intermediate;
            case Field<DTYPE,HCSTYPE>::BR_REFINE:
            refineTo(coord);
            return data[hcs.coord2index(coord)];
            case Field<DTYPE,HCSTYPE>::BR_NOTHING:
            return this->intermediate;
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
        Field<DTYPE, HCSTYPE>::get(coord, result, use_non_top);
        return result;
    }

    // Do coordinates exist in a higher level?
    bool isTop(coord_t coord) {
        return exists(coord) ? !exists(hcs.IncreaseLevel(coord, 0)) : false;
    }


    // Average all non-top coords from top-level
    void propagate() {

        uint32_t parts = hcs.parts;
        data_t inv_parts = 1. / data_t(parts);
        size_t idx = data.size() - 1;
        coord_t c = hcs.index2coord(idx);
        c -= c % parts;
        idx -= idx % parts;

        while (c > parts) {
            DTYPE sum = 0;
            for (size_t j = idx; j < idx + parts; j++)
                sum += data[j];
            sum *= inv_parts;
            data[hcs.coord2index(hcs.ReduceLevel(c))] = sum;
            hcs.decParts(c);
            idx -= parts;
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
    SparseField &operator=(const SparseField& f){
        //cout << "XCOPY\n";
        assert(("= Operator would alter structure. if this is intended, call takeStructure(x) first!",
                data.size() == f.data.size()));

        data = f.data;

        this->boundary_propagate = f.boundary_propagate;
        for (int i = 0; i < 64; i++)
            this->boundary[i] = this->boundary_propagate[i] ? f.boundary[i] : nullptr;
        return *this;
    };

    // Assignment operator requires equal structure, dirty-check with data.size()
    // isTop is not copied because of assumption of equal structure
    Field<DTYPE, HCSTYPE> &operator=(const Field<DTYPE, HCSTYPE>& f){
        this->operator =(dynamic_cast<SparseField<DTYPE, HCSTYPE> >(f));
        return *this;
    };

    SparseField &operator=(const DTYPE& f){
        Field<DTYPE,HCSTYPE>::operator =(f);
        return *this;
    }

    // Arithmetic Ops, preserving structure of current refinement.
    // Exampe: a * b keeps sparse structure of a and multiplies with (possible) interpolates from b
    // while b * a keeps sparse structure of b. A generic merge() can specify merged structure and arbitrary ops.
    SparseField<DTYPE, HCSTYPE> operator-() const { SparseField<DTYPE, HCSTYPE> result = *this; for (auto e : result) e.second = -e.second; return result;}

    // Clears the field and takes the same coordinate structure as the provided field, without copying their
    // values. The provided field may have a different DTYPE. The newly created coords are initialized with zero.
    template <typename DTYPE2>
    void takeStructure(SparseField<DTYPE2, HCSTYPE> &f) {
        if (sameStructure(f))
            return;
        tree = f.tree;
        data.resize(hcs.coord2index(tree.size() - 1), DTYPE(0));
    }

    // Tests if the provided field has the same structure.
    // The provided field may have a different DTYPE. The newly created coords are initialized with zero.
    template <typename DTYPE2>
    bool sameStructure(SparseField<DTYPE2, HCSTYPE> &f) {
        if (f.data.size() != data.size())
            return false;
        auto count = std::inner_product(std::begin(tree), std::end(tree), std::begin(f.tree), 0, std::plus<bool>(), std::equal_to<bool>());
        return count == tree.size();
    }

    // Empties all data
    void clear() {
        data.clear();
        data.resize(2, DTYPE(0));
        tree.clear();
        tree.resize(2, false);
        tree[1] = 1;
    }

};

// Other non-member arithmetic ops
template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator* (const SparseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result *= rhs;
    return result;
};

template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator* (const DTYPE& val, const SparseField<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = rhs;
    result *= val;
    return result;
};

template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator* (const SparseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result *= val;
    return result;
};

template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator/ (const SparseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result /= rhs;
    return result;
};
template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator/ (const DTYPE& val, const SparseField<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = rhs;
    for (auto e : result)
        e.second = val / e.second;
    return result;
}
template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator/ (const SparseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result /= val;
    return result;
}

template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator+ (const SparseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result += rhs;
    return result;
};
template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator+ (const DTYPE& val, const SparseField<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = rhs;
    result += val;
    return result;
}
template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator+ (const SparseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result += val;
    return result;
}

template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator- (const SparseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result -= rhs;
    return result;
};
template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator- (const DTYPE& val, const SparseField<DTYPE, HCSTYPE>& rhs) {
    SparseField<DTYPE, HCSTYPE> result = -rhs;
    result += val;
    return result;
}
template <typename DTYPE, typename HCSTYPE> SparseField<DTYPE, HCSTYPE> operator- (const SparseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    SparseField<DTYPE, HCSTYPE> result = lhs;
    result -= val;
    return result;
}

