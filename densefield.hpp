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
class DenseField : public Field<DTYPE, HCSTYPE> {

public:

    DenseField(HCSTYPE hcs_) : Field<DTYPE, HCSTYPE>(hcs_), max_level(1), max_coord(1) {
        clear();
        hcs = hcs_;
    }

    DenseField() : DenseField(HCSTYPE()) {}


    // The copy constructor, to make quick copies of the field and its structure
    // Field<??> a = b; Or Field<??> a(b);
    DenseField(const DenseField<DTYPE, HCSTYPE> &f) {
        //cout << "FCOPY\n"; // debug hint, this is an expensive op and can happen when you least expect it
        this->max_level = f.max_level;
        this->max_coord = f.max_coord;
        this->hcs = f.hcs;
        hcs = this->hcs;
        this->bracket_behavior = f.bracket_behavior;
        this->data = f.data;
        this->boundary_propagate = f.boundary_propagate;
        for (int i = 0; i < 64; i++)
            this->boundary[i] = this->boundary_propagate[i] ? f.boundary[i] : nullptr;
    }

    ~DenseField() {}

    // Any other type of Field is a friend.
    template <typename U, typename V>
    friend class DenseField;
    // Any other type of Field is a friend.
    template <typename U, typename V>
    friend class Field;

private:

    // The actual data is stored linear to coord for efficiency. Data storage is _not_ sparse!
    vector<DTYPE> data;

    // Indicates the current level in data.
    level_t max_level;
    coord_t max_coord;

public:
    HCSTYPE hcs;

    // C++ goodies, with this operator you can iterate over all existing coords in a field
    class DenseIterator : public Field<DTYPE, HCSTYPE>::CustomIterator {
    public:
        DenseIterator(DenseField<DTYPE, HCSTYPE>* field, bool top_only = false, int only_level = -1) : current(1), field(field), only_level(only_level), top_only(top_only), end_coord(1) {
            if (field == NULL)
                return;
            hcs = field->hcs;
            if (top_only && only_level > 0)
                throw range_error("Field iterator can only be top_only or only_level, not both.");
            this->at_end = (only_level > field->max_level);
            if (this->at_end)
                return;
            end_coord = hcs.CreateMaxLevel(field->max_level);
            if (only_level > 0) {
                current = hcs.CreateMinLevel(only_level);
                end_coord = hcs.CreateMaxLevel(only_level);
            } else if (top_only)
                current = hcs.CreateMinLevel(field->max_level);
            else
                current = 1;
        }

        virtual pair<coord_t, DTYPE&> getCurrentPair() {
            if (this->at_end)
                throw range_error("Iterator reached end and was queried for value!");

            const size_t idx = field->hcs.coord2index(current);
            return pair<coord_t, DTYPE&>(current, field->data[idx]);
        }

        virtual void increment() {
            if (hcs.inc(current))
                this->at_end = current > end_coord;
        }

        DenseIterator* clone() {
            DenseIterator* result = new DenseIterator(field, this->top_only, only_level);
            result->current = current;
            result->end_coord = end_coord;
            return result;
        }

    private:
        DenseField<DTYPE, HCSTYPE>* field;
        coord_t current, end_coord;
        int only_level;
        bool top_only;
        HCSTYPE hcs;
    };

    // Iterator methods & class
    //iterator dummy = iterator(this);
    typename Field<DTYPE, HCSTYPE>::Iterator begin(bool top_only = false, int only_level = -1) {
        return typename Field<DTYPE, HCSTYPE>::Iterator(new DenseIterator(this, top_only, only_level));
    }

    typename Field<DTYPE, HCSTYPE>::Iterator end() {	// Just dummy, the begin iterator determines termination
        return NULL;
    }


    // Returns the number of available elements for this field
    size_t nElements() {
        return this->data.size();
    }

    // Returns the number of top-level elements for this field
    size_t nElementsTop() {
        return hcs.CreateMaxLevel(max_level) - hcs.CreateMinLevel(max_level);
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
            throw bad_function_call();

            case Field<DTYPE,HCSTYPE>::BR_NOTHING:
            return this->intermediate;
            }
        }
        return data[hcs.coord2index(coord)];
    }

    bool exists(coord_t coord) {
        return coord <= max_coord;
    }

    // Does not query coefficients, throws if coord does not exist
    DTYPE& getDirect(coord_t coord) {
        //if (!this->exists(coord))
         //   throw range_error("[]: Coord does not exist");
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
        return hcs.GetLevel(coord) == max_level;
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

        max_level = level;
        max_coord = hcs.CreateMaxLevel(level);
        size_t level_end_idx = hcs.coord2index(max_coord);

        data.resize(level_end_idx + 1, DTYPE(0));
    }

    // Return highest stored coord-level. Could be faster.
    level_t getHighestLevel() {
        return max_level;
    }


    // Assignment operator requires equal structure, dirty-check with data.size()
    // isTop is not copied because of assumption of equal structure
    DenseField &operator=(const DenseField& f){
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
        this->operator =(dynamic_cast<DenseField<DTYPE, HCSTYPE> >(f));
        return *this;
    };

    DenseField &operator=(const DTYPE& f){
        fill(data.begin(), data.end(), f);
        //Field<DTYPE,HCSTYPE>::operator =(f);
        return *this;
    }

    // Arithmetic Ops, preserving structure of current refinement.
    // Exampe: a * b keeps sparse structure of a and multiplies with (possible) interpolates from b
    // while b * a keeps sparse structure of b. A generic merge() can specify merged structure and arbitrary ops.
    DenseField<DTYPE, HCSTYPE> operator-() const { DenseField<DTYPE, HCSTYPE> result = *this; for (auto e : result) e.second = -e.second; return result;}

    // Clears the field and takes the same coordinate structure as the provided field, without copying their
    // values. The provided field may have a different DTYPE. The newly created coords are initialized with zero.
    template <typename DTYPE2>
    void takeStructure(DenseField<DTYPE2, HCSTYPE> &f) {
        if (sameStructure(f))
            return;
        data.resize(f.data.size(), DTYPE(0));
        max_level = f.max_level;
        max_coord = f.max_coord;
    }

    // Tests if the provided field has the same structure.
    // The provided field may have a different DTYPE. The newly created coords are initialized with zero.
    template <typename DTYPE2>
    bool sameStructure(const DenseField<DTYPE2, HCSTYPE> &f) {
        return f.data.size() == data.size();
    }

    // Empties all data
    void clear() {
        data.clear();
        //data.resize(2, DTYPE(0));
        max_level = 0;
        max_coord = 0;
    }

    Field<DTYPE, HCSTYPE>& operator*= (const DenseField<DTYPE, HCSTYPE>& rhs) {
        if (sameStructure(rhs)) {
            for (size_t i = 0; i < data.size(); i++)
                data[i] *= rhs.data[i];
        } else
            return Field<DTYPE,HCSTYPE>::operator *=(rhs);
        return *this;
    }

    Field<DTYPE, HCSTYPE>& operator*= (const DTYPE& rhs) {
        for (size_t i = 0; i < data.size(); i++)
            data[i] *= rhs;
        return *this;
    }


    Field<DTYPE, HCSTYPE>& operator+= (const DenseField<DTYPE, HCSTYPE>& rhs) {
        if (sameStructure(rhs)) {
            for (size_t i = 0; i < data.size(); i++)
                data[i] += rhs.data[i];
        } else
            return Field<DTYPE,HCSTYPE>::operator +=(rhs);
        return *this;
    }

    Field<DTYPE, HCSTYPE>& operator+= (const Field<DTYPE, HCSTYPE>& rhs) {
        return Field<DTYPE,HCSTYPE>::operator +=(rhs);
    }

    Field<DTYPE, HCSTYPE>& operator+= (const DTYPE& rhs) {
        for (size_t i = 0; i < data.size(); i++)
            data[i] += rhs;
        return *this;
    }

    Field<DTYPE, HCSTYPE>& operator-= (const DenseField<DTYPE, HCSTYPE>& rhs) {
        if (sameStructure(rhs)) {
            for (size_t i = 0; i < data.size(); i++)
                data[i] -= rhs.data[i];
        } else
            return Field<DTYPE,HCSTYPE>::operator -=(rhs);
        return *this;
    }

    Field<DTYPE, HCSTYPE>& operator-= (const Field<DTYPE, HCSTYPE>& rhs) {
        return Field<DTYPE,HCSTYPE>::operator -=(rhs);
    }

    Field<DTYPE, HCSTYPE>& operator-= (const DTYPE& rhs) {
        for (size_t i = 0; i < data.size(); i++)
            data[i] -= rhs;
        return *this;
    }


    Field<DTYPE, HCSTYPE>& operator/= (const DenseField<DTYPE, HCSTYPE>& rhs) {
        if (sameStructure(rhs)) {
            for (size_t i = 0; i < data.size(); i++)
                data[i] /= rhs.data[i];
        } else
            return Field<DTYPE,HCSTYPE>::operator /=(rhs);
        return *this;
    }
    Field<DTYPE, HCSTYPE>& operator/= (const DTYPE& rhs) {
        for (size_t i = 0; i < data.size(); i++)
            data[i] /= rhs;
        return *this;
    }


};

// Other non-member arithmetic ops
template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator* (const DenseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result *= rhs;
    return result;
};

template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator* (const DenseField<DTYPE, HCSTYPE>& lhs, const DenseField<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result *= rhs;
    return result;
};


template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator* (const DTYPE& val, const DenseField<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = rhs;
    result *= val;
    return result;
};

template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator* (const DenseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result *= val;
    return result;
};

template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator/ (const DenseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result /= rhs;
    return result;
};
template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator/ (const DTYPE& val, const DenseField<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = rhs;
    for (auto e : result)
        e.second = val / e.second;
    return result;
}
template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator/ (const DenseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result /= val;
    return result;
}

template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator+ (const DenseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result += rhs;
    return result;
};

template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator+ (const DenseField<DTYPE, HCSTYPE>& lhs, const DenseField<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result += rhs;
    return result;
};

template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator+ (const DTYPE& val, const DenseField<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = rhs;
    result += val;
    return result;
}
template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator+ (const DenseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result += val;
    return result;
}

template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator- (const DenseField<DTYPE, HCSTYPE>& lhs, const Field<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result -= rhs;
    return result;
};
template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator- (const DTYPE& val, const DenseField<DTYPE, HCSTYPE>& rhs) {
    DenseField<DTYPE, HCSTYPE> result = -rhs;
    result += val;
    return result;
}
template <typename DTYPE, typename HCSTYPE> DenseField<DTYPE, HCSTYPE> operator- (const DenseField<DTYPE, HCSTYPE>& lhs, const DTYPE& val) {
    DenseField<DTYPE, HCSTYPE> result = lhs;
    result -= val;
    return result;
}

