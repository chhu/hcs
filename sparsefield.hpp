/*
 * sparsefield.hpp
 *
 * Sparse storage in HCS
 *
 * This is a sparse storage class for the H coordinate system.
 * - dedicated refinement / coarsening
 * - only complete "H"s exist
 * - lower-level coords always exist, but top-level get marked as such. (Top-Level-Coordinate TLC)
 * - iterator class that allows fast iteration over all top-level or all existing coords or
 *   existing coords of a specific level
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
class SparseField  : public Field<DTYPE, HCSTYPE> {

    // Forward declaration of sub-classes
private:
    class Bucket;   // The storage class (private)
 public:

    SparseField(HCSTYPE hcs_) : Field<DTYPE, HCSTYPE>(hcs_), _current(NULL), _level_current{NULL}, hcs(hcs_) {
        // Create single-value center bucket, the only coordinate that always exists. [0]
        data[1] = new Bucket(1, 1);
        data[1]->setTop(1, true);
    }

    SparseField() : SparseField(HCSTYPE()) {}

    SparseField(level_t level) : SparseField(HCSTYPE()) {
    	createEntireLevel(level);
    }



    // The copy constructor, to make quick copies of the field and its structure
    // Field<??> a = b; Or Field<??> a(b);
    SparseField(const SparseField<DTYPE, HCSTYPE> &f) {
        //cout << "FCOPY\n"; // debug hint, this is an expensive op and can happen when you least expect it
        this->hcs = f.hcs;
        this->bracket_behavior = f.bracket_behavior;
        this->data = f.data;
        this->boundary_propagate = f.boundary_propagate;
        for (int i = 0; i < 64; i++)
            this->boundary[i] = this->boundary_propagate[i] ? f.boundary[i] : nullptr;
        // The buckets are pointers, so in order to not get a reference to the values, we need to copy separately.
        for (auto & bucket : data) {
            bucket.second = new Bucket(*(bucket.second)); // Calls implicit copy constructor of Bucket.
        }
        _current = NULL;
        for (int i = 0; i < 64; i++)
            _level_current[i] = NULL;
    }

    // Because we "new"d Buckets, we need to release them.
    ~SparseField() {
        for (auto e : data)
            delete e.second;
    }

     // Any other type of Field is a friend.
     template <typename U, typename V>
     friend class SparseField;

    // The H-coordinate system to operate on
    HCSTYPE hcs;

    //  Used as reference for the [] operator if coord does not exist, see above
    DTYPE       intermediate;

 private:
    // The actual data and useful typedefs.
    //typedef typename map<coord_t, Bucket*, less<coord_t> >::iterator map_iter_rev_t;
    typedef map<coord_t, Bucket*, greater<coord_t> > map_t;
    typedef typename map_t::iterator map_iter_t;
    typedef typename map_t::reverse_iterator map_iter_t_rev;

    map_t       data;               // re-arrange key sort so we can use lower_bound().


    // The last successful bucket of a exists() query.
    // Saves a lot of calls to map.lower_bound() which is expensive
    Bucket* _current;
    Bucket* _level_current[64];

 public:
    class SparseIterator : public Field<DTYPE, HCSTYPE>::CustomIterator {
    public:
        SparseIterator(SparseField<DTYPE, HCSTYPE>* field, bool top_only = false, int only_level = -1) : field(field), bucket(NULL), bucket_index(0), only_level(only_level), top_only(top_only), current_pair(0, intermediate) {
            map_iter = field->data.begin();
            if (only_level >= 0) {
                while (map_iter != field->data.end()) {
                    coord_t start = map_iter->first;
                    if (field->hcs.GetLevel(start) == only_level)
                        break;
                    ++map_iter;
                }
            }
            this->at_end = !(map_iter != field->data.end());
            if (!this->at_end) {
                bucket = map_iter->second;

                // Skip eventual non-tops
                if (top_only)
                    while (!this->at_end && !bucket->top[bucket_index])
                        increment2();

            }
            this->currentCoord = bucket->start + bucket_index;
            this->currentValPtr = &bucket->data[bucket_index];
        }

        virtual pair<coord_t, DTYPE&>* getCurrentPairPtr() {
            if (this->at_end)
                throw range_error("Iterator reached end and was queried for value!");
            current_pair.~pair<coord_t, DTYPE&>();
			new(&current_pair) pair<coord_t, DTYPE&>(bucket->start + bucket_index, bucket->data[bucket_index]);
            return &current_pair;
        }

        virtual void increment() {
            if (this->top_only) {
                 do {
                     increment2();
                 } while (!this->at_end && !bucket->top[bucket_index]);
             } else
                 increment2();
            this->currentCoord = bucket->start + bucket_index;
            this->currentValPtr = &bucket->data[bucket_index];
        }

        SparseIterator* clone() {
            SparseIterator* result = new SparseIterator(field, this->top_only, only_level);
            result->bucket_index = bucket_index;
            result->map_iter = map_iter;
            result->bucket = bucket;
            return result;
        }

    private:
        void increment2() {
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
                     this->at_end = true;
                     return;
                 }
             }
         }

         bool top_only;
         map_iter_t map_iter;
         size_t bucket_index;
         int only_level;
         Bucket *bucket;
         DTYPE intermediate; // first ref in current_pair
         pair<coord_t, DTYPE&> current_pair;
         SparseField<DTYPE, HCSTYPE> *field;
    };

    // Iterator methods & class
    //iterator dummy = iterator(this);
    typename Field<DTYPE, HCSTYPE>::Iterator begin(bool top_only = false, int only_level = -1) {
        return typename Field<DTYPE, HCSTYPE>::Iterator(new SparseIterator(this, top_only, only_level));
    }

    typename Field<DTYPE, HCSTYPE>::Iterator end() {    // Just dummy, the begin iterator determines termination
        return NULL;//typename Field<DTYPE, HCSTYPE>::Iterator(new SparseIterator(NULL));//new SparseIterator(this, top_only, only_level));
    }

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

    // Do we have a value for this coord? And if yes, make sure it is in _current
    // A bucket's end coord is its last existing coord
    bool exists(coord_t coord) {
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
    DTYPE get(coord_t coord, bool use_non_top = true) {
        DTYPE result = 0;
        Field<DTYPE, HCSTYPE>::get(coord, result, use_non_top);
        return result;
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
            if (highest == 2) // 0-Bucket has only one coord that should have been filled.
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



    // .. and all levels below.
    // This routine DELETES everything in the field and is meant as an initializer.
    // If there are elements present, it throws.
    void createEntireLevel(level_t level) {
        if (data.size() > 1)
            throw range_error("Not empty!");
        data[1]->setTop(1, false);
        for (level_t l = 2; l <= level; l++) {
            coord_t level_start = hcs.CreateMinLevel(l);
            coord_t level_end = hcs.CreateMaxLevel(l);
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
            return; // ? nothing to do...
        Bucket* bucket = new Bucket(lower_corner, upper_corner);
        data[lower_corner] = bucket;
        fill(bucket->top.begin(), bucket->top.end(), true); // Mark as top
        fill(bucket->data.begin(), bucket->data.end(), coord_bucket->get(coord));   // Set values from orig coord
        // Now the original coord is not top anymore...
        coord_bucket->setTop(coord, false);
        _current = bucket;  // grant immediate access to new coords
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

    // Read-write access to existing coords. For (probably) non-existing, use get() and retrieve interpolated values.
    // the value set to bracket_behavior applies.
    DTYPE& operator[](coord_t coord) {
        if (!this->exists(coord)) {
            switch (this->bracket_behavior) {
            case Field<DTYPE, HCSTYPE>::BR_THROW:
                throw range_error("[]: Coord does not exist");
            case Field<DTYPE, HCSTYPE>::BR_INTERP:
                intermediate = get(coord);
                return intermediate;
            case Field<DTYPE, HCSTYPE>::BR_REFINE:
                refineTo(coord);
                return this->_current->get(coord);
            case Field<DTYPE, HCSTYPE>::BR_NOTHING:
                return intermediate;
            }
        }
        return this->_current->get(coord);
    }



    // Assignment operator requires equal structure, dirty-check with data.size()
    // isTop is not copied because of assumption of equal structure
    SparseField &operator=(const SparseField& f){
        //cout << "XCOPY\n";
        assert(("= Operator would alter structure. if this is intended, call takeStructure(x) first!",
                data.size() == f.data.size()));
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

    SparseField<DTYPE, HCSTYPE> operator-() const { SparseField<DTYPE, HCSTYPE> result = *this; for (auto e : result) e.second = -e.second; return result;}

    // Clears the field and takes the same coordinate structure as the provided field, without copying their
    // values. The provided field may have a different DTYPE. The newly created coords are initialized with zero.
    template <typename DTYPE2>
    void takeStructure(SparseField<DTYPE2, HCSTYPE> &f) {
        if (sameStructure(f))
            return;
        clear();
        for (auto e : f.data) {
            auto *b = e.second;
            Bucket *bn = new Bucket(b->start, b->end);
            bn->top = b->top;
            for (coord_t c = bn->start; c <= bn->end; c++)
                bn->get(c) = 0;
            data[b->start] = bn;
        }
    }

    template <typename DTYPE2>
    void takeStructure(Field<DTYPE2, HCSTYPE> &f) {
    	// could be slow...
    	auto old_bracket_behavior = this->bracket_behavior;
    	this->bracket_behavior = Field<DTYPE,HCSTYPE>::BR_REFINE;
    	for (auto it = f.begin(true); it != f.end(); ++it) {
    		refineTo(it->first);
    	}
    }

    // Tests if the provided field has the same structure.
    // The provided field may have a different DTYPE. The newly created coords are initialized with zero.
    template <typename DTYPE2>
    bool sameStructure(SparseField<DTYPE2, HCSTYPE> &f) {
        if (f.data.size() != data.size())
            return false;

        auto it_self = data.begin();
        auto it_foreign = f.data.begin();

        while (it_self != data.end() || it_foreign != f.data.end()) {
            if (it_self->first != it_foreign->first)
                return false;
            if ((it_self->second)->end != (it_foreign->second)->end)
                return false;
            ++it_self; ++it_foreign;
        }
        return true;
    }

    // Empties all data
    void clear() {
        for (auto e : data)
            delete e.second;
        data.clear();
        data[1] = new Bucket(1, 1);
        data[1]->setTop(1, true);
        _current = NULL;
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

        if (b_start == start && b_end == end) { // Bucket matches range
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

        } else {    // coords are within a bucket, need to split...
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
        friend class SparseField; // Field can access private

        coord_t         start, end;
        vector<DTYPE>   data;
        vector<char>    top;

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


