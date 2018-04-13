#include "includes.hpp"
#include "solver.hpp"

// Dimension setup

template<int dimensions>
data_t poisson_err(Field<data_t, HCS<dimensions> > &err, Field<data_t, HCS<dimensions> > &x, Field<data_t, HCS<dimensions> > &b, int level_only) {
	//err = x;
	//err = 0;
	data_t norm = 0;
	for (auto e = x.begin(false, level_only); e != x.end(); ++e) {
		coord_t c = e->first;
		data_t &value = e->second;
		//data_t dist = x.hcs.getDistance(c, 0);
		//data_t vol = pow(dist, dimensions); // volume of cell
		//data_t area = pow(dist, dimensions - 1); // area of cell-wall
		data_t sum = 0;
		data_t factor = x.hcs.parts;
		for (int ne = 0; ne < x.hcs.parts; ne++) {
			coord_t ne_coord = x.hcs.getNeighbor(c, ne);
			data_t ne_value = x.get(ne_coord, true);
			if (x.hcs.IsBoundary(ne_coord)) {
			    sum += 2 * ne_value;
			    factor++;
			} else {
                sum += ne_value;
			}
		}
		sum /= factor;
		err[c] = value + (b[c] - sum);// * vol * 0.5;
		if (x.hcs.GetLevel(c) < 4)
		    cout << x.hcs.toString(c) << " X " << x[c] << " ERR " << err[c] << endl;
		norm += err[c] * err[c];
	}
	//cout << "Residual norm: " << norm << endl;
	return norm;
}
// 1-0.25 + 3 * -0.25

template<class DTYPE, int dimensions>
void poisson(Field<DTYPE, HCS<dimensions> > &x, Field<DTYPE, HCS<dimensions> > &b) {
    HCS<dimensions> hcs = x.hcs;


    DenseField<DTYPE, HCS<dimensions> > err_b(b.getHighestLevel());
    DenseField<DTYPE, HCS<dimensions> > err_return(b.getHighestLevel());

    // Downpropagate b
    // Zero all non-top
    for (auto element : b)
        if (!b.isTop(element.first))
            element.second = 0;

    for (level_t level = b.getHighestLevel(); level > 0; level--) {
        for (auto it = b.begin(false, level); it != x.end(); ++it) {
            const coord_t &coord = (*it).first;
            DTYPE &value = (*it).second;
            auto coeffs = hcs.getCoeffs(coord);
            //b.correct_neumann(&coeffs[0]);
            //value = 0;
            for (auto coeff : coeffs) {
                coord_t xco = coeff.first;
                if (hcs.IsBoundary(xco))
                //    continue;
                   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
                b[xco] += coeff.second * value;
                //value += coeff.second * x.get(xco);
            }
        }
    }

    err_return = 0;
    err_b = 0;
    // Up
    for (level_t level = 0; level <= x.getHighestLevel(); level++) {
        data_t level_norm = 1;
        int iter = 0;
        while (level_norm > 1e-10) {
        	// Up propagation of x
			for (auto it = x.begin(false, level); it != x.end(); ++it) {
				const coord_t &coord = (*it).first;
				DTYPE &value = (*it).second;
				auto coeffs = hcs.getCoeffs(coord);
				//x.correct_neumann(&coeffs[0]);
				value = 0;
				for (auto coeff : coeffs) {
					coord_t xco = coeff.first;
					if (hcs.IsBoundary(xco))
					   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
					value += coeff.second * x.get(xco) + coeff.second * b.get(xco) + coeff.second * err_b.get(xco);
				}
			}
	        // Calc error on level
	        //level_norm = poisson_err<dimensions>(err_return, x, b, level);
	        //cout << "Level " << level << " Iteration " << iter++ << " Norm " << level_norm << endl;

	        // Correct error
			for (auto it = x.begin(false, level); it != x.end(); ++it) {
				const coord_t &coord = (*it).first;
				(*it).second -= err_return[coord];
			}
	        level_norm = poisson_err<dimensions>(err_return, x, b, level);

	        // Propagate error down one level
	        for (auto it = err_return.begin(false, level); it != x.end(); ++it) {
	            const coord_t &coord = (*it).first;
	            DTYPE &value = (*it).second;
	            auto coeffs = hcs.getCoeffs(coord);
	            //b.correct_neumann(&coeffs[0]);
	            //value = 0;
	            for (auto coeff : coeffs) {
	                coord_t xco = coeff.first;
	                if (hcs.IsBoundary(xco))
	                    continue;
	                //   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
	                err_b[xco] -= coeff.second * value;
	                //value += coeff.second * x.get(xco);
	            }
	        }
	        if (iter > 30)
	        	exit(1);
        }
    }
    cout << "CENTER: " << x[1] << endl;
}

template<class DTYPE, int dimensions>
void laplace(Field<DTYPE, HCS<dimensions> > &x) {
    HCS<dimensions> hcs = x.hcs;
    DenseField<DTYPE, HCS<dimensions> > b(x.getHighestLevel());
    b = 0.;

/*
    x[1]=0.25;
    for (level_t level = 1; level <= x.getHighestLevel(); level++) {
        for op(auto it = x.begin(false, level); it != x.end(); ++it) {
            const coord_t &coord = (*it).first;
            DTYPE &value = (*it).second;
            auto coeffs = hcs.getCoeffs2(coord);
            //x.correct_neumann(&coeffs[0]);
            value = 0;
            for (auto coeff : coeffs) {
                coord_t xco = coeff.first;
                //if (hcs.IsBoundary(xco))
                //   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
                value += coeff.second * x.get(xco);
            }
        }
    }
    cout << "CENTER: " << x[1] << endl;
    */
    poisson<DTYPE, dimensions>(x, b);

}


int main(int argc, char **argv) {

	// Resolution setup
	const int max_level = 9;	// Max level for "C" field

	DenseScalarField2 x(max_level), b(max_level), err(max_level);

    auto h = x.hcs;
	b.boundary[0] = x.boundary[0] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
	b.boundary[1] = x.boundary[1] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
	b.boundary[2] = x.boundary[2] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
	b.boundary[3] = x.boundary[3] = [](ScalarField2 *self, coord_t cc)->data_t { return 1;};

	coord_t c = h.createFromList({0,0});
	auto cc = h.getCoeffs(c);
	for (auto e : cc) {
		cout << h.toString(e.first) << " W: " << e.second << endl;
	}
//exit(1);
    // Laplace test
	laplace<data_t, 2>(x);

    // Poisson test

	b = 0;

    //b[h.createFromUnscaled(9, {256, 256})] = (512 *512);
    //poisson<data_t, 2>(x, b);



	poisson_err<2>(err, x, b, max_level);
    write_raw2("test9.raw", x);

	write_pgm("test9.pgm", x, max_level);

    write_pgm("test9_err.pgm", err, max_level);
    write_raw2("test9_err.raw", err);


}
