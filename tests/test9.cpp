#include "includes.hpp"
#include "solver.hpp"

// Dimension setup

template<int dimensions>
void poisson_err(Field<data_t, HCS<dimensions> > &err, Field<data_t, HCS<dimensions> > &x, Field<data_t, HCS<dimensions> > &b) {
	err = x;
	err = 0;
	data_t norm = 0;
	for (auto e = x.begin(false); e != x.end(); ++e) {
		coord_t c = e->first;
		data_t &value = e->second;
		data_t dist = x.hcs.getDistance(c, 0);
		data_t vol = pow(dist, dimensions); // volume of cell
		data_t area = pow(dist, dimensions - 1); // area of cell-wall
		data_t sum = 0;
		for (int ne = 0; ne < x.hcs.parts; ne++) {
			coord_t ne_coord = x.hcs.getNeighbor(c, ne);
			data_t ne_value = x.get(ne_coord, true);
			if (x.hcs.IsBoundary(ne_coord)) {
			    sum += ((ne_value  - value) * 2 * area) / (dist * vol);
			} else {
                sum += ((ne_value  - value) * area) / (dist * vol);
			}
		}
		err[c] = abs(b[c] - sum);
		if (x.hcs.GetLevel(c) < 4)
		    cout << x.hcs.toString(c) << " X " << x[c] << " ERR " << err[c] << endl;
		norm += sum * sum;
	}
	cout << "Residual norm: " << norm << endl;
	//err.propagate();
}
// 1-0.25 + 3 * -0.25

double CINT(double p[4], double x) {
    return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double bicubic(double p[4][4], double x, double y) {
    double arr[4];
    arr[0] = CINT(p[0], y);
    arr[1] = CINT(p[1], y);
    arr[2] = CINT(p[2], y);
    arr[3] = CINT(p[3], y);
    return CINT(arr, x);
}

//template<class DTYPE, int dimensions>
//void laplace(Field<DTYPE, HCS<dimensions> > &x) {
//    HCS<dimensions> hcs = x.hcs;
//
//    for (level_t level = 0; level <= x.getHighestLevel(); level++) {
//        for (auto it = x.begin(false, level); it != x.end(); ++it) {
//            const coord_t &coord = (*it).first;
//            DTYPE &value = (*it).second;
//            uint32_t bc_set=0;
//            auto bbox = hcs.getCoeffCoords(coord, bc_set);
//            coord_t bcc[4][4];
//            bcc[1][1] = bbox[0];
//            bcc[2][1] = bbox[1];
//            bcc[1][2] = bbox[2];
//            bcc[2][2] = bbox[3];
//
//
//            //value = 0;
//            /*for (auto coeff : coeffs) {
//                coord_t xco = coeff.first;
//                if (hcs.IsBoundary(xco))
//                   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
//                value += coeff.second * x.get(xco);
//            }*/
//        }
//    }
//    cout << "CENTER: " << x[1] << endl;
//}

template<class DTYPE, int dimensions>
void laplace(Field<DTYPE, HCS<dimensions> > &x) {
    HCS<dimensions> hcs = x.hcs;

    for (level_t level = 0; level <= x.getHighestLevel(); level++) {
        for (auto it = x.begin(false, level); it != x.end(); ++it) {
            const coord_t &coord = (*it).first;
            DTYPE &value = (*it).second;
            auto coeffs = hcs.getCoeffs(coord);
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
}

template<class DTYPE, int dimensions>
void poisson(Field<DTYPE, HCS<dimensions> > &x, Field<DTYPE, HCS<dimensions> > &b) {
    HCS<dimensions> hcs = x.hcs;

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
            b.correct_neumann(&coeffs[0]);
            //value = 0;
            for (auto coeff : coeffs) {
                coord_t xco = coeff.first;
                if (hcs.IsBoundary(xco))
                    continue;
                //   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
                b[xco] += coeff.second * value;
                //value += coeff.second * x.get(xco);
            }
        }
    }


    // Up
    for (level_t level = 0; level <= x.getHighestLevel(); level++) {
        for (auto it = x.begin(false, level); it != x.end(); ++it) {
            const coord_t &coord = (*it).first;
            DTYPE &value = (*it).second;
            auto coeffs = hcs.getCoeffs(coord);
            x.correct_neumann(&coeffs[0]);
            value = 0;
            for (auto coeff : coeffs) {
                coord_t xco = coeff.first;
                if (hcs.IsBoundary(xco))
                   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
                value += coeff.second * x.get(xco) + coeff.second * b.get(xco);
            }
        }
    }
    cout << "CENTER: " << x[1] << endl;
}


int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif



	// Resolution setup
	const int max_level = 9;	// Max level for "C" field



	DenseScalarField2 x(max_level), b(max_level), err(max_level);

    auto h = x.hcs;
/*
    uint32_t bcs = 0;
    coord_t cc = h.createFromList({1,2});
    cout << h.toString(cc) << endl<<endl;
    auto ac = h.getCoeffCoords(cc, bcs);
    for (auto e : ac)
        cout << h.toString(e) << endl;
return 0;
*/
	b.boundary[0] = x.boundary[0] = [](ScalarField2 *self, coord_t cc)->data_t {
	    return 0;
	    coord_t interior = self->hcs.removeBoundary(cc);
        double y = self->hcs.getPosition(interior)[1];
        //return y;
        //return y > 0.5;
        return sin(y * M_PI);
    };
	b.boundary[1] = x.boundary[1] = [](ScalarField2 *self, coord_t cc)->data_t {
        return 0;
        coord_t interior = self->hcs.removeBoundary(cc);
        double y = self->hcs.getPosition(interior)[1];
        return y;
        return cos(y * 2 * M_PI);
    };
	b.boundary[2] = x.boundary[2] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
	b.boundary[3] = x.boundary[3] = [](ScalarField2 *self, coord_t cc)->data_t { return 1;};

    coord_t test = h.createFromList({0});
    cout << "TARGET " << h.toString(test) << endl;
    auto coeffs = h.getCoeffs(test);
    //x.correct_neumann(&coeffs[0]);
    for (auto coeff : coeffs) {
        cout << h.toString(coeff.first) << " W: " << coeff.second << endl;
    };
    test = h.createFromList({0,0});
    cout << "TARGET " << h.toString(test) << endl;
    coeffs = h.getCoeffs(test);
    //x.correct_neumann(&coeffs[0]);
    for (auto coeff : coeffs) {
        cout << h.toString(coeff.first) << " W: " << coeff.second << endl;
    }
//return 0;

    // Laplace test
	laplace<data_t, 2>(x);

    // Poisson test
    b = 0;
    //b[h.createFromUnscaled(9, {256, 256})] = (512 *512);
    //poisson<data_t, 2>(x, b);



	poisson_err<2>(err, x, b);
    write_raw2("test9.raw", x);

	write_pgm("test9.pgm", x, max_level);

    write_pgm("test9_err.pgm", err, max_level);
    write_raw2("test9_err.raw", err);

    //cout << x[h.createFromPosition(9,{0.75,0.25})] << endl;

    //read_raw2("test9_sol.raw", x, 9);
    //cout << x[h.createFromPosition(9,{0.75,0.25})] << endl;
    //x.propagate();
    //poisson_err<2>(err, x, b);
    //write_raw2("test9_sol_err.raw", err);

	/*
	write_txt("test9e.txt", err,false);
	write_txt("test9.txt", x);
	write_txt("test9b.txt", b,false);

*/

}
