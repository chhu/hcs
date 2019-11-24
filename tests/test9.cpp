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

/*
data_t greens(data_t dist, int dim) {
	if (dist == 0 || dim <= 1)
		return 1;
	if (dim == 2)
		return (-1./ (2. * M_PI)) * log(dist);
	data_t unity_ball = pow(M_PI, data_t(dim) / 2) / tgamma(data_t(dim) / 2. + 1.);
	return (1. / (dim * (dim - 2) * unity_ball)) * 1. / pow(dist, dim - 2);
}

template <int dim>
data_t greens(coord_t source, coord_t dest) {
	HCS<dim> hcs;
	data_t dist = hcs.getDistance(source, dest);
	return greens(dist, dim);
}

template <int dim>
data_t greens(Tensor1<data_t, dim> &source,  Tensor1<data_t, dim> &dest) {
	data_t dist = (dest - source).length();
	return greens(dist, dim);
}
*/

// generic test version
template <unsigned dim>
void poisson(DenseScalarField<dim> &x, DenseScalarField<dim> &b) {
    typedef Tensor1<data_t, dim> Vec;

	HCS<dim> &hcs = x.hcs;

//    DenseVectorField<dim> source_vector(b.getHighestLevel());
//    DenseScalarField<dim> upscale(b.getHighestLevel());

    x.takeStructure(b);
    x = 0;

    // Downpropagate
    // Zero all non-top
    for (auto element : x)
        if (!x.isTop(element.first))
            element.second = 0;
    //b.pack();


//    for (level_t level = x.getHighestLevel(); level > 0; level--) {
//        for (auto it = x.begin(false, level); it != x.end(); ++it) {
//            const coord_t &coord = (*it).first;
//            data_t value = (*it).second;
//            //data_t div = hcs.getNeighborCount(coord, false);
//            auto coeffs = hcs.getCoeffs(coord);
//            for (auto coeff : coeffs) {
//                coord_t xco = coeff.first;
//                if (hcs.IsBoundary(xco)) {
//                	//mul-=1;
//                    //   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
//                    continue;
//                }
//                source_vector[xco] += c2v<dim>(hcs, coord, xco);
//                x[xco] += value * coeff.second;// / div;// * coeff.second * 1;
//
//             }
//        }
//    }
////return;
//    //    err_return = 0;
//    //    err_b = 0;
    // Up
    for (level_t level = 0; level <= x.getHighestLevel(); level++) {
		// Up propagation of x
		for (auto it = x.begin(false, level); it != x.end(); ++it) {
			const coord_t &coord = (*it).first;
			data_t &value = (*it).second;

			auto coeffs = hcs.getCoeffs(coord);

			data_t v = b[coord];
			for (auto coeff : coeffs) {
				coord_t xco = coeff.first;
			//	if (hcs.IsBoundary(xco)) {
			//		continue;
			//	}
				//vec += coeff.second * source_vector.get(xco);
				v += coeff.second * x.get(xco);
			}
			x[coord] = v;// / coeffs.size();
		}
	}
}

template<class DTYPE, int dimensions>
void laplace(Field<DTYPE, HCS<dimensions> > &x) {
    HCS<dimensions> hcs = x.hcs;
    DenseField<DTYPE, HCS<dimensions> > b(x.getHighestLevel());
    b = 0.;
    //poisson<DTYPE, dimensions>(x, b);
}


int main(int argc, char **argv) {

	// Resolution setup
	const int max_level = 9;	// Max level for "C" field

	DenseScalarField<2> x(max_level), b(max_level), err(max_level);

    auto h = x.hcs;
/*
    b.boundary[0] = x.boundary[0] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
	b.boundary[1] = x.boundary[1] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
	b.boundary[2] = x.boundary[2] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
	b.boundary[3] = x.boundary[3] = [](ScalarField2 *self, coord_t cc)->data_t { return 0;};
/*
	coord_t c = h.createFromList({0,0});
	auto cc = h.getCoeffs(c);
	for (auto e : cc) {
		cout << h.toString(e.first) << " W: " << e.second << endl;
	}
*/
	// Laplace test
	//laplace<data_t, 2>(x);

    // Poisson test

	//b = 1.;///512;//pow(1./512,2);
    b[h.createFromUnscaled(9, {255, 255})] = 512*512;
    //b[h.createFromUnscaled(9, {128})] = 1;
    //b[h.createFromUnscaled(1, {0,0})] = 1;
//    b[h.createFromUnscaled(1, {0,1})] = 1;
//    b[h.createFromUnscaled(1, {1,0})] = -1;
    //b[h.createFromUnscaled(1, {1,1})] = -1;
//    b[h.createFromPosition(3,{0.5,0.5})] = 1;
    //b[h.createFromPosition(2,{0.5,0.5})] = 1;

    b.propagate();
    //b.pack();
    //b[1]=1;

    /*
	b[h.createFromUnscaled(9, {256, 256})] = 1;//(512 *512);
    b[h.createFromUnscaled(9, {255, 255})] = 1;//(512 *512);
    b[h.createFromUnscaled(9, {256, 255})] = 1;//(512 *512);
    b[h.createFromUnscaled(9, {255, 256})] = 1;//(512 *512);
    */
    poisson<2>(x, b);
    b = x;
    //b.pack();
   // poisson<2>(x, b);


//	poisson_err<2>(err, x, b, max_level);
	x.write("test9_9.raw");
	x.write("test9_8.raw", 8);
	x.write("test9_7.raw", 7);
	x.write("test9_6.raw", 6);
	x.write("test9_5.raw", 5);
	x.write("test9_4.raw", 4);
	x.write("test9_3.raw", 3);
	x.write("test9_2.raw", 2);
	x.write("test9_1.raw", 1);
	b.write("btest9_9.raw");
	b.write("btest9_8.raw", 8);
	b.write("btest9_7.raw", 7);
	b.write("btest9_6.raw", 6);
	b.write("btest9_5.raw", 5);
	b.write("btest9_4.raw", 4);
	b.write("btest9_3.raw", 3);
	b.write("btest9_2.raw", 2);
	b.write("btest9_1.raw", 1);
	//write_pgm("test9.pgm", x, max_level);
	cout.precision(17);
	cout << "B[1] = " << b[1] << " X[1] = " << x[1] << endl;

	//DenseVectorField2 grad_x;
	//grad<2>(x, grad_x);
	//div<2>(grad_x, x);
	//x.write("test9_err.raw");



}
