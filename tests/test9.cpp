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

// non-generic test version for 2D only
void poisson2(DenseScalarField1 &x, DenseScalarField1 &b) {
    auto hcs = x.hcs;

    DenseVectorField2 grad_x(b.getHighestLevel());

    x.takeStructure(b);
    x = b;

    // Downpropagate, establish gradient
    // Zero all non-top
    for (auto element : x)
        if (!x.isTop(element.first))
            element.second = 0;
    for (level_t level = x.getHighestLevel(); level > 0; level--) {
        for (auto it = x.begin(false, level); it != x.end(); ++it) {
            const coord_t &coord = (*it).first;
            data_t dist = 1. / pow(2, level);
            data_t volume = dist * dist;
            data_t value = (*it).second;// * volume;
            //Vec2 grad = grad_x[coord];

            auto coeffs = hcs.getCoeffs(coord);
            data_t mul = 2;
            for (auto coeff : coeffs) {
                coord_t xco = coeff.first;

                if (hcs.IsBoundary(xco)) {
                	mul-=1;
                    continue;
                }
                //   xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
                x[xco] += coeff.second * value * 2;// * (1./2) * (level*2);//0.25; // 4 times the volume from source
             }
        }
    }
//return;
    //    err_return = 0;
    //    err_b = 0;
    // Up
    for (level_t level = 1; level <= x.getHighestLevel(); level++) {
		// Up propagation of x
		for (auto it = x.begin(false, level); it != x.end(); ++it) {
            data_t dist = 1. / pow(2, level);
            data_t volume = dist * dist;
			const coord_t &coord = (*it).first;
			data_t &value = (*it).second;
			auto coeffs = hcs.getCoeffs(coord);
			data_t v = 0;
			for (auto coeff : coeffs) {
				coord_t xco = coeff.first;
				if (hcs.IsBoundary(xco)) {
					continue;
					coord_t origin = hcs.removeBoundary(xco);
					v += coeff.second * -x.get(origin) * 0.5;
					//xco = (xco & ~(coord_t)0xFFFFFFFF) | coord;
				} else
					v += coeff.second  * x.get(xco);// / volume * 0.25;
			}
			value += v;
			//value *= 0.5;
		}
	}
//}
    cout << "CENTER: " << x[1] << endl;
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

	DenseScalarField1 x(max_level), b(max_level), err(max_level);

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
    b[h.createFromUnscaled(9, {255})] = 1;
    b[h.createFromUnscaled(9, {128})] = 1;
    /*
	b[h.createFromUnscaled(9, {256, 256})] = 1;//(512 *512);
    b[h.createFromUnscaled(9, {255, 255})] = 1;//(512 *512);
    b[h.createFromUnscaled(9, {256, 255})] = 1;//(512 *512);
    b[h.createFromUnscaled(9, {255, 256})] = 1;//(512 *512);
    */
    poisson2(x, b);

	//poisson_err<2>(err, x, b, max_level);
	x.write("test9_9.raw");
	x.write("test9_8.raw", 8);
	x.write("test9_7.raw", 7);
	x.write("test9_6.raw", 6);
	x.write("test9_5.raw", 5);
	x.write("test9_4.raw", 4);
	x.write("test9_3.raw", 3);
	x.write("test9_2.raw", 2);
	x.write("test9_1.raw", 1);
	//write_pgm("test9.pgm", x, max_level);

	//DenseVectorField2 grad_x;
	//grad<2>(x, grad_x);
	//div<2>(grad_x, x);
	//x.write("test9_err.raw");



}
