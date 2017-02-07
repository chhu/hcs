#include "includes.hpp"

template<typename DTYPE, typename FTYPE>
class Matrix {
public:
	function<DTYPE(coord_t, DTYPE, FTYPE&)> mul_stencil;
	void mul(FTYPE& x, FTYPE& result) {
		auto dual_it = x.begin(&result);
		while (dual_it != x.end()) {  // one row per while-loop
			coord_t coord = get<0>(*dual_it);
			get<2>(*dual_it) = mul_stencil(coord, get<1>(*dual_it), x);
			++dual_it;
		}
	}
};

template<typename DTYPE, typename FTYPE>
class Solver {
public:
	data_t dot(FTYPE& a, FTYPE& b) {
		data_t result = 0;
		for (auto iter = a.begin(&b); iter != a.end(); ++iter)
			result += std::get<1>(*iter) * std::get<2>(*iter);
		return result;
	}

	data_t norm(FTYPE& a) {
		data_t result = 0;
		for (auto e : a)
			result += e.second * e.second;
		return result;
	}

	int solve(Matrix<DTYPE,FTYPE>& M, FTYPE& x, FTYPE& b, int max_it, data_t r_tol, data_t a_tol) {

		typedef double T;
		long double rho_1, rho_2, alpha, beta, omega, norm_b, norm_r, start_t, end_t;
		//int n = x.size();
		//int n_sim = this->sim->grid->gross_points;
		int iter = 0;
		int debug = 2;
		int is_nan;

		//valarray<T> p(n), phat(n), s(n), shat(n), t(n), v(n), r(n), rtilde(n), Ap(n);
		T Ap = 1.;//1./(-4.);
		FTYPE p = x;
		p = 0.;
		FTYPE phat = p;
		FTYPE s = p;
		FTYPE shat = p;
		FTYPE t = p;
		FTYPE v = p;
		FTYPE r = p;
		FTYPE res = p;
		FTYPE rtilde = p;

		is_nan = 0;
		iter = 0;
		phat = 0;
		shat = 0;

		// r = b - A*x
		M.mul(x, r);
		r *= -1;
		r += b;

		norm_b = this->norm(b);
		norm_r = this->norm(r);

		if (isnan(norm_r) || isnan(norm_b)) {
			cout << "BICGS found nan solution (iter 0)" << norm_b << " " << norm_r << endl;
			norm_r = 1000.;
			exit(4);
		}

		if (norm_r < a_tol) return iter;

		rtilde = r;

		if (debug >= 2)
			printf("BiCGStab Solver started: B-Norm = %g; R-Norm = %g\n", T(norm_b), T(norm_r));

		while ((((norm_r > a_tol) && (norm_r / (norm_b == 0 ? norm_r : norm_b) > r_tol) && (iter < max_it)) || (iter < 1))) {

			rho_1 = this->dot(rtilde, r);
			if (rho_1 == T(0.)) {
				printf("bicg breakdown: r_tilde * r = 0; iter: %d norm_r: %g norm_b: %g \n", iter, norm_r, norm_b);
				exit(4);
			}

			if (iter == 0)
				p = r;
			else {
				beta = (rho_1 / rho_2) * (alpha / omega);
				p += v * T(-omega);
				p *= beta;
				p += r;
			}

			phat = p * Ap;

			M.mul(phat, v);

			alpha = rho_1 / this->dot(v, rtilde);

			s = r - (T(alpha) * v);
			shat = s * Ap;
			M.mul(shat, t);

			omega = this->dot(t, s) / this->norm(t);
			if (omega == T(0.)) {
				cerr << "bicg breakdown: omega = 0\n";
				exit(4);
			}

			x += T(alpha) * phat;
			x += T(omega) * shat;

			r = s;
			r -= T(omega) * t;

			rho_2 = rho_1;
			norm_r = this->norm(r);

			if (isnan(norm_r)) {
	           	cout << "BICGS found nan r-norm at iteration " << iter << endl;
	           	exit(8);
			}

			if (debug > 2) {
				printf(" iter  %d: res = %g\n", iter, T(norm_r / norm_b));
			} else if (debug > 1) {
				if ((iter % 100) == 0) {
					printf("Iteration %d: res = %g t = %g sec\n", iter, T(norm_r / norm_b), T(end_t - start_t));
				}
			}
			iter++;
		}
		if (iter >= max_it)
			printf("BiCGS: Maximum number of iteration reached! Best rtol: %g\n", T(norm_r / norm_b));

		return iter;
	}

	// Steepest Descent solver
	int solve2(Matrix<DTYPE,FTYPE>& M, FTYPE& x, FTYPE& b, int max_it, data_t r_tol, data_t a_tol) {
		int iter = 0;
		int debug = 3;

		typedef double T;
		T rho_1, rho_2, alpha, beta, omega, norm_b, norm_r, start_t, end_t;
		T Ap = 1.;//1./(-4.);
		FTYPE r = x;
		FTYPE res = x;

		iter = 0;

		// r = b - A*x
		M.mul(x, r);
		r *= -1;
		r += b;

		norm_b = norm(b);
		norm_r = norm(r);
		if (norm_b == 0)
			norm_b = norm_r;

		if (norm_r < a_tol) return iter;

		if (debug >= 2)
			printf("SD Solver started: B-Norm = %g; R-Norm = %g\n", T(norm_b), T(norm_r));
		while ((((norm_r > a_tol) && (norm_r / (norm_b == 0 ? norm_r : norm_b) > r_tol) && (iter < max_it)) || (iter < 1))) {

			norm_r = this->norm(r);
			M.mul(r, res); //res = Ar
			alpha = norm_r / this->dot(r, res); // alpha=r^Tr/r^TAr
			alpha *= 0.5;
			x += T(alpha) * r;
			r -= T(alpha) * res;
			iter++;
			printf(" iter  %d: res = %g\n", iter, T(norm_r / norm_b));
		}
		if (iter >= max_it)
			printf("Maximum number of iteration reached! Final: %g\n",T(norm_r / norm_b));

		// check for failed convergence...
		return iter;
	}
};

int main(int argc, char **argv) {

	#ifdef __BMI2__
	printf("Compiled with BMI2!\n");
	#else
	printf("NO BMI2\n");
	#endif

	H1 h1;
	H2 h2;
	H3 h3;
// 0.12475
	ScalarField1 x;
	x.createEntireLevel(8);
	cout << "Solver Test\n";
	//x[h1.createFromList({0,1,0,0,0})] = 1;
	//x[h1.createFromList({0})] = 1;
	//x[0] = 0.5;
	//x.propagate();
	//x.upAverage();
	Solver<data_t, ScalarField1> solver;
	cout << solver.norm(x) << " " << x.nElements() << " " << x.nElementsTop() << endl;
	Matrix<data_t, ScalarField1> M;
/*
 * bc_correct = 2 on boundary
	value = (0.5 * (var_k[i] + var_k[n_]) * f.area * bc_correct) / (f.distance * vol);
	A.add_element(n_,-value);
	A.add_element(i,  value);
*/
	M.mul_stencil = [](coord_t coord, data_t x_val, ScalarField1& x)->data_t {
		//if (!x.isTop(coord))
		//	return 0;
		coord_t level = (coord_t)1 << x.hcs.GetLevel(coord);
		data_t dist = 1. / level;  // Neighbor distance, assuming all directions with equal scale and scale = 1
		data_t vol = dist;
		data_t row_result = 0;//-2. * x_val; // This is the main diagonal entry
		for (int ne_idx = 0; ne_idx < x.hcs.parts; ne_idx++) {
			//Vec2 normal(x.hcs.getNeighborDirection(ne_idx));
			coord_t ne_coord = x.hcs.getNeighbor(coord, ne_idx);
			data_t ne_val = x.get(ne_coord);
			data_t coeff  = (1. ) / (dist * vol);

			if (x.hcs.IsBoundary(ne_coord)) {
				coeff *= 2;
				//ne_val = origin.get(ne_coord);
			}
			//data_t coeff  = (0.5 * (var_k[coord] + var_k[ne_coord]) * dist) / (dist * vol);
			row_result += coeff * x_val - coeff * ne_val;
		}
		return row_result;
	};

	x = 0;
	ScalarField1 b = x;
	b = 1.;
	b[h1.createFromPosition(8, {0.88})] = 1000;
	b.propagate();
	//x.boundary[0] = [](coord_t c)->data_t {
	//	return 1;
	//};
	cout << "ITER: " << solver.solve(M, x, b, 22000, 1e-7, 1e-12) << "\n";
	data_t max_val = 0;
	ofstream out("t6data.txt");
	for (auto e : x) {
		//if (!x.isTop(e.first))
		//	continue;
		//if (h1.GetLevel(e.first) != 1)
		//	continue;
		out << h1.getPosition(e.first)[0] << " " << e.second << endl;
		if (e.second > max_val)
			max_val = e.second;
	}
	out.close();
	cout.precision(17);
	cout << " Max value: " << max_val << endl;
	//M.mul(x, b);
	//write_pgm("test6.pgm", x, 8);

}
