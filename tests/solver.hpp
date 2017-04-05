/*
 * solver.hpp
 *
 *	A simple Matrix-free "matrix" and sparse-solver for the HCS tests
 *
 *  Created on: Apr 5, 2017
 *      Author: huettig
 */

#ifndef TESTS_SOLVER_HPP_
#define TESTS_SOLVER_HPP_

using namespace std;

template<typename DTYPE, typename FTYPE>
class Matrix {
public:
	function<DTYPE(coord_t, DTYPE, FTYPE&)> mul_stencil;	// (coord, value of coord, field) => mul result
	void mul(FTYPE& x, FTYPE& result) {
		for (auto dual_it = x.begin(&result, true); dual_it != x.end(); ++dual_it)
			get<2>(*dual_it) = mul_stencil(get<0>(*dual_it), get<1>(*dual_it), x);
	}
};

template<typename DTYPE, typename FTYPE>
class Solver {
public:
	data_t dot(FTYPE& a, FTYPE& b) {
		data_t result = 0;
		for (auto iter = a.begin(&b, true); iter != a.end(); ++iter) {	// dual iterator
			result += std::get<1>(*iter) * std::get<2>(*iter);
		}
		return result;
	}

	data_t norm(FTYPE& a) {
		data_t result = 0;
		for (auto e = a.begin(true); e != a.end(); ++e)
			result += (*e).second * (*e).second;
		return result;
	}

	// A BiCGStab implementation
	int solve(Matrix<DTYPE,FTYPE>& M, FTYPE& x, FTYPE& b, int max_it, data_t r_tol, data_t a_tol) {

		typedef double T;
		double rho_1, rho_2, alpha, beta, omega, norm_b, norm_r, start_t, end_t;
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

		if (::isnan(norm_r) || ::isnan(norm_b)) {
			cout << "BICGS found nan solution (iter 0)" << norm_b << " " << norm_r << endl;
			norm_r = 1000.;
			exit(4);
		}

		if (norm_r < a_tol) return iter;

		rtilde = r;

		if (debug >= 2)
			printf("BiCGStab Solver started: B-Norm = %g; R-Norm = %g\n", T(norm_b), T(norm_r));

		norm_b = norm_r;
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

			if (::isnan(norm_r)) {
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

	// Steepest Descent solver (alternative to BiCGStab)
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



#endif /* TESTS_SOLVER_HPP_ */
