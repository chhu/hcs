/*
 * tensor.hpp
 *	
 *	An quick-and-dirty arithmetic wrapper class for vectors (tensors of rank 1, therefore called Tensor1)
 *	The first template argument specifies the data (float) type, the second the dimensions.
 *
 *	Example:
 *	  Tensor1<double, 3> t1({1,2,3});
 *	  Tensor1<double, 3> t2({9,8,7});
 *	  cout << t1.length() << endl;
 *	  cout << t1 * t2 << endl; // dot product
 *	  cout << t1 ^ t2 << endl; // cross product
 *	  cout << t2 * 5  << endl; // scale (5 * t2 does not work because double does not know how to multiply Tensor1)
 *	  cout << t2.normalized() << " " << t2.normalized().length() << endl; // normalized
 *	
 *  Created on: Dec 27, 2016
 *      Author: Christian Huettig
 */
#pragma once

using namespace std;

template <typename T, unsigned char D>
class Tensor1 {
// Attributes
public:
	union {
		struct {
			T x;
			T y;
			T z;
		};
		struct {
			T u;
			T v;
			T w;
		};
		struct {
			T r;
			T phi;
			T theta;
		};
		array<T, D> value{};  // braces guarantee zero init
	};

	// zeroes a vector of Point<T>
	static void ClearVector(vector<Tensor1<T, D> > &vec) {
		for (int i = 0; i < vec.size(); i++)
			vec[i] = 0;
	}

	// Constructors
	Tensor1<T, D>() {}
	Tensor1<T, D>(initializer_list<T> il) { copy(il.begin(), il.end(), value.begin()); }
	Tensor1<T, D>(T s) { for (auto &e : value) e = s;}


	T norm()			    { return (*this) * (*this);		}
	T length()			{ return sqrt(norm());			}
	void normalize()			{ (*this) *= 1./length(); 		}
	Tensor1<T, D> normalized()	{ return Tensor1<T, D>((*this) /= length()); }

	/* TODO
	Tensor1<T, n> Sph2Cart();
	Tensor1<T, n> Cart2Sph();
	*/

	Tensor1<T, D>& operator=(const Tensor1<T, D> &P)	{ value = P.value;	return *this;	} // assign
	Tensor1<T, D>& operator=(const T &f)  			{ for (auto &e : value) e = f;	return *this;		} // assign

	Tensor1<T, D> operator+(Tensor1<T, D> P) {
		Tensor1<T, D> result;
		for (int i = 0; i < D; i++)
			result.value[i] = value[i] + P.value[i];
		return result;
	} // add

	Tensor1<T, D>& operator+=(Tensor1<T, D> P){
		for (int i = 0; i < D; i++)
			value[i] += P.value[i];
		return *this;
	} // add +=

	Tensor1<T, D>& operator+=(T f)  			{
		for (int i = 0; i < D; i++)
			value[i] += f;
		return *this;
	} // add +=

	Tensor1<T, D> operator-(Tensor1<T, D> P) {
		Tensor1<T, D> result;
		for (int i = 0; i < D; i++)
			result.value[i] = value[i] - P.value[i];
		return result;
	} // subtract

	Tensor1<T, D> operator-()				{
		Tensor1<T, D> result;
		for (int i = 0; i < D; i++)
			result.value[i] = -value[i];
		return result;
	} // unary -

	Tensor1<T, D>& operator-=(Tensor1<T, D> P){
		for (int i = 0; i < D; i++)
			value[i] -= P.value[i];
		return *this;
	} // subtract -=


	Tensor1<T, D>& operator-=(T f)  			{
		for (int i = 0; i < D; i++)
				value[i] -= f;
		return *this;
	} // subtract -=

	T        	 operator*(Tensor1<T, D> P) {
		T total = 0;
		for (int i = 0; i < D; i++)
			total += value[i] * P.value[i];
		return total;
	} // dot product

	/* cannot overload for some reason
	inline Tensor1<T, n> operator*(Tensor1<T, n> P) {
		Tensor1<T, n> result;
		for (int i = 0; i < n; i++)
			result[i] = value[i] * P.value[i];
		return result;
	} // dot product
	 */

	Tensor1<T, D> operator*(T f)     		{
		Tensor1<T, D> result;
		for (int i = 0; i < D; i++)
			result.value[i] = value[i] * f;
		return result;
	} // scalar product

	Tensor1<T, D>& operator*=(T f)    		{
		for (int i = 0; i < D; i++)
			value[i] *= f;
		return *this;
	} // scalar mult *=

	Tensor1<T, D>& operator*=(Tensor1<T, D> P){
		for (int i = 0; i < D; i++)
			value[i] *= P[i];
		return *this;
	} // scalar mult *=

	Tensor1<T, D> operator/(T f)     		{
		Tensor1<T, D> result;
		for (int i = 0; i < D; i++)
			result[i] = value[i] / f;
		return result;
	} // scalar div

	Tensor1<T, D>& operator/=(T f)    		{
		for (int i = 0; i < D; i++)
			value[i] /= f;
		return *this;
	} // scalar div /=

	Tensor1<T, D> operator/(Tensor1<T, D> P) {
		Tensor1<T, D> result;
		for (int i = 0; i < D; i++)
			result[i] = value[i] / P[i];
		return result;
	} // scalar div

	Tensor1<T, D>& operator/=(Tensor1<T, D> P){
		for (int i = 0; i < D; i++)
			value[i] /= P[i];
		return *this;
	} // scalar div /=

	inline Tensor1<T, 3> operator^(Tensor1<T, D> P) {
		if (D != 2 && D !=3)
			throw range_error("Tensor1 cross product N mismatch");
		Tensor1<T, 3> result;
		result.z = x * P.y - y * P.x;
		if (D == 3) {
			result.x = y * P.z - z * P.y;
			result.y = z * P.x - x * P.z;
		}
		return result;
	} // cross product
	
	inline T 	operator[] (unsigned i)	{ return value[i];	} // index access to x, y, z
	
	bool    operator==(Tensor1<T, D> P) {
		bool result = true;
		for (int i = 0; i < D; i++)
			result &= value[i] == P[i];
		return result;
	}

	bool    operator==(T f) {
		bool result = true;
		for (int i = 0; i < D; i++)
			result &= value[i] == f;
		return result;
	}

	bool    operator!=(Tensor1<T, D>& P) { return !(*this == P); 			} // is not equal to?

	friend ostream& operator<< (std::ostream& os, const Tensor1<T, D>& t) {
		os << "T1(";
		for (auto e : t.value)
			os << e << ", ";
		os << "[" << D << "]) ";
		return os;
	}

};



// matmul nxn dense matrix as valarray in A with x = result; no checks.
template <typename T>
void gemv(int n, valarray<T> &A, valarray<T> &x, valarray<T> &result) {
	result.resize(n, T(0));
	for (int row = 0; row < n; row++)
		for (int col = 0; col < n; col++)
			result[row] += A[col + row * n] * x[col];
}

template <typename T>
int dmat_solve(int n, int rhs_num, T a[]) {
  T apivot, factor, temp;
  int i, ipivot, j, k;

  for ( j = 0; j < n; j++ )  {
    ipivot = j;
    apivot = a[j+j*n];
    for ( i = j; i < n; i++ ) 
      if ( abs ( apivot ) < abs ( a[i+j*n] ) ) {
        apivot = a[i+j*n];
        ipivot = i;
      }
                                                                                                                                                             
    if ( apivot == 0.0 ) 
      return j;
//
//  Interchange.
//
    for ( i = 0; i < n + rhs_num; i++ ) {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
//
//  A(J,J) becomes 1.
//
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
        a[j+k*n] = a[j+k*n] / apivot;
//
//  A(I,J) becomes 0.
//
    for ( i = 0; i < n; i++ )
      if ( i != j ) {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
            a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
      }
  }                                                                                                                                                             
  return 0;
}



#endif /* TENSOR_HPP_ */
