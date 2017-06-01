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
 *	  cout << t2.normalized() << " Length: " << t2.normalized().length() << endl; // normalized
 *	
 *  Created on: Dec 27, 2016
 *      Author: Christian Huettig
 */
#pragma once

using namespace std;

// T == elemental data type, D == dimensions
template <typename T, unsigned char D>
class Tensor1 {
// Attributes
public:
	// This wastes storage space for D < 3... But I would miss the comfort.
	// Better Ideas?
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

	// Constructors
	Tensor1<T, D>() {}
	Tensor1<T, D>(const array<T, D> &ia) { value = ia; }
	Tensor1<T, D>(initializer_list<T> il) { copy(il.begin(), il.end(), value.begin()); }
	Tensor1<T, D>(T s) { for (auto &e : value) e = s;}


	T norm()			  		{ return (*this) * (*this);		}
	T length()					{ return sqrt(norm());			}
	T magnitude()				{ return sqrt(norm());			}
	T mag()						{ return sqrt(norm());			}
	void normalize()			{ (*this) *= 1./length(); 		}
	Tensor1<T, D> normalized()	{ return Tensor1<T, D>((*this) /= length()); }

	/* TODO
	Tensor1<T, n> Sph2Cart();
	Tensor1<T, n> Cart2Sph();
	*/

	Tensor1<T, D>& operator=(const Tensor1<T, D> &P)		{ value = P.value;	return *this;	} // assign
	Tensor1<T, D>& operator=(const T &val)  				{ for (auto &e : value) e = val;	return *this;		} // assign

	Tensor1<T, D>& operator*= (const Tensor1<T, D>& rhs)	{ for (int i = 0; i < D; i++) value[i] *= rhs.value[i]; return *this;}
	Tensor1<T, D>& operator/= (const Tensor1<T, D>& rhs)	{ for (int i = 0; i < D; i++) value[i] /= rhs.value[i]; return *this;}
	Tensor1<T, D>& operator+= (const Tensor1<T, D>& rhs)	{ for (int i = 0; i < D; i++) value[i] += rhs.value[i]; return *this;}
	Tensor1<T, D>& operator-= (const Tensor1<T, D>& rhs)	{ for (int i = 0; i < D; i++) value[i] -= rhs.value[i]; return *this;}

	Tensor1<T, D>& operator*= (const T& val) 				{ for (auto &e : value) e *= val;	return *this;		}
	Tensor1<T, D>& operator/= (const T& val) 				{ for (auto &e : value) e /= val;	return *this;		}
	Tensor1<T, D>& operator+= (const T& val) 				{ for (auto &e : value) e += val;	return *this;		}
	Tensor1<T, D>& operator-= (const T& val) 				{ for (auto &e : value) e -= val;	return *this;		}

	T& 			   operator[] (const unsigned i)			{ return value[i];	} // index access to x, y, z

	// cross-product for 2 or 3 dimensions
	Tensor1<T, 3> operator^(const Tensor1<T, D>& P) {
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
	
	// for use as output stream cout << myT1...
	friend ostream& operator<< (std::ostream& os, const Tensor1<T, D>& t) {
		os << "T1(";
		for (auto e : t.value)
			os << e << ", ";
		os << "[" << (int)D << "]) ";
		return os;
	}

};

// All non-member operator overloads to allow lhs/rhs exchange (3 * t1 same as t1 * 3)

// "Multiply" of two Tensor1 leads to dot-product. Component-wise only with *=
template <typename T, unsigned char D> T operator* (const Tensor1<T, D>& lhs, const Tensor1<T, D>& rhs)
{ T result = 0; for (int i = 0; i < D; i++) result += lhs.value[i] * rhs.value[i]; return result;}

// Compiler will complain about ambiguity
//template <typename T, unsigned char D> Tensor1<T, D> operator* (const Tensor1<T, D>& lhs, const Tensor1<T, D>& rhs)
//{ Tensor1<T, D> result = lhs; result *= rhs; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator* (const T& val, const Tensor1<T, D>& rhs)
{ Tensor1<T, D> result = rhs; result *= val; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator* (const Tensor1<T, D>& lhs, const T& val)
{ Tensor1<T, D> result = lhs; result *= val; return result;}

template <typename T, unsigned char D> Tensor1<T, D> operator/ (const Tensor1<T, D>& lhs, const Tensor1<T, D>& rhs)
{ Tensor1<T, D> result = lhs; result /= rhs; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator/ (const T& val, const Tensor1<T, D>& rhs)
{ Tensor1<T, D> result = rhs; result /= val; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator/ (const Tensor1<T, D>& lhs, const T& val)
{ Tensor1<T, D> result = lhs; result /= val; return result;}

template <typename T, unsigned char D> Tensor1<T, D> operator+ (const Tensor1<T, D>& lhs, const Tensor1<T, D>& rhs)
{ Tensor1<T, D> result = lhs; result += rhs; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator+ (const T& val, const Tensor1<T, D>& rhs)
{ Tensor1<T, D> result = rhs; result += val; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator+ (const Tensor1<T, D>& lhs, const T& val)
{ Tensor1<T, D> result = lhs; result += val; return result;}

template <typename T, unsigned char D> Tensor1<T, D> operator- (const Tensor1<T, D>& lhs, const Tensor1<T, D>& rhs)
{ Tensor1<T, D> result = lhs; result -= rhs; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator- (const T& val, const Tensor1<T, D>& rhs)
{ Tensor1<T, D> result = rhs; result -= val; return result;}
template <typename T, unsigned char D> Tensor1<T, D> operator- (const Tensor1<T, D>& lhs, const T& val)
{ Tensor1<T, D> result = lhs; result -= val; return result;}

/*
template <typename T, unsigned char D> valarray<bool> operator== (const Tensor1<T, D>& lhs, const Tensor1<T, D>& rhs);
template <typename T, unsigned char D> valarray<bool> operator== (const T& val, const Tensor1<T, D>& rhs);
template <typename T, unsigned char D> valarray<bool> operator== (const Tensor1<T, D>& lhs, const T& val);

template <typename T, unsigned char D> valarray<bool> operator!= (const Tensor1<T, D>& lhs, const Tensor1<T, D>& rhs);
template <typename T, unsigned char D> valarray<bool> operator!= (const T& val, const Tensor1<T, D>& rhs);
template <typename T, unsigned char D> valarray<bool> operator!= (const Tensor1<T, D>& lhs, const T& val);
*/
