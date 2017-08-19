/*
 * numerics.hpp
 *	
 *	Numerical building blocks that operate on Field() objects.
 *	Also defines useful typedefs to avoid template usage. IE DenseScalarField2 is a 2D non-sparse scalar field.
 *
 *	Example:
 *	  DenseVectorField3 grad_of_x;
 *	  grad_of_x.createEntireLevel(8);
 *	  grad<3>(x);  // x is instance of ScalarField3
 *	
 *  Created on: Dec 27, 2016
 *      Author: Christian Huettig
 */
#pragma once

// Own includes
/*
#include "hcs.hpp"
#include "tensor.hpp"
#include "field.hpp"
#include "sparsefield.hpp"
#include "densefield.hpp"
*/

typedef HCS<1> H1;  // 1D
typedef HCS<2> H2;  // 2D
typedef HCS<3> H3;  // 3D
typedef HCS<4> H4;  // 4D
typedef HCS<5> H5;  // 5D

typedef Tensor1<data_t, 2> Vec2;	// Single "vector" in 2D / 3D. Similar to old Point<T>
typedef Tensor1<data_t, 3> Vec3;
typedef Tensor1<data_t, 4> Vec4;
typedef Tensor1<data_t, 5> Vec5;

typedef Field<data_t, H1> ScalarField1; // 1D scalar field type
typedef Field<data_t, H2> ScalarField2;
typedef Field<data_t, H3> ScalarField3;
typedef Field<data_t, H4> ScalarField4;
typedef Field<data_t, H5> ScalarField5;

typedef Field<data_t, H1> VectorField1; // 1D scalar field type
typedef Field<Vec2, H2> VectorField2;
typedef Field<Vec3, H3> VectorField3;
typedef Field<Vec4, H4> VectorField4;
typedef Field<Vec5, H5> VectorField5;

typedef DenseField<data_t, H1> DenseScalarField1; // 1D scalar field type
typedef DenseField<data_t, H2> DenseScalarField2;
typedef DenseField<data_t, H3> DenseScalarField3;
typedef DenseField<data_t, H4> DenseScalarField4;
typedef DenseField<data_t, H5> DenseScalarField5;

typedef DenseField<data_t, H1> DenseVectorField1; // 1D scalar field type
typedef DenseField<Vec2, H2> DenseVectorField2;
typedef DenseField<Vec3, H3> DenseVectorField3;
typedef DenseField<Vec4, H4> DenseVectorField4;
typedef DenseField<Vec5, H5> DenseVectorField5;

typedef SparseField<data_t, H1> SparseScalarField1; // 1D scalar field type
typedef SparseField<data_t, H2> SparseScalarField2;
typedef SparseField<data_t, H3> SparseScalarField3;
typedef SparseField<data_t, H4> SparseScalarField4;
typedef SparseField<data_t, H5> SparseScalarField5;

typedef SparseField<data_t, H1> SparseVectorField1; // 1D scalar field type
typedef SparseField<Vec2, H2> SparseVectorField2;
typedef SparseField<Vec3, H3> SparseVectorField3;
typedef SparseField<Vec4, H4> SparseVectorField4;
typedef SparseField<Vec5, H5> SparseVectorField5;




// Dimension-independent gradient operator
// Takes a scalar field as input and returns the vector field
// Example:  VectorField3 grad_of_x = grad<3>(x);  // x is instance of ScalarField3
template<int dimension>
void grad(Field<data_t, HCS<dimension> > &source, Field<Tensor1<data_t, dimension>, HCS<dimension> > &result) {
	typedef Tensor1<data_t, dimension> Vec;
	typedef HCS<dimension> HX;
	typedef Field<Vec, HX> VectorField;
    typedef Field<data_t, HX> ScalarField;

	// strange calling convention because template depends now on another template
	result.template convert<data_t>(source, [](coord_t c, ScalarField &source)->Vec {

		HX &h = source.hcs;
		Vec gradient;

		// Gradient calculation with finite-difference stencil
		level_t l = h.GetLevel(c);
		data_t dist = 4 * (h.scales[0] / data_t(1U << l)); // neighbor distance at that level, assuming all scales equal, * 2 because gradient is 2nd order
		for (int n_idx = 0; n_idx < h.parts; n_idx++) {  // Traverse all neighbors
			coord_t c_ne = h.getNeighbor(c, n_idx);
			data_t ne_val = source.get(c_ne); // will respect boundary condition
			gradient[n_idx >> 1] += (n_idx & 1 ? -ne_val : ne_val) / dist;
		}
		return gradient;
	});
	result.propagate();
}

// Dimension-independent divergence operator
// Takes a scalar field as input and returns the vector field
// Example:  VectorField3 grad_of_x = grad<3>(x);  // x is instance of ScalarField3
template<int dimension>
void div(Field<Tensor1<data_t, dimension>, HCS<dimension> >  &source, Field<data_t, HCS<dimension> > &result) {
	typedef Tensor1<data_t, dimension> Vec;
	typedef HCS<dimension> HX;
	typedef Field<data_t, HX> ScalarField;
    typedef Field<Vec, HX> VectorField;

	// strange calling convention because template depends now on another template
	result.template convert<Vec>(source, [](coord_t c, VectorField &source)->data_t {

		HX &h = source.hcs;
		data_t divergence = 0;

		// Gradient calculation with finite-difference stencil
		level_t l = h.GetLevel(c);
		data_t dist = 4 * (h.scales[0] / data_t(1U << l)); // neighbor distance at that level, assuming all scales equal, * 2 because gradient is 2nd order
		for (int n_idx = 0; n_idx < h.parts; n_idx++) {  // Traverse all neighbors
			coord_t c_ne = h.getNeighbor(c, n_idx);
			data_t ne_val = source.get(c_ne)[n_idx >> 1]; // will respect boundary condition
			divergence += (n_idx & 1 ? -ne_val : ne_val) / dist;
		}
		return divergence;
	});
	result.propagate();
}


