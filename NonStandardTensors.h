#ifndef NONSTANDARDTENSORS_H
#define NONSTANDARDTENSORS_H

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/elasticity/kinematics.h>
#include <iostream>
#include <cmath>

using namespace dealii; 

namespace NonStandardTensors

{
   template<int dim> 
	Tensor<4,dim> OverlineDyads (Tensor<2,dim> &A, Tensor<2,dim> &B)
	  {
	     Tensor <4,dim> C;
              for(unsigned int i=0; i<dim ; ++i)
		  for(unsigned int j=0; j<dim ; ++j)
		     for(unsigned int k=0; k<dim ; ++k)
			  for(unsigned int l=0; l<dim ; ++l)
				C[i][j][k][l] = A[i][k]*B[j][l];

		return C;
	    }


  template<int dim> 
	Tensor<4,dim> underlineDyads (Tensor<2,dim> &A, Tensor<2,dim> &B)
	  {
	     Tensor <4,dim> C ;
              for(unsigned int i=0; i<dim ; ++i)
		  for(unsigned int j=0; j<dim ; ++j)
		     for(unsigned int k=0; k<dim ; ++k)
			  for(unsigned int l=0; l<dim ; ++l)
				C[i][j][k][l] = A[i][l]*B[j][k];

		return C;
	    }

  template<int dim>
        Tensor<2,dim> fourth_second_orders_contrac (Tensor<4,dim> A, Tensor<2,dim> B)
	  {
             Tensor<2,dim> C =Tensor<2,dim>();

	      for(unsigned int i=0; i<dim ; ++i)
		  for(unsigned int j=0; j<dim ; ++j)
		      for(unsigned int k=0; k<dim ; ++k)
			   for(unsigned int l=0; l<dim ; ++l)
				C[i][j] += A[i][j][k][l]*B[k][l];
		
		    return C;
		}

}

#endif
