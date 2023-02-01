#ifndef NEOHOOKEANMATERIAL_H
#define NEOHOOKEANMATERIAL_H


#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/elasticity/kinematics.h>
#include <iostream>
#include <cmath>
#include "NonStandardTensors.h"

using namespace dealii;


template<int dim>
class NeoHookeanMaterial 
{
    public:
				
        NeoHookeanMaterial(const std::string stiffness_case, const double shear_modulud_cortex, const double Poisson, const double stiffness_ratio ,const double max_cell_density, const double subcortix_raduis):
      st_case(stiffness_case), mu_cmax(shear_modulud_cortex), mu_s(shear_modulud_cortex/stiffness_ratio), nu(Poisson) ,c_max(max_cell_density), R_c(subcortix_raduis),
        F_e( Physics::Elasticity::StandardTensors< dim >::I), J_e(1.0)
        {}
        
        
        ~NeoHookeanMaterial(){}
         
	 void update_material_data(const Tensor<2, dim> &F_elastic, const Point<dim> &p, const double &c)
	    {
         double mu_c = 0, dmuc_dc=0;
         if (st_case == "Varying"){
              mu_c = mu_s + ((mu_cmax - mu_s)/(c_max-200)) * ((c > c_max)? (c_max-200):(c-200)) * (c<200? 0:1);
              dmuc_dc = ((mu_cmax - mu_s)/(c_max-200)) * ((c > c_max)? 0:1) * (c<200? 0:1);
         }
         else if (st_case == "Constant"){
             mu_c = mu_cmax;
             dmuc_dc = 0;
         }

         double r = 0;
         if (dim ==2)
             r = p.distance(Point<dim>(0.0,0.0));
         else if (dim ==3)
             r = p.distance(Point<dim>(0.0,0.0, 0.0));

	      H    = std::exp((r-R_c)*20)/(1+std::exp((r-R_c)*20));
              mu        = mu_s+((mu_c-mu_s)*H);
	      dmu_dc    = dmuc_dc*H;
              Lamda     = (2*mu*nu)/(1-(2*nu));
	      dLamda_dc = (2*nu* dmu_dc)/(1-(2*nu));

	      F_e = F_elastic;
	      b_e = Physics::Elasticity::Kinematics::b(F_e);
	      J_e = determinant(F_e);

	     Assert(J_e > 0, ExcInternalError());
	    }
    
   				
	  SymmetricTensor<2, dim> get_Cauchy_stress() const
	    {
	      
	      SymmetricTensor<2, dim> First_term;
	      SymmetricTensor<2, dim> Second_term;
	      double cof = (Lamda* log(J_e))-mu;
	      First_term = cof* Physics::Elasticity::StandardTensors< dim >::I;
	      Second_term = mu* b_e;
	      
	      return (1/J_e)*(First_term + Second_term);
	    }


           double get_J_e() const
	    {
	      return J_e;
	    }

          SymmetricTensor<4, dim> get_tangent_tensor() const
	    {

	     
	      double cof =2*( mu-(Lamda * std::log(J_e)))/J_e;

	      return (cof *  Physics::Elasticity::StandardTensors< dim >::S + (Lamda/J_e) * Physics::Elasticity::StandardTensors< dim >::IxI);
			     
	    }

          Tensor<4 ,dim> get_tangent_tensor_elastic() const 
            {
               Tensor<4 ,dim> c;
               Tensor<2 ,dim> inv_F_e = invert(F_e);
               Tensor<2 ,dim> tr_inv_F_e = transpose(inv_F_e);     
	       Tensor<2 ,dim> I = Physics::Elasticity::StandardTensors< dim >::I;       
               double cof = mu-(Lamda*std::log(J_e));

               c = mu* NonStandardTensors::OverlineDyads<dim>(I,I) + cof* NonStandardTensors::underlineDyads<dim>(tr_inv_F_e,inv_F_e) + Lamda* outer_product(tr_inv_F_e, tr_inv_F_e); 

                   return c;
            } // dp_e / dF_e

	Tensor<2, dim> get_stress_dervitive_wrt_cell() const
	   {
	       Tensor<2 ,dim> tr_inv_F_e = transpose(invert(F_e));
	       double cof = dLamda_dc * std::log(J_e);
		
		return (cof * tr_inv_F_e + dmu_dc * F_e - dmu_dc * tr_inv_F_e);
		} // dp_e / d_c

           double get_elastic_modulus() const {return mu;}

    protected:

    const std::string st_case;
    const double mu_cmax;
	const double mu_s;
    const double nu;
	const double c_max;
    double r;
    const double R_c;
    double H;

        double mu, dmu_dc, Lamda, dLamda_dc;

        Tensor<2, dim> F_e;
        SymmetricTensor <2,dim> b_e;
        double J_e;
       

};


#endif
