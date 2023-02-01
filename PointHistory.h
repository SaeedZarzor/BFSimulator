#ifndef POINTHISTORY_H
#define POINTHISTORY_H


#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/physics/elasticity/standard_tensors.h>		
#include "NeoHookeanMaterial.h"
#include "NonStandardTensors.h"
#include "CellDensity.h"
#include "Parameter.h"
#include "Growth.h"
#include <iostream>


using namespace dealii;


template <int dim>
class PointHistory
{
public:
	PointHistory():
        q(Tensor<1, dim>()),
	    F(Physics::Elasticity::StandardTensors< dim >::I),
        inv_F_g(Physics::Elasticity::StandardTensors< dim >::I),
        J(1.0)
	{}
	~PointHistory(){
	       delete material;
	       material = NULL;

               delete growth;
               growth= NULL;

               delete density;
               density= NULL;
             
                       }


	   void setup_lqp (const Parameter::GeneralParameters &parameter, const Point<dim> & position, const double delta_t)
	    {
             
          p = position;
          d_t = delta_t;
	      double R_c = (1-parameter.cortex_thickness)*parameter.initial_radius;
	      double r_v = parameter.ventricular_raduis*parameter.initial_radius;
	      double r_ivz = parameter.subventricular_raduis*parameter.initial_radius;
          double r_osvz = 0.662*parameter.initial_radius;
	      

	      material = new NeoHookeanMaterial<dim>(parameter.stiffness_case, parameter.shear_modulud_cortex, parameter.Poisson, parameter.stiffness_ratio ,parameter.max_cell_density, R_c);

	      growth = new Growth<dim>(parameter.growth_rate, parameter.growth_ratio, parameter.growth_exponent,
                  parameter.cell_migration_threshold,  parameter.damention_ratio, R_c);

              density = new CellDensity<dim>(parameter.cell_dvision_rate_v, parameter.cell_dvision_rate_ovz ,parameter.cell_migration_speed,parameter.diffusivity, parameter.cell_migration_threshold, parameter.exponent, parameter.damention_ratio,
                  r_v,R_c, r_ivz,r_osvz, parameter.MST_factor, parameter.OSVZ_divion_distr);

       
                 update_values (Tensor<2, dim>(), Tensor<1, dim>(), 0 , Tensor<2, dim>(), Tensor<2, dim>(), 0, 0 
                                 ,0 , false, std::vector<double>{1.0,1.0,1.0});
	    }

          void update_values (const Tensor<2, dim> &Grad_u, const Tensor<1, dim> &Grad_c, const double &c, const Tensor<2, dim> &Grad_u_n ,    
                                const Tensor<2, dim> &Grad_u_n_1,const double &c_n, const double &c_n_1,
                                const double &t, bool update, const std::vector<double> &stretch_max)
    {

           F = (Physics::Elasticity::StandardTensors< dim >::I + Grad_u);
	       Tensor<2, dim> F_n   = (Physics::Elasticity::StandardTensors< dim >::I + Grad_u_n);
	       Tensor<2, dim> F_n_1 = (Physics::Elasticity::StandardTensors< dim >::I + Grad_u_n_1);

          J     =  determinant(F);
		  J_n   =  determinant(F_n);
          J_n_1 =  determinant(F_n_1);

          dcc_r = density-> get_dcc_r(p);
          cell_density = c;
          velocity = density -> compute_speed(F,c,p);
		  old_cell_density = c_n;
		  old_old_cell_density = c_n_1;
          old_velocity_values = density-> compute_speed (F_n, c_n, p);
          old_old_velocity_values = density-> compute_speed (F_n_1, c_n_1, p);

                     if(update){ 
                         growth-> update_growth(p, c);
                         Tensor<2,dim> F_g= growth->get_growth_tensor();
        		
       			  inv_F_g=invert(F_g);
                        }
   
                      F_e = F * inv_F_g;
         
                     density-> update_flux(F, Grad_c, c, p);
            
                     material-> update_material_data(F_e, p, c_n);

                     
                    q              = J*(density-> get_flux());
                    q_migration    = J*(density-> get_first_flux_term());
                    dq_dc          = J*(density->get_flux_derivative());
                    dq_dgrad_c     = J*(density-> get_diffusion_tensor());
                    dq_dF_Ft       = J*(density->get_flux_deformation_derivative());
                    tau            = J* (material-> get_Cauchy_stress());
                    elastic_tensor = J*(material-> get_tangent_tensor());
                    dP_e_dF_e      = material-> get_tangent_tensor_elastic();
		            dP_e_dc        = (update ? material-> get_stress_dervitive_wrt_cell():Tensor<2,dim>());
             	    G              = (update ? growth-> get_G():Tensor<2,dim>());
              
                 F_c   = J*(density-> compute_denisty_source(p, stretch_max[0]));
		         F_c_2 = J*(density-> compute_denisty_source_ORGC(p, t, stretch_max[0]));
      
         }


             Tensor<2 ,dim> get_dtau_dc() const 
            {
                Tensor<2 ,dim> dtau_dc;
                double         J_g         = 1/determinant(inv_F_g);
                          
                Tensor<2 ,dim> tr_inv_F_g  = transpose(inv_F_g);
                Tensor<2, dim> trm_t       = transpose(F_e * G * invert(F));  
		        Tensor<2, dim> src         = F_e * G * inv_F_g;

		src = NonStandardTensors::fourth_second_orders_contrac<dim>(dP_e_dF_e, src);
  
                dtau_dc  = scalar_product(tr_inv_F_g, G) * tau;
                dtau_dc -= tau * trm_t;
                dtau_dc -= J_g * src * transpose(F_e);
		        dtau_dc += J_g * dP_e_dc * transpose(F_e);
		
		
                return dtau_dc;
              }

               /*Tensor<2 ,dim> get_dtau_dc() const 
            {
                Tensor<2 ,dim> dtau_dc;
                //double         J_g         = 1/determinant(inv_F_g);
                          
                Tensor<2 ,dim> tr_inv_F_g  = transpose(inv_F_g);
                Tensor<2, dim> trm         = F_e * G * invert(F);  
		Tensor<2, dim> src         = F_e * G * inv_F_g;
  
                dtau_dc  = scalar_product(tr_inv_F_g, G) * material-> get_Cauchy_stress();
                dtau_dc -= trm * material-> get_Cauchy_stress();
                dtau_dc -= NonStandardTensors::fourth_second_orders_contrac<dim>((material-> get_tangent_tensor()), src);
		
                return (J*dtau_dc);
              }*/

           
           SymmetricTensor<4, dim> get_elastic_tensor() const {return elastic_tensor;}
           SymmetricTensor<2, dim> get_tau() const {return tau;}
    
           Tensor<3 ,dim> get_Jdq_dF_Ft() const {return dq_dF_Ft;}

           Tensor<2, dim> get_inv_F() const {return invert(F);}
           Tensor<2 ,dim> get_Jdq_dgrad_c() const {return dq_dgrad_c;}
    
           Tensor<1, dim> get_grad_c_spatial() const {return (density->get_grad_c_s());}
           Tensor<1 ,dim> get_Jq() const {return q;}
           Tensor<1 ,dim> get_Jq_migration() const {return q_migration;}
           Tensor<1, dim> get_Jdq_dc() const {return dq_dc;}
           Tensor<1, dim> get_velocity() const {return velocity;}
           Tensor<1, dim> get_old_velocity_values() const {return old_velocity_values;}
           Tensor<1, dim> get_old_old_velocity_values() const {return old_old_velocity_values;}

           double  get_dcc_r() const{return dcc_r;}
           double  get_J() const {return J;}
	       double  get_J_old() const {return J_n;}
           double  get_J_old_old() const {return J_n_1;}
           double  get_c() const {return cell_density;}
     	   double  get_c_old() const {return old_cell_density;}
	       double  get_c_old_old() const {return old_old_cell_density;}
           double  get_F_c() const {return F_c;}
           double  get_F_c_2() const {return F_c_2;}

           double  get_growth_factor_t() const {return (growth-> get_growth_factor_tangent());}
           double  get_growth_factor_r() const {return (growth-> get_growth_factor_radius());}
           double  get_growth_norm() const {  
   
                 if(dim == 3)
                   return ((growth-> get_growth_tensor().norm())/std::sqrt(3));

                 else if(dim == 2)
                   return ((growth-> get_growth_tensor().norm())/std::sqrt(2));
                  }

          
          double get_elastic_modulus() const {return (material-> get_elastic_modulus());}

	private:
	NeoHookeanMaterial<dim> *material;
    Growth<dim>             *growth;
    CellDensity<dim>        *density;

	Point<dim>  p;
    Tensor<1, dim> q;
    Tensor<1, dim> q_migration;
    Tensor<1, dim> dq_dc;
    Tensor<1, dim> velocity;
    Tensor<1, dim> old_velocity_values;
    Tensor<1, dim> old_old_velocity_values;

	Tensor<2, dim> F;
    Tensor<2 ,dim> inv_F_g;
    Tensor<2, dim> F_e;
    Tensor<2 ,dim> dq_dgrad_c;
	Tensor<2 ,dim> dP_e_dc;
	Tensor<2, dim> G ;
	
    Tensor<3, dim> dq_dF_Ft;
    SymmetricTensor<2 ,dim> tau;
    SymmetricTensor<4, dim> elastic_tensor;
    Tensor<4, dim> dP_e_dF_e;


    double dcc_r;
    double d_t;
	double J;
	double J_n;
	double J_n_1;
    double F_c;
	double F_c_2;
    double cell_density;
	double old_cell_density;
	double old_old_cell_density;

        

};


#endif
