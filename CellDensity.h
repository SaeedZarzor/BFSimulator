#ifndef CELLDENSITY_H
#define CELLDENSITY_H

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace dealii;

template<int dim>
class CellDensity
{
     public:
       
	CellDensity(const double &cell_dvision_rate_v, const double &cell_dvision_rate_ovz, const double &cell_migration_speed,  
    const double &diffusivity, const double &cell_migration_threshold, const double &exponent, const double &damention_ratio, const double &V_raduis, const double &subcortix_raduis,
        const double &IVZ_raduis, const double &OVZ_raduis, const double &MST_factor, const std::string &OSVZ_divion_distr):
	G_c_v(cell_dvision_rate_v), G_c_ovz(cell_dvision_rate_ovz),
        v(cell_migration_speed), 
	d_cc(diffusivity), c_0(cell_migration_threshold), 
        gamma(exponent), alpha(damention_ratio), 
        r_v(V_raduis), R_c(subcortix_raduis), 
	r_ivz(IVZ_raduis),r_osvz(OVZ_raduis),c_mst(MST_factor),
    OSVZ_varying(OSVZ_divion_distr),
    first_cell_flux_term(Tensor<1, dim>()),
    second_cell_flux_term(Tensor<1, dim>())
        {}

        ~CellDensity(){}


        void update_flux(const Tensor<2, dim> &F, const Tensor<1, dim> &grad_c, const double &c,const Point<dim> &p )
          {
            
               Tensor<1 ,dim> N = direction_vector(p);
               Tensor<1 ,dim> n = F*N;
               Tensor<1 ,dim> speed = compute_speed(F,c, p);
               Tensor<2 ,dim> I = Physics::Elasticity::StandardTensors< dim >::I;
               grad_c_s = grad_c * invert(F);
               
                        //for(unsigned int i =0; i< dim; i++ )
                          //    speed[i] *= ((i==(dim-1))? alpha : 1);

		       double cof = heaviside_function((c-c_0), gamma) + c * heaviside_function_derivative((c-c_0), gamma);
	
               diffusion_tensor = get_dcc_r(p) * I; //dq_dgrad_c
	
        		           

               first_cell_flux_term  = -c * speed;
               second_cell_flux_term = diffusion_tensor* grad_c_s; // q

		       cell_flux_derivative = -cof * (speed/heaviside_function((c-c_0),gamma));        // dq_dc
                    

               {
	  	      double cof_2 = -c * heaviside_function((c-c_0), gamma) * (get_v_r(p)/n.norm());
		      double a = 1/std::pow((n.norm()),2);

		      Tensor<2, dim> nxn = outer_product(n,n);
		      Tensor<2, dim> direction_derivative = I - (a * nxn);
		      Tensor<3, dim> seconed_term = - outer_product(grad_c_s, diffusion_tensor);

		      cell_flux_deformation_derivative = cof_2 * outer_product(direction_derivative, n) + seconed_term;
                    } // compute dq_dF_Ft
           }
      
    double compute_denisty_source(const Point<dim> &p, const double s)
         {
            double r = 0;
            double G_c_t = 0 ;
        
            if (dim ==2) {
                r= p.distance(Point<dim>(0.0,0.0));
                G_c_t = G_c_v - G_c_v * ((s<=1.8)? (s-1):0.8);
            }
            if (dim==3){
                r = p.distance(Point<dim>(0.0,0.0, 0.0));
                G_c_t = G_c_v;
            }

        
        double G_c_r = G_c_t*(1-heaviside_function((r-r_v),50));
               
        return G_c_r;
      
               }

  double compute_denisty_source_ORGC(const Point<dim> &p, const double time, const double s)
    {
      double alpha = std::abs(std::atan(p[1]/p[0]));
      double sin=1.0;
      if (OSVZ_varying == "Constant")
          sin = 1;
      
      else if (OSVZ_varying == "Linear-gradient" )
          sin = std::sin(alpha);
      
      else if (OSVZ_varying =="Quadratic-gradient")
          sin = std::pow(alpha,2);
      
      else if (OSVZ_varying == "Random1" )
          sin = std::sin(10*alpha) * (alpha);
      
      else if (OSVZ_varying == "Random2")
          sin = (std::sin(20*alpha) * std::log(1/alpha)) - (std::sin(10*alpha) *std::abs(alpha));
      
      else
          ExcMessage("'The OSVZ regional variation' Entery unknown ");

      if (sin < 0)
          sin = 0;
      
      double r= 0;
      double G_c_t= 0;
      double r_ovz_t = r_ivz+ c_mst * ((time > 0)? time:0);
      r_ovz_t = ((r_ovz_t> r_osvz)?  r_osvz:r_ovz_t);
      if (dim ==2){
          r = p.distance(Point<dim>(0.0,0.0));
          G_c_t = G_c_ovz - G_c_ovz * ((s<=1.8)? (s-1):0.8);
      }
      else if (dim ==3){
          r = p.distance(Point<dim>(0.0,0.0, 0.0));
          G_c_t = G_c_ovz;
      }
      double G_c_r = G_c_t * sin  * (heaviside_function((r-r_ivz),50)-heaviside_function((r-r_ovz_t),50));
                  return G_c_r;
     
               }
    
  Tensor<1, dim> compute_speed(const Tensor<2, dim> &F ,const double &c ,const Point<dim> &p )
    {
        
        double v_r    = get_v_r(p);

       Tensor<1 ,dim> N = direction_vector(p);
       Tensor<1 ,dim> n = F*N;
       Tensor<1 ,dim> speed = v_r * (n/n.norm());
        return (heaviside_function((c-c_0),gamma) * speed);
    }
        
        double  get_v_r(const Point<dim> &p)  {
            double r =0;
            if (dim ==2)
             r = p.distance(Point<dim>(0.0,0.0));
            else if (dim ==3)
                r = p.distance(Point<dim>(0.0,0.0, 0.0));
            return (v*(1-heaviside_function((r-R_c),10)));}
    
        double  get_dcc_r(const Point<dim> &p)  {
            double r =0;
            if (dim ==2)
             r = p.distance(Point<dim>(0.0,0.0));
            else if (dim ==3)
                r = p.distance(Point<dim>(0.0,0.0, 0.0));
            return (d_cc*(heaviside_function((r-R_c),10)));}
        
        Tensor<1, dim> get_first_flux_term() const {return first_cell_flux_term;}
        Tensor<1, dim> get_second_flux_term() const {return second_cell_flux_term;}
        Tensor<1, dim> get_flux() const {return (first_cell_flux_term + second_cell_flux_term);}
        Tensor<1 ,dim> get_flux_derivative() const {return cell_flux_derivative;}
        Tensor<1, dim> get_grad_c_s() const {return grad_c_s;}
        Tensor<2 ,dim> get_diffusion_tensor() const {return diffusion_tensor;}
        Tensor<3 ,dim> get_flux_deformation_derivative() const {return cell_flux_deformation_derivative;}


     protected:
	
        double G_c_v;
        double G_c_ovz;
        double v;
        double d_cc;
        double c_0;
        double gamma;
        double alpha;
        double r_v;
        double R_c;
        double r_ivz;
        double r_osvz;
        double c_mst;

        std::string OSVZ_varying;
        Tensor<1 ,dim> grad_c_s;
        Tensor<1, dim> first_cell_flux_term;
        Tensor<1, dim> second_cell_flux_term;
        Tensor<1 ,dim> cell_flux_derivative;              // dq / dc
        Tensor<2 ,dim> diffusion_tensor;                  // d^cc Tensor
        Tensor<3, dim> cell_flux_deformation_derivative;  // dq / dF * F^t      


        Tensor<1 ,dim> direction_vector(const Point<dim> &p)
          {
               return Tensor<1, dim>(p/p.norm());
              }       // n= x/|x|


       double heaviside_function(const double &x, const double c)  // H(c-c_0; gamma)
        {
             
             double a = std::exp(c*x);
             double b = 1.+std::exp(c*x);

            return a/b;
        
          }

      double heaviside_function_derivative(const double &x, const double c)  // dH_dc
        {
            
              double a = c * std::exp(c*x);
              double b = std::pow((1.+std::exp(c*x)),2);
                return a/b;
           
             }


};

#endif
