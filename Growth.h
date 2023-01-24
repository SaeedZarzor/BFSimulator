#ifndef GROWTH_H
#define GROWTH_H

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <iostream>
#include <cmath>

using namespace dealii;

template <int dim>
  class StandardTensors
  {
    
  public:
    // $\mathbf{I}$
    static const SymmetricTensor<2, dim> I;
    // $\mathbf{I} \otimes \mathbf{I}$
    static const SymmetricTensor<4, dim> IxI;
    // $\mathcal{S}$, note that as we only use this fourth-order unit tensor
    // to operate on symmetric second-order tensors.  To maintain notation
    // consistent with Holzapfel (2001) we name the tensor $\mathcal{I}$
    static const SymmetricTensor<4, dim> II;
    // Fourth-order deviatoric tensor such that
    // $\textrm{dev} \{ \bullet \} = \{ \bullet \} -
    //  [1/\textrm{dim}][ \{ \bullet\} :\mathbf{I}]\mathbf{I}$
    static const SymmetricTensor<4, dim> dev_P;
    
    Tensor<1,dim> Na_initialize() ;
    Tensor<1,dim> Nb_initialize() ;
    Tensor<1,dim> Nc_initialize() ;
    
  };
 

  template <int dim>
  const SymmetricTensor<2, dim>
  StandardTensors<dim>::I = unit_symmetric_tensor<dim>();

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::IxI = outer_product(I, I);

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::II = identity_tensor<dim>();

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::dev_P = deviator_tensor<dim>();

  
  template <int dim>
  Tensor<1,dim> StandardTensors<dim>::Na_initialize() 
  {
      Point<dim> value;
      for(unsigned int i=0; i<dim ;++i)
          value[i]=((i==0)? 1:0);
      return value;
  }
      
    template <int dim>
  Tensor<1,dim> StandardTensors<dim>::Nb_initialize() 
  {
      Point<dim> value;
      for(unsigned int i=0; i<dim ;++i)
          value[i]=((i==1)? 1:0);
      return value;
  }
  
    template <int dim>
  Tensor<1,dim> StandardTensors<dim>::Nc_initialize()
  {
      Point<dim> value;
      for(unsigned int i=0; i<dim ;++i)
          value[i]=((i==2 && dim==3)? 1:0);
      return value;
  }


template <int dim>
class Growth {
   public:
	Growth(const double &growth_rate, const double &growth_ratio, const double &growth_exponent ,const double &cell_threshold, const double &damention_ratio, const double &subcortix_raduis):
        k_s(growth_rate), b(growth_ratio), alpha_g(growth_exponent), c_0(cell_threshold), 
        alpha(damention_ratio),R_c(subcortix_raduis), v_t(1.0),v_r(1.0), 
        G(Tensor<2, dim>()), F_g(Physics::Elasticity::StandardTensors< dim >::I){}
        ~Growth(){}




      void update_growth( const Point<dim> &p , const double &c)
         {
             
   
              //StandardTensors<dim>   StandardTensor;
              //Tensor<1,dim> Na = StandardTensor.Na_initialize();
	      //Tensor<1,dim> Nb = StandardTensor.Nb_initialize();
	      //Tensor<2,dim> tangent_direction = outer_product(Na,Na);
              //Tensor<2, dim> normal_direction = outer_product(Nb,Nb);

	      Tensor<1, dim> N = (Tensor<1, dim>(p/p.norm())); 
              Tensor<2, dim> normal_direction = outer_product(N,N);
              Tensor<2, dim> tangent_direction = Physics::Elasticity::StandardTensors< dim >::I - normal_direction;
               
	       //double k_s0 = 1e-5;
	       //double t_0 = 0;
	       //double H_time  = std::exp((t-t_0)*0.008)/(1+std::exp((t-t_0)*0.008));
               //double k_s_time = k_s0 +((k_s-k_s0)*H_time);
               double k_t =0 , k_r =0;
               double dv_t_dc = 0;
               double dv_r_dc = 0;
               //double value = c* (std::exp((c-c_0)*0.008)/(1+std::exp((c-c_0)*0.008)));

          double r = 0;
          if (dim==2)
              r = p.distance(Point<dim>(0.0,0.0));
          else if (dim ==2)
              r = p.distance(Point<dim>(0.0,0.0, 0.0));

               //double R_c = 0.5-t_c;
               double H  = std::exp((r-R_c)*20)/(1+std::exp((r-R_c)*20));  
            
   
                k_t = k_s + (k_s *(b-1)*H);
                k_r = k_s + (k_s *((1/b)-1)*H); 

                v_t = std::pow((1+ (k_t* c)), alpha_g);
                v_r = std::pow((1+ (k_r* c)), alpha_g); 

                dv_t_dc = alpha_g* std::pow((1+ (k_t* c)),(alpha_g -1))* k_t;
		dv_r_dc = alpha_g* std::pow((1+ (k_r* c)),(alpha_g -1))* k_r;

           	F_g = v_t * tangent_direction + v_r * normal_direction;
                G   =  dv_t_dc * tangent_direction + dv_r_dc * normal_direction;

               //inv_F_g = 1/v_t*Physics::Elasticity::StandardTensors< dim >::I +((v_t-v_r)/(v_t*v_r))*normal_direction;
               //  inv_F_g= std::cbrt(1/v_t)*Physics::Elasticity::StandardTensors< dim >::I;
           }

      
      double get_growth_factor_tangent() const {return v_t;}
      double get_growth_factor_radius() const {return v_r;}

     Tensor<2, dim> get_growth_tensor() const { return F_g;}
     Tensor<2, dim> get_growt_tensor_invert() const { return inv_F_g;}
     Tensor<2, dim> get_G() const {return G;}
     


   private:
	
	 const double k_s;
         const double b;
         const double alpha_g;
         const double c_0;
         const double alpha;
         const double R_c;
         double v_t;
         double v_r;
         
         Tensor<2 ,dim> G;
         Tensor<2, dim> F_g;
         Tensor<2, dim> inv_F_g;

};


#endif
