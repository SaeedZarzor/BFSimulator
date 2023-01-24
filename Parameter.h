#ifndef PARAMETER_H
#define PARAMETER_H

#include <deal.II/base/parameter_handler.h>

/*! A namespace to group some functions to declare
 * and parse Parameters into instances of
 * ParameterHandler - READ Tutorial-19
 */



using namespace dealii;

namespace Parameter
{

struct GeneralParameters
{
	
	GeneralParameters(const std::string &input_filename);
	std::string output_file_name;
    std::string solver_type;
    std::string stiffness_case;
    std::string OSVZ_divion_distr;
	double tolerance_residual_u;
    double tolerance_residual_c;
	unsigned int global_refinements;
    double cortex_thickness;
    double ventricular_raduis;
	double subventricular_raduis;
	double Poisson;
	double total_time;
    double time_step;
	unsigned int degree;
    double scale;
    unsigned int max_number_newton_iterations;
    unsigned int multiplier_max_iterations_linear_solver;
    double growth_rate;
    double growth_ratio;
    double initial_radius;
    double cell_dvision_rate_v;
	double cell_dvision_rate_ovz;
	double dvision_value;
    double cell_migration_speed;
    double diffusivity;
    double cell_migration_threshold;
    double exponent;
    double growth_exponent;
    double tol_u;
    double damention_ratio;
    double shear_modulud_cortex;
	double stiffness_ratio;
	double max_cell_density;
	double MST_factor;
    double Betta;
	double c_k;
    
	
	
	void declare_parameters(ParameterHandler &prm);
	void parse_parameters(ParameterHandler &prm);
};

GeneralParameters::GeneralParameters(const std::string &input_filename)
{
	ParameterHandler prm;
	declare_parameters(prm);
	prm.parse_input(input_filename);
	parse_parameters(prm);
}

void GeneralParameters::declare_parameters(ParameterHandler &prm)
{
	prm.enter_subsection ("General");
	{
		prm.declare_entry ("Output file name", "output",
							Patterns::Anything(),
							"The name of the output file to be generated");
		prm.declare_entry ("Poly degree","1",
						   Patterns::Integer(),
				"The polynomial degree of the FE");
		prm.declare_entry ("Number global refinements","1",
						   Patterns::Integer(),
				 "The number of mesh global refinements");
                prm.declare_entry ("Grid scale","1",
						   Patterns::Double(),
				 "Grid scale");
                prm.declare_entry ("Damention ratio","1",
						   Patterns::Double(),
				 "Damention ratio");
                prm.declare_entry ("Initial radius","0.5",
						   Patterns::Double(),
					"Initial radius [mm]");
                prm.declare_entry ("Cortex thickness","0.1",
						   Patterns::Double(),
					"Initial cortex thickness as a ratio to initial radius (value between 0 and 1)");
                prm.declare_entry ("Ventricular zone raduis","0.25",
						   Patterns::Double(),
					"Ventricular zone raduis as a ratio to initial radius (value between 0.2 and 0.6)");

		prm.declare_entry ("Subventricular zone raduis","0.4",
						   Patterns::Double(),
					"Subventricular zone raduis as a ratio to initial radius (value between ventricular zone raduis and 0.6)");

	 	prm.declare_entry ("Mitotic somal translocation factor","0.05",
						   Patterns::Double(),
					"Mitotic somal translocation factor (value <= 0.1)");

		prm.declare_entry ("Total time","1",
						   Patterns::Double(),
						   "Total time");
                prm.declare_entry ("Multiplier max iterations linear solver","10",
						   Patterns::Integer(),
						   "Multiplier max iterations linear solver");
                prm.declare_entry ("Max number newton iterations","10",
						   Patterns::Integer(),
						   "Max number newton iterations");
                prm.declare_entry ("Time step size","0.1",
						   Patterns::Double(),
						   "Time step size");
		prm.declare_entry ("Tolerance residual deformation","1e-5",
						   Patterns::Double(),
						   "The tolerance wrt the normalised residual deformation");
                prm.declare_entry ("Tolerance residual diffusion","1e-6",
						   Patterns::Double(),
						   "The tolerance wrt the normalised residual diffusion");
                prm.declare_entry ("Tolerance update","1e-6",
						   Patterns::Double(),
						   "The tolerance wrt the normalised residual norm");
        prm.declare_entry("The state of the stiffness","Varying",
                          Patterns::Anything(),
                          "The state of the stiffness Constant or Varying");
		prm.declare_entry ("Poisson's ratio","0.45",
						   Patterns::Double(),
						   "The value of the Poisson's ratio");
                prm.declare_entry ("The shear modulus of conrtex","2",
						   Patterns::Double(),
						   "The shear modulus of conrtex layer [KPa]");
		prm.declare_entry ("The ratio of stiffness","3.5",
						   Patterns::Double(),
						   "The ratio of stiffness between cortex and subcortex  mu_cmax/mu_smax");
        
		prm.declare_entry ("The max cell density","400",
						   Patterns::Double(),
						   "The max cell density where the stiffness still constant c_max");
		prm.declare_entry ("Growth rate","1e-3",
						   Patterns::Double(),
					"Growth rate k_s mm^2");
                prm.declare_entry ("Growth exponent","1",
						   Patterns::Double(),
					"Growth exponent alpha ");
                prm.declare_entry ("Growth ratio","0.5",
						   Patterns::Double(),
					"Growth ratio bitta_k");

               prm.declare_entry ("Cell dvision rate of RGCs","120",
						   Patterns::Double(),
					"Cell dvision rate in ventricular zone G^c_v [1/(mm^2 wk)]");
	       prm.declare_entry ("Cell dvision rate of Outer RGCs","0",
						   Patterns::Double(),
					"Cell dvision rate in outer-subventricular zone G^c_ovz [1/(mm^2 wk)]");
        prm.declare_entry ("The OSVZ regional variation","0",
                        Patterns::Anything(),
                    "Constant, Linear-gradient, Random1, Random2, Random3");
	       prm.declare_entry ("Cell dvision intial value","5",
						   Patterns::Double(),
					"Cell dvision intial value [1/(mm^2)]");
	       prm.declare_entry ("Cell migration speed","5",
 						   Patterns::Double(),
					"Cell migration speed  [mm/wk]");
               prm.declare_entry ("Diffusivity","0.1",
 						   Patterns::Double(),
					"Diffusivity  d^cc [mm^2/wk]");
               prm.declare_entry ("Cell migration threshold","400",
 						   Patterns::Double(),
					"Cell migration threshold  c_0 [1/mm^3]");
               prm.declare_entry ("Heaviside function exponent","20",
 						   Patterns::Double(),
					"Heaviside function exponent  gamma");
               prm.declare_entry ("Linear solver type", "CG",
					      Patterns::Anything(),
						"Linear solver type (CG or Direct)");
               prm.declare_entry ("Stabilization constant", "0.017",
                           Patterns::Double(),
                           "Stabilization constant Betta");
               prm.declare_entry ("c_k factor", "0.25",
                           Patterns::Double(),
                           "c_k factor (CFL condition)");

 
	}
	prm.leave_subsection ();
}

void GeneralParameters::parse_parameters (ParameterHandler &prm)
{
	prm.enter_subsection("General");
    {
        stiffness_case = prm.get("The state of the stiffness");
		tolerance_residual_u=prm.get_double("Tolerance residual deformation");
                tolerance_residual_c=prm.get_double("Tolerance residual diffusion");
		global_refinements =prm.get_integer("Number global refinements");
		Poisson=prm.get_double("Poisson's ratio");
		shear_modulud_cortex=prm.get_double("The shear modulus of conrtex");
	        stiffness_ratio=prm.get_double("The ratio of stiffness");
		max_cell_density=prm.get_double("The max cell density");
		total_time = prm.get_double("Total time");
                time_step = prm.get_double("Time step size");
		output_file_name = prm.get("Output file name");
		degree = prm.get_integer("Poly degree");
                scale = prm.get_double("Grid scale");
                max_number_newton_iterations = prm.get_integer("Max number newton iterations");
                multiplier_max_iterations_linear_solver = prm.get_integer("Multiplier max iterations linear solver");
                growth_rate = prm.get_double("Growth rate");
                growth_ratio = prm.get_double("Growth ratio");
                initial_radius = prm.get_double("Initial radius");
                cortex_thickness = prm.get_double("Cortex thickness");
                ventricular_raduis = prm.get_double("Ventricular zone raduis");
		subventricular_raduis = prm.get_double("Subventricular zone raduis");
		MST_factor = prm.get_double("Mitotic somal translocation factor");
                cell_dvision_rate_v = prm.get_double("Cell dvision rate of RGCs");
 		cell_dvision_rate_ovz = prm.get_double("Cell dvision rate of Outer RGCs");
        OSVZ_divion_distr = prm.get("The OSVZ regional variation");
		dvision_value = prm.get_double("Cell dvision intial value");
                cell_migration_speed = prm.get_double("Cell migration speed");
                diffusivity = prm.get_double("Diffusivity");
                cell_migration_threshold = prm.get_double("Cell migration threshold");
                exponent = prm.get_double("Heaviside function exponent");
                growth_exponent = prm.get_double("Growth exponent");
                solver_type = prm.get("Linear solver type");
                tol_u = prm.get_double("Tolerance update");
                damention_ratio = prm.get_double("Damention ratio");
    	        Betta = prm.get_double("Stabilization constant");
		c_k = prm.get_double("c_k factor");


    }
	prm.leave_subsection();
}


}//END namespace
#endif
