import subprocess
import fileinput
import sys
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import  messagebox
import customtkinter
from PIL import Image
from PIL.Image import Resampling
import time
import psutil
import webbrowser


#================================================== the main windows =================================================
Image.MAX_IMAGE_PIXELS = 1000000000
customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("dark-blue")

root = customtkinter.CTk()
root.title("Brain model parameters")
root.geometry("1430x730")

#=============================== define the variables ================================================================
OSVZ_varying_options = ["Constant","Linear-gradient","Quadratic-gradient","Random1","Random2"]
def_vz_raduis = tk.StringVar(root, "")
def_svz_raduis = tk.StringVar(root, "")
def_cr_thickness = tk.StringVar(root, "")
def_intial_raduis = tk.StringVar(root, "")
def_MST_factor = tk.StringVar(root, "")
def_ridial_rate = tk.StringVar(root, "")
def_outer_ridial_rate = tk.StringVar(root, "")
def_intial_dvision = tk.StringVar(root, "")
def_migration_speed = tk.StringVar(root, "")
def_diffusivity = tk.StringVar(root, "")
def_migration_threshold = tk.StringVar(root, "")
def_HV_exp = tk.StringVar(root, "")
def_shear_modulus = tk.StringVar(root, "")
def_stiffness_ratio = tk.StringVar(root, "")
def_poisson_ratio = tk.StringVar(root, "")
def_max_density = tk.StringVar(root, "")
def_stiffness_case = tk.StringVar(root, "")
def_case = tk.StringVar(root, "")
def_refinement = tk.StringVar(root, "")
def_degree = tk.StringVar(root, "")
def_stability_con = tk.StringVar(root, "")
def_c_k = tk.StringVar(root, "")
def_total_time = tk.StringVar(root, "")
def_delt_t = tk.StringVar(root, "")
def_nonlinear_it = tk.StringVar(root, "")
def_tol_u = tk.StringVar(root, "")
def_tol_c = tk.StringVar(root, "")
def_update_u = tk.StringVar(root, "")
def_solver_type = tk.StringVar(root, "")
def_linear_it = tk.StringVar(root, "")
def_k_growth = tk.StringVar(root, "")
def_growth_exp = tk.StringVar(root, "")
def_growth_ratio = tk.StringVar(root, "")
def_ORG_variation_case =  tk.StringVar(root, "   ")

Error_1 = BooleanVar(root, value=False)
Error_2 = BooleanVar(root, value=False)
Error_3 = BooleanVar(root, value=False)
Error_4 = BooleanVar(root, value=False)
Error_5 = BooleanVar(root, value=False)
Error_6 = BooleanVar(root, value=False)
Error_7 = BooleanVar(root, value=False)
Error_8 = BooleanVar(root, value=False)
Error_9 = BooleanVar(root, value=False)
Error_10 = BooleanVar(root, value=False)
Error_11 = BooleanVar(root, value=False)
Error_12 = BooleanVar(root, value=False)
Error_13 = BooleanVar(root, value=False)
Error_14 = BooleanVar(root, value=False)

#================================= Open images ===========================================================

light_backround = Image.open('Images/Untitled-3.png')
dark_backround = Image.open('Images/Untitled-3-dark.png')
light_gemetry_vz=Image.open('Images/gemotry_vz_light.png')
dark_gemetry_vz=Image.open('Images/gemotry_vz_dark.png')
light_gemetry_isvz=Image.open('Images/gemotry_isvz_light.png')
dark_gemetry_isvz=Image.open('Images/gemotry_isvz_dark.png')
light_gemetry_tc=Image.open('Images/gemotry_tc_light.png')
dark_gemetry_tc=Image.open('Images/gemotry_tc_dark.png')
light_gemetry_R=Image.open('Images/gemotry_R_light.png')
dark_gemetry_R=Image.open('Images/gemotry_R_dark.png')
light_gemetry_mst=Image.open('Images/gemotry_mst_light.png')
dark_gemetry_mst=Image.open('Images/gemotry_mst_dark.png')
strain_image_mu_light = Image.open('Images/strain_Energy_mu_light.png')
strain_image_mu_dark = Image.open('Images/strain_Energy_mu_dark.png')
strain_image_nu_light = Image.open('Images/strain_Energy_nu_light.png')
strain_image_nu_dark = Image.open('Images/strain_Energy_nu_dark.png')
strain_image_ratio_light = Image.open('Images/strain_Energy_ratio_light.png')
strain_image_ratio_dark = Image.open('Images/strain_Energy_ratio_dark.png')
growth_image_ratio_light = Image.open('Images/growth_eqautions_ratio_light.png')
growth_image_ratio_dark = Image.open('Images/growth_eqautions_ratio_dark.png')
growth_image_ks_light = Image.open('Images/growth_eqautions_ks_light.png')
growth_image_ks_dark = Image.open('Images/growth_eqautions_ks_dark.png')
growth_image_exp_light = Image.open('Images/growth_eqautions_exp_light.png')
growth_image_exp_dark = Image.open('Images/growth_eqautions_exp_dark.png')
diffusion_image_gamma_light = Image.open('Images/adv_dif_eq_gamma_light.png')
diffusion_image_gamma_dark = Image.open('Images/adv_dif_eq_gamma_dark.png')
diffusion_image_c0_light = Image.open('Images/adv_dif_eq_c0_light.png')
diffusion_image_c0_dark = Image.open('Images/adv_dif_eq_c0_dark.png')
diffusion_image_v_light = Image.open('Images/adv_dif_eq_v_light.png')
diffusion_image_v_dark = Image.open('Images/adv_dif_eq_v_dark.png')
diffusion_image_d_light = Image.open('Images/adv_dif_eq_d_light.png')
diffusion_image_d_dark = Image.open('Images/adv_dif_eq_d_dark.png')
diffusion_image_Gvz_light = Image.open('Images/adv_dif_eq_G_vz_light.png')
diffusion_image_Gvz_dark = Image.open('Images/adv_dif_eq_G_vz_dark.png')
diffusion_image_Gosvz_light = Image.open('Images/adv_dif_eq_G_osvz_light.png')
diffusion_image_Gosvz_dark = Image.open('Images/adv_dif_eq_G_osvz_dark.png')
non_image_light = Image.open('Images/non.png')
non_image_light_2 = Image.open('Images/non_2.png')
two_d_light = Image.open('Images/2D_light.png')
two_d_dark = Image.open('Images/2D_dark.png')
three_d_light = Image.open('Images/3D_light.png')
three_d_dark = Image.open('Images/3D_dark.png')
ref_2d_light =  Image.open('Images/ref_2d_light.png')
ref_2d_dark =  Image.open('Images/ref_2d_dark.png')
varying_light = Image.open('Images/varying.png')
varying_dark = Image.open('Images/varying_dark.png')
varying_cmax_light = Image.open('Images/varying_cmax_light.png')
varying_cmax_dark = Image.open('Images/varying_cmax_dark.png')
logo =  Image.open('Images/Logo_BRAINIACS.png')
author =  Image.open('Images/my_photo.png')
OSVZ_constant = Image.open('Images/OSVZ_Constant.png')
OSVZ_Linear_gradient = Image.open('Images/OSVZ_Linear_gradient.png')
OSVZ_Quadratic_gradient = Image.open('Images/OSVZ_Quadratic_gradient.png')
OSVZ_Random1 = Image.open('Images/OSVZ_Random1.png')
OSVZ_Random2 = Image.open('Images/OSVZ_Random2.png')
#OSVZ_Random3 = Image.open('Images/OSVZ_Random3.png')
OSVZ_Linear_gradient_curve = Image.open('Images/Linear_gradient_curve.png')
OSVZ_Linear_gradient_curve_dark = Image.open('Images/Linear_gradient_curve_dark.png')
OSVZ_Quadratic_gradient_curve = Image.open('Images/Quadratic_gradient_curve.png')
OSVZ_Quadratic_gradient_curve_dark = Image.open('Images/Quadratic_gradient_curve_dark.png')
OSVZ_Random1_curve_dark = Image.open('Images/OSVZ_Random1_curve_dark.png')
OSVZ_Random1_curve = Image.open('Images/OSVZ_Random1_curve.png')
OSVZ_Random2_curve_dark = Image.open('Images/OSVZ_Random2_curve_dark.png')
OSVZ_Random2_curve = Image.open('Images/OSVZ_Random2_curve.png')



width, height = light_backround.size;
new_height  = 200
new_width = new_height * width / height

logo_width, logo_height = logo.size;
new_logo_height  = 140
new_logo_width = new_logo_height * logo_width / logo_height

author_width, author_height = author.size;
new_author_height  = 200
new_author_width = new_author_height * author_width / author_height

#=========================================== the Frames =========================================================

frame_backround = customtkinter.CTkFrame(master=root,
                               width=new_width+10,
                               height=new_height+10,
                               corner_radius=10)
                               
frame_backround.place(x=10, y=10, anchor='nw')

frame_Gemotry = customtkinter.CTkFrame(master=root,
                               width=new_width+10,
                               height=210,
                               corner_radius=10)
                               
frame_Gemotry.place(x=new_width+30, y=10, anchor='nw')

frame_diffusion = customtkinter.CTkFrame(master=root,
                               width=new_width+10,
                               height=330,
                               corner_radius=10)
                               
frame_diffusion.place(x=10, y=new_height+30, anchor='nw')

frame_stiffness = customtkinter.CTkFrame(master=root,
                               width=new_width+10,
                               height=210,
                               corner_radius=10)

frame_stiffness.place(x=2*new_width+50, y=10)


frame_mesh = customtkinter.CTkFrame(master=root,
                               width=new_width+10,
                               height=280,
                               corner_radius=10)
frame_mesh.place(x=new_width+30, y=new_height+30, anchor='nw')

frame_solver = customtkinter.CTkFrame(master=root,
                               width=new_width+10,
                               height=280,
                               corner_radius=10)
                               
frame_solver.place(x =2*new_width+50, y= 210+20 , anchor='nw')

fram_growth = customtkinter.CTkFrame(master=root,
                               width=new_width+10,
                               height=150,
                               corner_radius=10)
                               
fram_growth.place(x =10, y= new_height+330+40 , anchor='nw')

photo_info_fram = customtkinter.CTkFrame(master=root,
                               width=new_width+10+300,
                               height=200,
                               corner_radius=10)
                               
photo_info_fram.place(x=new_width+30, y=new_height+280+40, anchor='nw')


# ================================== instert images =========================================================

backround_image = customtkinter.CTkImage(light_image=light_backround, dark_image=dark_backround, size=(new_width, new_height))
label_img1 = customtkinter.CTkLabel(frame_backround, image = backround_image, text="")
label_img1.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

gemotry_image_vz = customtkinter.CTkImage(light_image=light_gemetry_vz, dark_image=dark_gemetry_vz, size=dark_gemetry_vz.size)
gemotry_image_isvz = customtkinter.CTkImage(light_image=light_gemetry_isvz, dark_image=dark_gemetry_isvz, size=dark_gemetry_isvz.size)
gemotry_image_tc = customtkinter.CTkImage(light_image=light_gemetry_tc, dark_image=dark_gemetry_tc, size=dark_gemetry_tc.size)
gemotry_image_R = customtkinter.CTkImage(light_image=light_gemetry_R, dark_image=dark_gemetry_R, size=dark_gemetry_R.size)
gemotry_image_mst = customtkinter.CTkImage(light_image=light_gemetry_mst, dark_image=dark_gemetry_mst, size=dark_gemetry_mst.size)
strain_image_ratio = customtkinter.CTkImage(light_image=strain_image_ratio_light, dark_image=strain_image_ratio_dark, size=strain_image_ratio_dark.size)
strain_image_mu = customtkinter.CTkImage(light_image=strain_image_mu_light, dark_image=strain_image_mu_dark, size=strain_image_mu_dark.size)
strain_image_nu = customtkinter.CTkImage(light_image=strain_image_nu_light, dark_image=strain_image_nu_dark, size=strain_image_nu_dark.size)
diffusion_image_gamma = customtkinter.CTkImage(light_image=diffusion_image_gamma_light, dark_image=diffusion_image_gamma_dark, size=diffusion_image_gamma_dark.size)
diffusion_image_c0 = customtkinter.CTkImage(light_image=diffusion_image_c0_light, dark_image=diffusion_image_c0_dark, size=diffusion_image_c0_dark.size)
diffusion_image_v = customtkinter.CTkImage(light_image=diffusion_image_v_light, dark_image=diffusion_image_v_dark, size=diffusion_image_v_dark.size)
diffusion_image_d = customtkinter.CTkImage(light_image=diffusion_image_d_light, dark_image=diffusion_image_d_dark, size=diffusion_image_d_dark.size)
diffusion_image_Gvz = customtkinter.CTkImage(light_image=diffusion_image_Gvz_light, dark_image=diffusion_image_Gvz_dark, size=diffusion_image_Gvz_dark.size)
diffusion_image_Gosvz = customtkinter.CTkImage(light_image=diffusion_image_Gosvz_light, dark_image=diffusion_image_Gosvz_dark, size=diffusion_image_Gosvz_dark.size)
growth_image_ks = customtkinter.CTkImage(light_image=growth_image_ks_light, dark_image=growth_image_ks_dark, size=growth_image_ks_light.size)
growth_image_ratio = customtkinter.CTkImage(light_image=growth_image_ratio_light, dark_image=growth_image_ratio_dark, size=growth_image_ratio_light.size)
growth_image_exp = customtkinter.CTkImage(light_image=growth_image_exp_light, dark_image=growth_image_exp_dark, size=growth_image_exp_light.size)
non_image = customtkinter.CTkImage(light_image=non_image_light, dark_image=non_image_light, size=non_image_light.size)
non_image_2 = customtkinter.CTkImage(light_image=non_image_light_2, dark_image=non_image_light_2, size=non_image_light_2.size)
two_d = customtkinter.CTkImage(light_image=two_d_light, dark_image=two_d_dark, size=two_d_light.size)
three_d = customtkinter.CTkImage(light_image=three_d_light, dark_image=three_d_dark, size=three_d_dark.size)
ref_2d = customtkinter.CTkImage(light_image=ref_2d_light, dark_image=ref_2d_dark, size=ref_2d_dark.size)
varying_image = customtkinter.CTkImage(light_image=varying_light, dark_image=varying_dark, size=varying_dark.size)
varying_image_cmax = customtkinter.CTkImage(light_image=varying_cmax_light, dark_image=varying_cmax_dark, size=varying_cmax_dark.size)
logo_image = customtkinter.CTkImage(light_image=logo, dark_image=logo, size=(new_logo_width, new_logo_height))
author_image = customtkinter.CTkImage(light_image=author, dark_image=author, size=(new_author_width, new_author_height))
OSVZ_constant_image = customtkinter.CTkImage(light_image=OSVZ_constant, dark_image=OSVZ_constant, size=OSVZ_constant.size)
OSVZ_Linear_gradient_image = customtkinter.CTkImage(light_image=OSVZ_Linear_gradient, dark_image=OSVZ_Linear_gradient, size=OSVZ_Linear_gradient.size)
OSVZ_Quadratic_gradient_image = customtkinter.CTkImage(light_image=OSVZ_Quadratic_gradient, dark_image=OSVZ_Quadratic_gradient, size=OSVZ_Quadratic_gradient.size)
OSVZ_Random1_image = customtkinter.CTkImage(light_image=OSVZ_Random1, dark_image=OSVZ_Random1, size=OSVZ_Random1.size)
OSVZ_Random2_image = customtkinter.CTkImage(light_image=OSVZ_Random2, dark_image=OSVZ_Random2, size=OSVZ_Random2.size)
#OSVZ_Random3_image = customtkinter.CTkImage(light_image=OSVZ_Random3, dark_image=OSVZ_Random3, size=OSVZ_Random2.size)
OSVZ_Linear_gradient_curve_image = customtkinter.CTkImage(light_image=OSVZ_Linear_gradient_curve, dark_image=OSVZ_Linear_gradient_curve_dark, size=OSVZ_Linear_gradient_curve.size)
OSVZ_Quadratic_gradient_curve_image = customtkinter.CTkImage(light_image=OSVZ_Quadratic_gradient_curve, dark_image=OSVZ_Quadratic_gradient_curve_dark, size=OSVZ_Quadratic_gradient_curve.size)
OSVZ_Random1_curve_image = customtkinter.CTkImage(light_image=OSVZ_Random1_curve, dark_image=OSVZ_Random1_curve_dark, size=OSVZ_Random1_curve.size)
OSVZ_Random2_curve_image = customtkinter.CTkImage(light_image=OSVZ_Random2_curve, dark_image=OSVZ_Random2_curve_dark, size=OSVZ_Random2_curve.size)




label_img2 = customtkinter.CTkLabel(photo_info_fram, text="")
label_img3 = customtkinter.CTkLabel(photo_info_fram, text="")


#================================================ update parameters function =============================================
def update_parameters():
    if  not vz_raduis.get() or not svz_raduis.get() or not cr_thickness.get() or not intial_raduis.get() or not MST_factor.get() or not ridial_rate.get() or not outer_ridial_rate.get() or not intial_dvision.get() or not migration_speed.get() or not diffusivity.get() or not migration_threshold.get() or not HV_exp.get() or not def_stiffness_case.get() or not poisson_ratio.get() or not shear_modulus.get() or not stiffness_ratio.get() or not max_density.get() or not refinement.get() or not degree.get() or not total_time.get() or not delt_t.get() or not stability_con.get() or not k_growth.get() or not growth_ratio.get() or not growth_exp.get() or not c_k.get() or (def_ORG_variation_case.get() == "   "):
        messagebox.showerror("","One or more field is empty!")
        return
  
    if  Error_1.get() or Error_2.get() or Error_3.get() or Error_4.get() or Error_5.get() or Error_6.get() or Error_7.get() or Error_8.get() or Error_9.get() or Error_10.get() or Error_11.get() or Error_12.get() or Error_13.get() or Error_14.get():
        messagebox.showerror("","One or more enrty values not correct!")
        return
 
    for line in fileinput.input("Parameters.prm", inplace=1):
    
        if "set Ventricular zone raduis" in line:
            line = "        set Ventricular zone raduis                       = "+vz_raduis.get()+" \n"
            
        if "set Subventricular zone raduis"in line:
            line = "        set Subventricular zone raduis                    = "+svz_raduis.get()+" \n"

        if "set Cortex thickness" in line:
            line = "        set Cortex thickness                              = "+cr_thickness.get()+" \n"

        if "set Initial radius" in line:
            line = "        set Initial radius                                = "+intial_raduis.get()+" \n"
            
        if "set Mitotic somal translocation factor" in line:
            line = "        set Mitotic somal translocation factor            = "+MST_factor.get()+" \n"
            
        if "set Cell dvision rate of RGCs" in line:
            line = "        set Cell dvision rate of RGCs                    = "+ridial_rate.get()+" \n"

        if "set Cell dvision rate of Outer RGCs" in line:
            line = "        set Cell dvision rate of Outer RGCs              = "+outer_ridial_rate.get()+" \n"

        if "set Cell dvision intial value" in line:
            line = "        set Cell dvision intial value                    = "+intial_dvision.get()+" \n"

        if "set Cell migration speed" in line:
            line = "        set Cell migration speed                         = "+migration_speed.get()+" \n"

        if "set Diffusivity" in line:
            line = "        set Diffusivity                                  = "+diffusivity.get()+" \n"

        if "set Cell migration threshold" in line:
            line = "        set Cell migration threshold                     = "+migration_threshold.get()+" \n"

        if "set Heaviside function exponent" in line:
            line = "        set Heaviside function exponent                  = "+HV_exp.get()+" \n"

        if "set The state of the stiffness" in line:
            line = "        set The state of the stiffness                   = "+def_stiffness_case.get()+" \n"

        if "set Poisson's ratio" in line:
            line = "        set Poisson's ratio                              = "+poisson_ratio.get()+" \n"

        if "set The shear modulus of conrtex" in line:
            line = "        set The shear modulus of conrtex                 = "+shear_modulus.get()+" \n"

        if "set The ratio of stiffness" in line:
            line = "        set The ratio of stiffness                       = "+stiffness_ratio.get()+" \n"

        if "set The max cell density" in line:
            if (def_stiffness_case.get()=='Constant'):
                def_max_density.set(700)
            line = "        set The max cell density                         = "+def_max_density.get()+" \n"

        if "set Number global refinements" in line:
            line = "        set Number global refinements                     = "+refinement.get()+" \n"

        if "set Poly degree" in line:
            line = "        set Poly degree                                   = "+degree.get()+" \n"

        if "set Total time" in line:
            line = "        set Total time                                    = "+total_time.get()+" \n"

        if "set Time step size" in line:
            line = "        set Time step size                                = "+delt_t.get()+" \n"

        if "set Multiplier max iterations linear solver" in line:
            if (def_solver_type.get()=='Direct'):
                def_linear_it.set(100)
            line = "        set Multiplier max iterations linear solver       = "+def_linear_it.get()+" \n"

        if "set Max number newton iterations" in line:
            line = "        set Max number newton iterations                  = "+nonlinear_it.get()+" \n"

        if " set Tolerance residual deformation" in line:
            line = "        set Tolerance residual deformation                =" +tol_u.get()+" \n"

        if "set Tolerance residual diffusion" in line:
            line = "        set Tolerance residual diffusion                  = "+tol_c.get()+" \n"

        if "set Tolerance update" in line:
            line = "        set Tolerance update                              = "+update_u.get()+" \n"

        if "set Growth rate" in line:
            line = "        set Growth rate                                   = "+k_growth.get()+" \n"

        if "set Growth exponent" in line:
            line = "        set Growth exponent                               = "+growth_exp.get()+" \n"

        if "set  Growth ratio" in line:
            line = "        set  Growth ratio                                 = "+growth_ratio.get()+" \n"

        if "set Linear solver type" in line:
            line = "        set Linear solver type                            = "+def_solver_type.get()+" \n"
 
        if "set Stabilization constant" in line:
            line = "        set Stabilization constant                        = "+stability_con.get()+" \n"

        if "set c_k factor" in line:
            line = "        set c_k factor                                   = "+c_k.get()+" \n"
            
        if "set The distribution shape of Outer RGCs" in line:
            line = "        set The OSVZ regional variation                  = "+ def_ORG_variation_case.get()+"\n"
            
        sys.stdout.write(line)

    make_run = subprocess.Popen(['Python3', 'make_run.py', def_case.get()])
    root.destroy()
    
#=============================== set_default_values function ================================================================
def messageWindow():
    win2 = Toplevel(root)
    x_position2 = 565+root.winfo_x()
    y_position2 = 325+root.winfo_y()
    win2.geometry(f"300x80+{x_position2}+{y_position2}")
    win2.title('Choose')
    message2 = "Which case do you want to set?"
    L2 = Label(win2, text=message2)
    L2.place(relx=0.5, rely=0.3, anchor=tk.CENTER)
    B3 = Button(win2, text='2D', command=lambda:[ set_2D_default(), win2.destroy()])
    B3.place(relx=0.3, rely=0.7, anchor=tk.CENTER)
    B4 = Button(win2, text='3D', command=lambda:[ set_3D_default(), win2.destroy()])
    B4.place(relx=0.7, rely=0.7, anchor=tk.CENTER)
    

def set_default_values():
    win = Toplevel(root)
    x_position = 565+root.winfo_x()
    y_position = 325+root.winfo_y()
    win.geometry(f"350x80+{x_position}+{y_position}")
    win.title('Confirmation')
    message = "Are you sure that you want to set all values to  defaults?"
    L = Label(win, text=message)
    L.place(relx=0.5, rely=0.3, anchor=tk.CENTER)
    B1 = Button(win, text='Yes', command=lambda:[win.destroy(), messageWindow()])
    B1.place(relx=0.3, rely=0.7, anchor=tk.CENTER)
    B2 = Button(win, text='No', command = lambda:[win.destroy()])
    B2.place(relx=0.7, rely=0.7, anchor=tk.CENTER)
        
def set_2D_default():
    def_vz_raduis.set(0.25)
    def_svz_raduis.set(0.4)
    def_cr_thickness.set(0.1)
    def_intial_raduis.set(2)
    def_MST_factor.set(0.02)
    def_ridial_rate.set(60)
    def_outer_ridial_rate.set(10)
    def_intial_dvision.set(0)
    def_migration_speed.set(5)
    def_diffusivity.set(0.11)
    def_migration_threshold.set(500)
    def_HV_exp.set(0.008)
    def_shear_modulus.set(2.07)
    def_stiffness_ratio.set(3)
    def_max_density.set(700)
    def_stiffness_case.set('Varying')
    def_poisson_ratio.set(0.38)
    def_case.set('2')
    def_refinement.set(3)
    def_degree.set(2)
    def_total_time.set(1000)
    def_delt_t.set(0.1)
    def_stability_con.set(0.0335)
    def_c_k.set(0.33334)
    def_nonlinear_it.set(8)
    def_tol_c.set(1.0e-8)
    def_tol_u.set(1.0e-8)
    def_update_u.set(1.0e-4)
    def_solver_type.set('Direct')
    def_k_growth.set(4.7e-4)
    def_growth_exp.set(1.65)
    def_growth_ratio.set(1.5)
    def_ORG_variation_case.set(OSVZ_varying_options[0])

    
def set_3D_default():
    def_vz_raduis.set(0.25)
    def_svz_raduis.set(0.4)
    def_cr_thickness.set(0.1)
    def_intial_raduis.set(2)
    def_MST_factor.set(0.02)
    def_ridial_rate.set(600)
    def_outer_ridial_rate.set(100)
    def_intial_dvision.set(0)
    def_migration_speed.set(5)
    def_diffusivity.set(0.11)
    def_migration_threshold.set(500)
    def_HV_exp.set(0.008)
    def_shear_modulus.set(2.07)
    def_stiffness_ratio.set(3)
    def_max_density.set(700)
    def_stiffness_case.set('Varying')
    def_poisson_ratio.set(0.38)
    def_case.set('3')
    def_refinement.set(2)
    def_degree.set(2)
    def_total_time.set(1000)
    def_delt_t.set(0.1)
    def_stability_con.set(0.0335)
    def_c_k.set(0.33334)
    def_nonlinear_it.set(8)
    def_tol_c.set(1.0e-8)
    def_tol_u.set(1.0e-8)
    def_update_u.set(1.0e-4)
    def_solver_type.set('Direct')
    def_k_growth.set(4.7e-5)
    def_growth_exp.set(1.65)
    def_growth_ratio.set(2)
    def_ORG_variation_case.set(OSVZ_varying_options[0])


#======================================== Other functions ======================================

def enable_entry_it():
   linear_it.config(state= "")
   
def disable_entry_it():
   linear_it.config(state= "disabled")
   def_linear_it.set("")
   
def enable_entry_cmax():
   max_density.config(state= "")
   
def disable_entry_cmax():
   max_density.config(state= "disabled")
   def_max_density.set("")
   
def About_author():
    label_img3.configure(image=non_image_2)
    label_img3.place(relx=0.99, rely=0.99, anchor=tk.CENTER)
    label_img2.configure(image=author_image)
    label_img2.place(relx=0.8, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Mohammad Saeed Zarzor \n\n-PhD candidate in the field of Biomechanics.\n-Scientific employee in Institute of Applied Mechanics, \nFriedrich-Alexander-University Erlangen-Nürnberg. \n-Master of Science in Computational Engineering from FAU University. \n-Bachelor of Mechanical Engineering from Damascus University. ")
    info_label.place(relx=0.02, rely=0.07, anchor='nw')
    web_label.configure(text="contact details ")
    web_label.place(relx=0.02, rely=0.7, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://www.ltm.tf.fau.eu/person/zarzor-mohammad-saeed-m-sc/"))
    
def About_programm():
    label_img3.configure(image=non_image_2)
    label_img3.place(relx=0.99, rely=0.99, anchor=tk.CENTER)
    label_img2.configure(image=logo_image)
    label_img2.place(relx=0.85, rely=0.45, anchor=tk.CENTER)
    info_label.configure(text="This work is part of the BRAINIACS Project at the Institute of  Applied \nMechanics, Friedrich-Alexander-University Erlangen-Nürnberg, under \nthe supervision of Dr. Silvia Budday and in cooperation with  Prof. Dr. \nmed Ingmar Blümcke from Neuropathological Institute, University Ho-\nspitals Erlangen. We gratefully acknowledge the funding by the Deut-\nsche Forschungsgemeinschaft (DFG, German Research Foundation). \nThis work based on the following paper: ")
    info_label.place(relx=0.02, rely=0.07, anchor='nw')
    web_label.configure(text="Exploring the role of the outer subventricular zone during cortical fold-\ning through a physics-based model")
    web_label.place(relx=0.02, rely=0.7, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://www.biorxiv.org/content/10.1101/2022.09.25.509401v1.abstract"))
    
def Copy_right():
    label_img3.configure(image=non_image_2)
    label_img3.place(relx=0.99, rely=0.99, anchor=tk.CENTER)
    web_label.configure(text=" CC-BY 4.0 International license.", text_color=('red'))
    web_label.place(relx=0.02, rely=0.3, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://creativecommons.org/licenses/by/4.0/"))
    label_img2.configure(image=non_image)
    label_img2.place(relx=0.5, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="The copyright holder for this preprint is the author/funder, who has granted bioRxiv a license to display \nthe preprint in perpetuity. It is made available under a Copyright:")
    info_label.place(relx=0.02, rely=0.1, anchor='nw')
    


    
def callback(url):
    webbrowser.open_new(url)
#=================================== info functions ==============================================

def VZ_info(event):
    label_img2.configure(image=gemotry_image_vz)
    label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Ventricular zone raduis as a ratio to initial radius \nshould take a value between 0.2 and 0.4.")
    info_label.place(relx=0.4, rely=0.1, anchor='nw')
    
        
def SVZ_info(event):
    label_img2.configure(image=gemotry_image_isvz)
    label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Inner subventricular zone raduis as a ratio to initial radius \nshould take a value between ventricular zone raduis and 0.5.")
    info_label.place(relx=0.4, rely=0.1, anchor='nw')

        
def cr_thickness_info(event):
    label_img2.configure(image=gemotry_image_tc)
    label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Initial cortex thickness as a ratio to initial radius \nshould take a value between 0.01 and 0.35.")
    info_label.place(relx=0.4, rely=0.1, anchor='nw')



def intial_raduis_info(event):
    label_img2.configure(image=gemotry_image_R)
    label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Initial fetal brain radius at gestational week 24 in [mm].")
    info_label.place(relx=0.4, rely=0.1, anchor='nw')

def MST_factor_info(event):
    label_img2.configure(image=gemotry_image_mst)
    label_img2.place(relx=0.29, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Mitotic somal translocation factor of ORGCs \nshould take a value smaller than 0.1.")
    info_label.place(relx=0.5, rely=0.1, anchor='nw')


        
def ridial_rate_info(event):
    label_img2.configure(image=diffusion_image_Gvz)
    label_img2.place(relx=0.3, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Dvision rate in ventricular zone \nin [1/(mm^2 wk)]. \nThis factor mimics RGCs proliferation.")
    info_label.place(relx=0.6, rely=0.1, anchor='nw')
    web_label.configure(text="For more details click here")
    web_label.place(relx=0.6, rely=0.4, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://en.wikipedia.org/wiki/Radial_glial_cell"))

def outer_ridial_rate_info(event):
    label_img2.configure(image=diffusion_image_Gosvz)
    label_img2.place(relx=0.3, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Dvision rate in outer-subventricular zone \nin [1/(mm^2 wk)]. \nThis factor mimics ORGCs proliferation.")
    info_label.place(relx=0.6, rely=0.1, anchor='nw')
    
def ORG_variation_case_info(event):
    if (def_ORG_variation_case.get() == OSVZ_varying_options[0]):
        label_img2.configure(image=OSVZ_constant_image)
        label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
        info_label.configure(text="The variable controllers the regional variation of ORGCs. \nConstant mean no regional variation.")
        info_label.place(relx=0.4, rely=0.1, anchor='nw')
    if (def_ORG_variation_case.get() == OSVZ_varying_options[1]):
        label_img2.configure(image=OSVZ_Linear_gradient_image)
        label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
        info_label.configure(text="The variable controllers the regional variation of ORGCs. \nThe ORGCs division rate increases linearly between \nthe angle 0 and 90.")
        info_label.place(relx=0.4, rely=0.1, anchor='nw')
        label_img3.configure(image=OSVZ_Linear_gradient_curve_image)
        label_img3.place(relx=0.55, rely=0.7, anchor=tk.CENTER)

    if (def_ORG_variation_case.get() == OSVZ_varying_options[2]):
        label_img2.configure(image=OSVZ_Quadratic_gradient_image)
        label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
        info_label.configure(text="The variable controllers the regional variation of ORGCs. \nThe ORGCs division rate increases quadratically \nbetween the angle 0 and 90.")
        info_label.place(relx=0.4, rely=0.1, anchor='nw')
        label_img3.configure(image=OSVZ_Quadratic_gradient_curve_image)
        label_img3.place(relx=0.55, rely=0.7, anchor=tk.CENTER)
        
        
    if (def_ORG_variation_case.get() == OSVZ_varying_options[3]):
        label_img2.configure(image=OSVZ_Random1_image)
        label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
        info_label.configure(text="The variable controllers the regional variation of ORGCs. \nThe random variation of ORGCs division rate as shown in Figure.")
        info_label.place(relx=0.4, rely=0.1, anchor='nw')
        label_img3.configure(image=OSVZ_Random1_curve_image)
        label_img3.place(relx=0.55, rely=0.7, anchor=tk.CENTER)
        
    if (def_ORG_variation_case.get() == OSVZ_varying_options[4]):
        label_img2.configure(image=OSVZ_Random2_image)
        label_img2.place(relx=0.2, rely=0.5, anchor=tk.CENTER)
        info_label.configure(text="The variable controllers the regional variation of ORGCs. \nThe random variation of ORGCs division rate as shown in Figure.")
        info_label.place(relx=0.4, rely=0.1, anchor='nw')
        label_img3.configure(image=OSVZ_Random2_curve_image)
        label_img3.place(relx=0.55, rely=0.7, anchor=tk.CENTER)


def intial_dvision_info(event):
    info_label.configure(text="Cell density intial value in ventricular zone in [1/(mm^2)].")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')


        
def migration_speed_info(event):
    label_img2.configure(image=diffusion_image_v)
    label_img2.place(relx=0.3, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Cell migration speed in [mm/wk]. \nThe cells migrate along RGC fibers \ni.e. radial direction.")
    info_label.place(relx=0.6, rely=0.1, anchor='nw')


def diffusivity_info(event):
    label_img2.configure(image=diffusion_image_d)
    label_img2.place(relx=0.3, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Diffusivity in cortex in [mm^2/wk]. \nIn this model, we consider isotropic diffusion \nwhich  means the diffusion is equal in all \ndirections.")
    info_label.place(relx=0.6, rely=0.1, anchor='nw')


        
def migration_threshold_info(event):
    label_img2.configure(image=diffusion_image_c0)
    label_img2.place(relx=0.3, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Cell migration threshold in [1/mm^3].")
    info_label.place(relx=0.6, rely=0.1, anchor='nw')


def HV_exp_info(event):
    label_img2.configure(image=diffusion_image_gamma)
    label_img2.place(relx=0.3, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Heaviside function exponent. \nFor a smooth solution, should take a value \nsmaller than 0.01.")
    info_label.place(relx=0.6, rely=0.1, anchor='nw')


        
def shear_modulus_info(event):
    label_img2.configure(image=strain_image_mu)
    label_img2.place(relx=0.35, rely=0.4, anchor=tk.CENTER)
    info_label.configure(text="The shear modulus of the cortical layer in [KPa]. The recommended value according to literature is 2.07 KPa.")
    info_label.place(relx=0.05, rely=0.8, anchor='nw')

    
def stiffness_ratio_info(event):
    label_img2.configure(image=strain_image_ratio)
    label_img2.place(relx=0.35, rely=0.4, anchor=tk.CENTER)
    info_label.configure(text="The ratio of stiffness between cortex and subcortex. The recommended values 3 and 5.")
    info_label.place(relx=0.05, rely=0.8, anchor='nw')

    
def poisson_ratio_info(event):
    label_img2.configure(image=strain_image_nu)
    label_img2.place(relx=0.35, rely=0.4, anchor=tk.CENTER)
    info_label.configure(text="The value of the Poisson's ratio, have to take a vlaue between 0.0 and 0.5.")
    info_label.place(relx=0.05, rely=0.75, anchor='nw')
    web_label.configure(text="For more details click here")
    web_label.place(relx=0.05, rely=0.85, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://en.wikipedia.org/wiki/Poisson%27s_ratio"))


    
def max_density_info(event):
    label_img2.configure(image=varying_image_cmax)
    label_img2.place(relx=0.35, rely=0.4, anchor=tk.CENTER)
    info_label.configure(text="The max cell density, after this value the stiffness becomes constant and equals to value set in \n(shear modulus of cortex). Here the c_min was set to 200, so the c_max have to be bigger than 200.")
    info_label.place(relx=0.05, rely=0.8, anchor='nw')

    
def stiffness_case_info(event):
    info_label.configure(text="The state of cortical stiffness is constant.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')

    
def stiffness_case_varying_info(event):
    label_img2.configure(image=varying_image)
    label_img2.place(relx=0.35, rely=0.4, anchor=tk.CENTER)
    info_label.configure(text="The state of cortical stiffness. Varying means exists a positive relation with cell density value.")
    info_label.place(relx=0.05, rely=0.8, anchor='nw')

    
def case_info_2d(event):
    label_img2.configure(image=two_d)
    label_img2.place(relx=0.35, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="")

    
def case_info_3d(event):
    label_img2.configure(image=three_d)
    label_img2.place(relx=0.35, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="")

    
def c_k_info(event):
    info_label.configure(text="c_k factor that satisfies CFL condition, this value matter only in case of Newton Raphson method not \nconverged. In this case, the solver repeat solving the not converged time-step with considering \na smaller time-step size according to the value of c_k.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')
    web_label.configure(text="For more details click here")
    web_label.place(relx=0.05, rely=0.5, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition"))

    
    
def refinement_info(event):
    label_img2.configure(image=ref_2d)
    label_img2.place(relx=0.35, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="The number of mesh global \nrefinements. With increasing this \nvalue the mesh becomes finer \nand the solution smoother but \nwith a longer solving time. \nThe recommended value for the \n2D case is 3 and for the 3D \ncase is 2.")
    info_label.place(relx=0.7 , rely=0.1, anchor='nw')

    
def degree_info(event):
    info_label.configure(text="The shape function polynomial degree of the FE. The shape function is used to approximate the solution. \nwith increasing the degree value the shape functions become smoother and the solution should become \nmore accurate but that increases the solving time.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')
    web_label.configure(text="For more details click here")
    web_label.place(relx=0.05, rely=0.4, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://en.wikipedia.org/wiki/Hp-FEM"))

    
def total_time_info(event):
    info_label.configure(text="Total run time. If you do not know the exact time, write 1000. Thus the solver \nwill automatically stop when it reaches to mechanical instability point.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')

    
def delt_t_info(event):
    info_label.configure(text="Time step size, should take a value smaller than 1.0.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')


def stability_con_info(event):
    info_label.configure(text="Stabilization constant Beta of advection-diffusion equation. \nIn this model, a numerical stabilization method is applied. \nThis value should be smaller than 0.1. The recommended value is 0.03334.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')
    web_label.configure(text="For more details click here")
    web_label.place(relx=0.05, rely=0.4, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://www.dealii.org/current/doxygen/deal.II/step_31.html"))


def nonlinear_it_info(event):
    info_label.configure(text="Max number of nonlinear iterations allowed. \nHere the Newton-Raphson method is used to solve the nonlinear problem.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')
    web_label.configure(text="For more details click here")
    web_label.place(relx=0.05, rely=0.3, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://en.wikipedia.org/wiki/Newton%27s_method"))

    
def tol_u_info(event):
    info_label.configure(text="Force residual error tolerance. The recommended value is 1.0e-8.\nThe smallest allowed value is 1.0e-4.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')

    
def tol_c_info(event):
    info_label.configure(text="Advection-diffusion residual error tolerance. The recommended value is 1.0e-8.\nThe smallest allowed value is 1.0e-4.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')


def update_u_info(event):
    info_label.configure(text="Displacement & cell-density update error tolerance. The recommended value is 1.0e-4.\nThe smallest allowed value is 1.0e-3.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')

    
def solver_type_info(event):
    info_label.configure(text="Type of solver used to solve the linear system. \nIn case of choosing CG solver you should enter the number of linear solver iterations.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')
    web_label.configure(text="For more details click here")
    web_label.place(relx=0.05, rely=0.3, anchor='nw')
    web_label.bind("<Button 1>", lambda e: callback("https://en.wikipedia.org/wiki/Conjugate_gradient_method"))

    
def linear_it_info(event):
    info_label.configure(text="The number of max iterations CG linear solver.")
    info_label.place(relx=0.05, rely=0.1, anchor='nw')

    
def k_growth_info(event):
    label_img2.configure(image=growth_image_ks)
    label_img2.place(relx=0.35, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Growth rate factor. This constant \ncoefficient controls the amount of \nisotropic growth in the subcortical \nlayer. To serve numerical stability \nrequirements, this factor has to \nbe smaller than 1.0e-3 for the 2D \ncase and 1.0e-4 for the 3D case.")
    info_label.place(relx=0.67, rely=0.1, anchor='nw')


def growth_ratio_info(event):
    label_img2.configure(image=growth_image_ratio)
    label_img2.place(relx=0.35, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Growth ratio. This ratio controls \nthe tangential and radial growth \namount in the cortical layer. With \nincreasing this value the growth \nvarying between tangential and \nradial growth increases. \nThe recommended value 1.5 and 3.")
    info_label.place(relx=0.67, rely=0.1, anchor='nw')

    
def growth_exp_info(event):
    label_img2.configure(image=growth_image_exp)
    label_img2.place(relx=0.35, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="Growth exponent.")
    info_label.place(relx=0.67, rely=0.1, anchor='nw')

def info_dis(event):
    label_img3.configure(image=non_image_2)
    label_img3.place(relx=0.99, rely=0.99, anchor=tk.CENTER)
    label_img2.configure(image=non_image)
    label_img2.place(relx=0.5, rely=0.5, anchor=tk.CENTER)
    info_label.configure(text="")
    info_label.place(relx=1, rely=1, anchor=tk.CENTER)
    web_label.configure(text="")
    web_label.place(relx=1, rely=1, anchor=tk.CENTER)

#============================================== Check values function =========================================
def check_vz_raduis(*args):
    if(def_vz_raduis.get() != ''):
        if (float(def_vz_raduis.get()) > 0.4) or (float(def_vz_raduis.get())< 0.2):
            style.configure("vz_style.TEntry",background="red2")
            Error_1.set(value= True)
        else:
            style.configure("vz_style.TEntry",background="systemTextBackgroundColor")
            Error_1.set(value = False)

def check_svz_raduis(*args):
    if(def_svz_raduis.get() != ''):
        if (float(def_svz_raduis.get()) > 0.5) or (float(def_svz_raduis.get())< float(vz_raduis.get())):
            style.configure("svz_style.TEntry",background="red2")
            Error_2.set(value= True)
        else:
            style.configure("svz_style.TEntry",background="systemTextBackgroundColor")
            Error_2.set(value= False)
        
def check_cr_thickness(*args):
    if(def_cr_thickness.get() != ''):
        if (float(def_cr_thickness.get()) < 0.01) or (float(def_cr_thickness.get())> 0.35):
            style.configure("cr_style.TEntry",background="red2")
            Error_3.set(value= True)
        else:
            style.configure("cr_style.TEntry",background="systemTextBackgroundColor")
            Error_3.set(value= False)
        
def check_MST_factor(*args):
    if(def_MST_factor.get() != ''):
        if (float(def_MST_factor.get()) > 0.1):
            style.configure("MST_style.TEntry",background="red2")
            Error_4.set(value= True)
        else:
            style.configure("MST_style.TEntry",background="systemTextBackgroundColor")
            Error_4.set(value= False)

def check_nonlinear_it(*args):
    if(def_nonlinear_it.get() != ''):
        if (float(def_nonlinear_it.get()) < 3):
            style.configure("nonlinear_style.TEntry",background="red2")
            Error_5.set(value= True)
        else:
            style.configure("nonlinear_style.TEntry",background="systemTextBackgroundColor")
            Error_5.set(value= False)
        
def check_tol_u(*args):
    if(def_tol_u.get() != ''):
        if (float(def_tol_u.get()) > 1.0e-4):
            style.configure("tol_u_style.TEntry",background="red2")
            Error_6.set(value= True)
        else:
            style.configure("tol_u_style.TEntry",background="systemTextBackgroundColor")
            Error_6.set(value= False)
        
def check_tol_c(*args):
    if(def_tol_c.get() != ''):
        if (float(def_tol_c.get()) > 1.0e-4):
            style.configure("tol_c_style.TEntry",background="red2")
            Error_7.set(value= True)
        else:
            style.configure("tol_c_style.TEntry",background="systemTextBackgroundColor")
            Error_7.set(value= False)
        
def check_tol_update(*args):
    if(def_update_u.get() != ''):
        if (float(def_update_u.get()) > 1.0e-3):
            style.configure("tol_update_style.TEntry",background="red2")
            Error_8.set(value= True)
        else:
            style.configure("tol_update_style.TEntry",background="systemTextBackgroundColor")
            Error_8.set(value= False)
        
def check_refinement(*args):
    if(def_refinement.get() != ''):
        try:
            int(def_refinement.get())
        except ValueError:
            style.configure("refinementf_style.TEntry",background="red2")
            Error_9.set(value= True)
        else:
            style.configure("refinementf_style.TEntry",background="systemTextBackgroundColor")
            Error_9.set(value= False)
        
def check_degree(*args):
    if(def_degree.get() != ''):
        try:
            integerdegree = int(def_degree.get())
        except ValueError:
            style.configure("degree_style.TEntry",background="red2")
            Error_10.set(value= True)
        else:
            if (integerdegree == 0):
                style.configure("degree_style.TEntry",background="red2")
                Error_10.set(value= True)
            else:
                style.configure("degree_style.TEntry",background="systemTextBackgroundColor")
                Error_10.set(value= False)
        

def check_delt_t(*args):
    if(def_delt_t.get() != ''):
        if (float(def_delt_t.get()) >= 1.0) or (float(def_delt_t.get()) == 0.0):
            style.configure("delt_t_style.TEntry",background="red2")
            Error_11.set(value= True)
        else:
            style.configure("delt_t_style.TEntry",background="systemTextBackgroundColor")
            Error_11.set(value= False)

def check_total_time(*args):
    if(def_total_time.get() != ''):
        try:
            integertime = int(def_total_time.get())
        except ValueError:
            style.configure("total_time_style.TEntry",background="red2")
            Error_12.set(value= True)
        else:
            if (integertime == 0):
                style.configure("total_time_style.TEntry",background="red2")
                Error_12.set(value= True)
            else:
                style.configure("total_time_style.TEntry",background="systemTextBackgroundColor")
                Error_12.set(value= False)
            
def check_poisson_ratio(*args):
    if(def_poisson_ratio.get() != ''):
        if (float(def_poisson_ratio.get()) >0.5) or (float(def_poisson_ratio.get())<0.0):
            style.configure("poisson_ratio_style.TEntry", background = "red2")
            Error_13.set(value= True)
        else:
            style.configure("poisson_ratio_style.TEntry",background="systemTextBackgroundColor")
            Error_13.set(value= False)
            
def cehck_max_density(*args):
    if(def_max_density.get() != ''):
        if (float(def_max_density.get()) <= 200):
            style.configure("max_density_style.TEntry", background = "red2")
            Error_14.set(value= True)
        else:
            style.configure("max_density_style.TEntry",background="systemTextBackgroundColor")
            Error_14.set(value= False)
            
def cehck_stability_con(*args):
    if(def_stability_con.get() != ''):
        if (float(def_stability_con.get()) > 0.1):
            style.configure("stability_con_style.TEntry", background = "red2")
            Error_14.set(value= True)
        else:
            style.configure("stability_con_style.TEntry",background="systemTextBackgroundColor")
            Error_14.set(value= False)
            
def check_c_k(*args):
    if(def_c_k.get() != ''):
        if (float(def_c_k.get()) > 1):
            style.configure("c_k_style.TEntry", background = "red2")
            Error_14.set(value= True)
        else:
            style.configure("c_k_style.TEntry",background="systemTextBackgroundColor")
            Error_14.set(value= False)
            
def check_k_growth(*args):
    if(def_k_growth.get() != ''):
        if (float(def_case.get())==2):
            if (float(def_k_growth.get()) > 1.0e-3):
                style.configure("k_growth_style.TEntry", background = "red2")
                Error_14.set(value= True)
            else:
                style.configure("k_growth_style.TEntry",background="systemTextBackgroundColor")
                Error_14.set(value= False)
        if (float(def_case.get())==3):
            if (float(def_k_growth.get()) > 1.0e-4):
                style.configure("k_growth_style.TEntry", background = "red2")
                Error_14.set(value= True)
            else:
                style.configure("k_growth_style.TEntry",background="systemTextBackgroundColor")
                Error_14.set(value= False)

#============================================== Gemotry frame ============================================
style = ttk.Style()

label_Gemotry = customtkinter.CTkLabel(master=frame_Gemotry, text="Gemotry Parameters", font=("Roboto", 20, "bold"), text_color=("darkgreen"))
label_Gemotry.place(relx=0.5, rely=0.1, anchor=tk.CENTER)

label_vz = customtkinter.CTkLabel(master=frame_Gemotry, text="Ventricular zone raduis:", font=("Roboto", 16)).place(relx=0.05, rely=0.225)

vz_raduis = ttk.Entry(master=frame_Gemotry,font =("Aral",10) ,textvariable = def_vz_raduis , width=15 ,style="vz_style.TEntry")
vz_raduis.place(relx=0.6, rely=0.225)
def_vz_raduis.trace('w', check_vz_raduis)
vz_raduis.bind('<FocusIn>', VZ_info)
vz_raduis.bind('<FocusOut>', info_dis)


label_svz = customtkinter.CTkLabel(master=frame_Gemotry, text="Subventricular zone raduis:", font=("Roboto", 16)).place(relx=0.05, rely=0.375)

svz_raduis = ttk.Entry(master=frame_Gemotry, font =("Aral",10), textvariable = def_svz_raduis ,width=15, style="svz_style.TEntry")
svz_raduis.place(relx=0.6, rely=0.375)
def_svz_raduis.trace('w', check_svz_raduis)
svz_raduis.bind('<FocusIn>', SVZ_info)
svz_raduis.bind('<FocusOut>', info_dis)

label_cth = customtkinter.CTkLabel(master=frame_Gemotry, text="Cortex thickness:", font=("Roboto", 16)).place(relx=0.05, rely=0.525)

cr_thickness = ttk.Entry(master=frame_Gemotry,font =("Aral",10) ,textvariable = def_cr_thickness ,width=15, style = "cr_style.TEntry")
cr_thickness.place(relx=0.6, rely=0.525)
def_cr_thickness.trace('w',check_cr_thickness)
cr_thickness.bind('<FocusIn>', cr_thickness_info)
cr_thickness.bind('<FocusOut>', info_dis)

label_raduis = customtkinter.CTkLabel(master=frame_Gemotry, text="Initial brain radius:", font=("Roboto", 16)).place(relx=0.05, rely=0.675)

intial_raduis = ttk.Entry(master=frame_Gemotry, font =("Aral",10), textvariable = def_intial_raduis ,width=15)
intial_raduis.place(relx=0.6, rely=0.675)
intial_raduis.bind('<FocusIn>', intial_raduis_info)
intial_raduis.bind('<FocusOut>', info_dis)

label_mst = customtkinter.CTkLabel(master=frame_Gemotry, text="Mitotic translocation factor:", font=("Roboto", 16)).place(relx=0.05, rely=0.825)

MST_factor = ttk.Entry(master=frame_Gemotry, font =("Aral",10),textvariable = def_MST_factor  ,width=15, style="MST_style.TEntry")
MST_factor.place(relx=0.6, rely=0.825)
def_MST_factor.trace('w',check_MST_factor)
MST_factor.bind('<FocusIn>', MST_factor_info)
MST_factor.bind('<FocusOut>', info_dis)

#=================================================== Diffusion frame ===========================================

label_diffusion = customtkinter.CTkLabel(master=frame_diffusion, text="Advection diffusion Parameters", font=("Roboto", 20, "bold"), text_color=("darkgreen")).place(relx=0.5, rely=0.08, anchor=tk.CENTER)

label_RG = customtkinter.CTkLabel(master=frame_diffusion, text="Cell dvision rate of RGCs:", font=("Roboto", 16)).place(relx=0.05, rely=0.18)

ridial_rate = ttk.Entry(master=frame_diffusion ,font =("Aral",10) ,textvariable = def_ridial_rate  ,width=15)
ridial_rate.place(relx=0.6, rely=0.18)
ridial_rate.bind('<FocusIn>', ridial_rate_info)
ridial_rate.bind('<FocusOut>', info_dis)

label_ORG = customtkinter.CTkLabel(master=frame_diffusion, text="Cell dvision rate of Outer RGCs:", font=("Roboto", 16)).place(relx=0.05, rely=0.28)

outer_ridial_rate = ttk.Entry(master=frame_diffusion,font =("Aral",10) ,textvariable = def_outer_ridial_rate  , width=15)
outer_ridial_rate.place(relx=0.6, rely=0.28)
outer_ridial_rate.bind('<FocusIn>', outer_ridial_rate_info)
outer_ridial_rate.bind('<FocusOut>', info_dis)

label_ORG_variation = customtkinter.CTkLabel(master=frame_diffusion, text="The OSVZ regional variation:", font=("Roboto", 16)).place(relx=0.05, rely=0.38)

ORG_variation_case = ttk.Combobox(frame_diffusion , textvariable=def_ORG_variation_case ,values = OSVZ_varying_options, width = 11, state='readonly')
ORG_variation_case.place(relx=0.6, rely=0.38)
ORG_variation_case.bind('<FocusIn>', ORG_variation_case_info)
ORG_variation_case.bind('<FocusOut>', info_dis)

label_divsion = customtkinter.CTkLabel(master=frame_diffusion, text="Cell density intial value:", font=("Roboto", 16)).place(relx=0.05, rely=0.48)

intial_dvision = ttk.Entry(master=frame_diffusion, font =("Aral",10) ,textvariable = def_intial_dvision , width=15)
intial_dvision.place(relx=0.6, rely=0.48)
intial_dvision.bind('<FocusIn>', intial_dvision_info)
intial_dvision.bind('<FocusOut>', info_dis)

label_speed = customtkinter.CTkLabel(master=frame_diffusion, text="Cell migration speed:", font=("Roboto", 16)).place(relx=0.05, rely=0.58)

migration_speed = ttk.Entry(master=frame_diffusion,font =("Aral",10), textvariable = def_migration_speed  , width=15)
migration_speed.place(relx=0.6, rely=0.58)
migration_speed.bind('<FocusIn>', migration_speed_info)
migration_speed.bind('<FocusOut>', info_dis)

label_diff = customtkinter.CTkLabel(master=frame_diffusion, text="Diffusivity:", font=("Roboto", 16)).place(relx=0.05, rely=0.68)

diffusivity = ttk.Entry(master=frame_diffusion,font =("Aral",10), textvariable = def_diffusivity , width=15)
diffusivity.place(relx=0.6, rely=0.68)
diffusivity.bind('<FocusIn>', diffusivity_info)
diffusivity.bind('<FocusOut>', info_dis)

label_threshold = customtkinter.CTkLabel(master=frame_diffusion, text="Cell migration threshold:", font=("Roboto", 16)).place(relx=0.05, rely=0.78)

migration_threshold = ttk.Entry(master=frame_diffusion,font =("Aral",10), textvariable = def_migration_threshold  , width=15)
migration_threshold.place(relx=0.6, rely=0.78)
migration_threshold.bind('<FocusIn>', migration_threshold_info)
migration_threshold.bind('<FocusOut>', info_dis)

label_Hvexp = customtkinter.CTkLabel(master=frame_diffusion, text="Heaviside function exponent:", font=("Roboto", 16)).place(relx=0.05, rely=0.88)

HV_exp = ttk.Entry(master=frame_diffusion,font =("Aral",10) , textvariable = def_HV_exp ,width=15)
HV_exp.place(relx=0.6, rely=0.88)
HV_exp.bind('<FocusIn>', HV_exp_info)
HV_exp.bind('<FocusOut>', info_dis)

#=========================================== Stiffness frame ==============================================

label_stiffness = customtkinter.CTkLabel(master=frame_stiffness, text="Mechanical properties Parameters", font=("Roboto", 20, "bold"), text_color=("darkgreen")).place(relx=0.5, rely=0.1, anchor=tk.CENTER)

label_stiffness_case = customtkinter.CTkLabel(master=frame_stiffness, text="The state of the stiffness:", font=("Roboto", 16)).place(relx=0.05, rely=0.25)

stiffness_varying_case = ttk.Radiobutton(master=frame_stiffness ,text= "Varying", variable = def_stiffness_case, value='Varying', command=enable_entry_cmax)
stiffness_constant_case = ttk.Radiobutton(master=frame_stiffness ,text= "Constant", variable = def_stiffness_case, value='Constant', command = disable_entry_cmax)
stiffness_varying_case.place(relx=0.6, rely=0.25)
stiffness_constant_case.place(relx=0.76, rely=0.25)
stiffness_varying_case.bind('<FocusIn>', stiffness_case_varying_info)
stiffness_constant_case.bind('<FocusIn>', stiffness_case_info)
stiffness_varying_case.bind('<FocusOut>', info_dis)
stiffness_constant_case.bind('<FocusOut>', info_dis)



label_shear_modulus = customtkinter.CTkLabel(master=frame_stiffness, text="The shear modulus of conrtex:", font=("Roboto", 16)).place(relx=0.05, rely=0.4)

shear_modulus = ttk.Entry(master=frame_stiffness ,font =("Aral",10) ,textvariable = def_shear_modulus  ,width=15)
shear_modulus.place(relx=0.6, rely=0.4)
shear_modulus.bind('<FocusIn>', shear_modulus_info)
shear_modulus.bind('<FocusOut>', info_dis)

label_stiffness_ratio = customtkinter.CTkLabel(master=frame_stiffness, text="The ratio of stiffness:", font=("Roboto", 16)).place(relx=0.05, rely=0.55)

stiffness_ratio = ttk.Entry(master=frame_stiffness ,font =("Aral",10) ,textvariable = def_stiffness_ratio  ,width=15)
stiffness_ratio.place(relx=0.6, rely=0.55)
stiffness_ratio.bind('<FocusIn>', stiffness_ratio_info)
stiffness_ratio.bind('<FocusOut>', info_dis)

labe_poisson_ratio = customtkinter.CTkLabel(master=frame_stiffness, text="Poisson's ratio:", font=("Roboto", 16)).place(relx=0.05, rely=0.7)

poisson_ratio = ttk.Entry(master=frame_stiffness ,font =("Aral",10) ,textvariable = def_poisson_ratio  ,width=15, style="poisson_ratio_style.TEntry")
poisson_ratio.place(relx=0.6, rely=0.7)
def_poisson_ratio.trace('w', check_poisson_ratio)
poisson_ratio.bind('<FocusIn>', poisson_ratio_info)
poisson_ratio.bind('<FocusOut>', info_dis)

max_density_ratio = customtkinter.CTkLabel(master=frame_stiffness, text="The max cell density:", font=("Roboto", 16)).place(relx=0.05, rely=0.85)

max_density = ttk.Entry(master=frame_stiffness ,font =("Aral",10) ,textvariable = def_max_density ,width=15, style ="max_density_style.TEntry" ,state = "disabled")
max_density.place(relx=0.6, rely=0.85)
def_max_density.trace('w', cehck_max_density)
max_density.bind('<FocusIn>', max_density_info)
max_density.bind('<FocusOut>', info_dis)

#=================================================== mesh frame ===========================================

label_mesh = customtkinter.CTkLabel(master=frame_mesh, text="Discretization Parameters", font=("Roboto", 20, "bold"), text_color=("darkgreen")).place(relx=0.5, rely=0.08, anchor=tk.CENTER)

label_case = customtkinter.CTkLabel(master=frame_mesh, text="State of gemotry:", font=("Roboto", 16)).place(relx=0.05, rely=0.18)

d2_case = ttk.Radiobutton(master=frame_mesh ,text= "2D", variable = def_case, value='2')
d3_case = ttk.Radiobutton(master=frame_mesh ,text= "3D", variable = def_case, value='3')
d2_case.place(relx=0.6, rely=0.18)
d3_case.place(relx=0.75, rely=0.18)
d2_case.bind('<FocusIn>', case_info_2d)
d2_case.bind('<FocusOut>', info_dis)
d3_case.bind('<FocusIn>', case_info_3d)
d3_case.bind('<FocusOut>', info_dis)

label_refinement = customtkinter.CTkLabel(master=frame_mesh, text="Number global refinements:", font=("Roboto", 16)).place(relx=0.05, rely=0.295)

refinement = ttk.Entry(master=frame_mesh, font =("Aral",10) ,textvariable = def_refinement , width=15, style = "refinementf_style.TEntry")
refinement.place(relx=0.6, rely=0.295)
def_refinement.trace('w', check_refinement)
refinement.bind('<FocusIn>', refinement_info)
refinement.bind('<FocusOut>', info_dis)

label_degree = customtkinter.CTkLabel(master=frame_mesh, text="Polynomial degree:", font=("Roboto", 16)).place(relx=0.05, rely=0.41)

degree = ttk.Entry(master=frame_mesh,font =("Aral",10), textvariable = def_degree  , width=15, style = "degree_style.TEntry")
degree.place(relx=0.6, rely=0.41)
def_degree.trace('w', check_degree)
degree.bind('<FocusIn>', degree_info)
degree.bind('<FocusOut>', info_dis)

label_total_time = customtkinter.CTkLabel(master=frame_mesh, text="Total time:", font=("Roboto", 16)).place(relx=0.05, rely=0.525)

total_time = ttk.Entry(master=frame_mesh,font =("Aral",10), textvariable = def_total_time , width=15, style = "total_time_style.TEntry")
total_time.place(relx=0.6, rely=0.525)
def_total_time.trace('w', check_total_time)
total_time.bind('<FocusIn>', total_time_info)
total_time.bind('<FocusOut>', info_dis)

label_delt_t = customtkinter.CTkLabel(master=frame_mesh, text="Time step size:", font=("Roboto", 16)).place(relx=0.05, rely=0.64)

delt_t = ttk.Entry(master=frame_mesh,font =("Aral",10), textvariable = def_delt_t  , width=15, style = "delt_t_style.TEntry")
delt_t.place(relx=0.6, rely=0.64)
def_delt_t.trace('w', check_delt_t)
delt_t.bind('<FocusIn>', delt_t_info)
delt_t.bind('<FocusOut>', info_dis)

label_stability_con = customtkinter.CTkLabel(master=frame_mesh, text="Stabilization constant:", font=("Roboto", 16)).place(relx=0.05, rely=0.755)

stability_con = ttk.Entry(master=frame_mesh,font =("Aral",10) , textvariable = def_stability_con ,width=15, style = "stability_con_style.TEntry")
stability_con.place(relx=0.6, rely=0.755)
def_stability_con.trace('w', cehck_stability_con)
stability_con.bind('<FocusIn>', stability_con_info)
stability_con.bind('<FocusOut>', info_dis)

label_c_k= customtkinter.CTkLabel(master=frame_mesh, text="c_k factor:", font=("Roboto", 16)).place(relx=0.05, rely=0.87)

c_k = ttk.Entry(master=frame_mesh,font =("Aral",10) ,textvariable = def_c_k  , width=15, style="c_k_style.TEntry")
c_k.place(relx=0.6, rely=0.87)
def_c_k.trace('w', check_c_k)
c_k.bind('<FocusIn>', c_k_info)
c_k.bind('<FocusOut>', info_dis)

#============================================= Solver Frame ==================================================


label_solver = customtkinter.CTkLabel(master=frame_solver, text="Numerical solver Parameters", font=("Roboto", 20, "bold"), text_color=("darkgreen")).place(relx=0.5, rely=0.08, anchor=tk.CENTER)


label_newton = customtkinter.CTkLabel(master=frame_solver, text="Max number newton iterations:", font=("Roboto", 16)).place(relx=0.05, rely=0.21)

nonlinear_it = ttk.Entry(master=frame_solver,font =("Aral",10), textvariable = def_nonlinear_it  , width=15, style="nonlinear_style.TEntry")
nonlinear_it.place(relx=0.6, rely=0.21)
def_nonlinear_it.trace('w', check_nonlinear_it)
nonlinear_it.bind('<FocusIn>', nonlinear_it_info)
nonlinear_it.bind('<FocusOut>', info_dis)

label_tol_u = customtkinter.CTkLabel(master=frame_solver, text="Tolerance residual deformation:", font=("Roboto", 16)).place(relx=0.05, rely=0.34)

tol_u = ttk.Entry(master=frame_solver,font =("Aral",10), textvariable = def_tol_u  , width=15, style="tol_u_style.TEntry")
tol_u.place(relx=0.6, rely=0.34)
def_tol_u.trace('w', check_tol_u)
tol_u.bind('<FocusIn>', tol_u_info)
tol_u.bind('<FocusOut>', info_dis)

label_tol_c = customtkinter.CTkLabel(master=frame_solver, text="Tolerance residual diffusion :", font=("Roboto", 16)).place(relx=0.05, rely=0.47)

tol_c = ttk.Entry(master=frame_solver,font =("Aral",10), textvariable = def_tol_c  , width=15, style="tol_c_style.TEntry")
tol_c.place(relx=0.6, rely=0.47)
def_tol_c.trace('w', check_tol_c)
tol_c.bind('<FocusIn>', tol_c_info)
tol_c.bind('<FocusOut>', info_dis)

label_update_u = customtkinter.CTkLabel(master=frame_solver, text="Tolerance update:", font=("Roboto", 16)).place(relx=0.05, rely=0.6)

update_u = ttk.Entry(master=frame_solver,font =("Aral",10), textvariable = def_update_u  , width=15, style="tol_update_style.TEntry")
update_u.place(relx=0.6, rely=0.6)
def_update_u.trace('w', check_tol_update)
update_u.bind('<FocusIn>', update_u_info)
update_u.bind('<FocusOut>', info_dis)

label_solver_type = customtkinter.CTkLabel(master=frame_solver, text="Linear solver type:", font=("Roboto", 16)).place(relx=0.05, rely=0.73)

solver_type_direct= ttk.Radiobutton(master=frame_solver ,text= "Direct", variable = def_solver_type, value='Direct', command=disable_entry_it)
solver_type_cg = ttk.Radiobutton(master=frame_solver ,text= "CG", variable = def_solver_type, value='CG', command=enable_entry_it)
solver_type_direct.place(relx=0.6, rely=0.73)
solver_type_cg.place(relx=0.75, rely=0.73)
solver_type_direct.bind('<FocusIn>', solver_type_info)
solver_type_direct.bind('<FocusOut>', info_dis)
solver_type_cg.bind('<FocusIn>', solver_type_info)
solver_type_cg.bind('<FocusOut>', info_dis)

label_linear_it = customtkinter.CTkLabel(master=frame_solver, text="Iterations linear solver:", font=("Roboto", 16)).place(relx=0.05, rely=0.86)

linear_it = ttk.Entry(master=frame_solver,font =("Aral",10), textvariable = def_linear_it  , width=15, state = "disabled")
linear_it.place(relx=0.6, rely=0.86)
linear_it.bind('<FocusIn>', linear_it_info)
linear_it.bind('<FocusOut>', info_dis)

#=============================================== Growth Frame =================================================
label_growth = customtkinter.CTkLabel(master=fram_growth, text="Growth Parameters", font=("Roboto", 20, "bold"), text_color=("darkgreen")).place(relx=0.5, rely=0.12, anchor=tk.CENTER)

label_k_growth = customtkinter.CTkLabel(master=fram_growth, text="Growth rate:", font=("Roboto", 16)).place(relx=0.05, rely=0.32)

k_growth = ttk.Entry(master=fram_growth,font =("Aral",10), textvariable = def_k_growth  , width=15, style="k_growth_style.TEntry")
k_growth.place(relx=0.6, rely=0.32)
def_k_growth.trace('w', check_k_growth)
k_growth.bind('<FocusIn>', k_growth_info)
k_growth.bind('<FocusOut>', info_dis)

label_growth_ratio = customtkinter.CTkLabel(master=fram_growth, text="Growth ratio:", font=("Roboto", 16)).place(relx=0.05, rely=0.52)

growth_ratio = ttk.Entry(master=fram_growth,font =("Aral",10), textvariable = def_growth_ratio , width=15)
growth_ratio.place(relx=0.6, rely=0.52)
growth_ratio.bind('<FocusIn>', growth_ratio_info)
growth_ratio.bind('<FocusOut>', info_dis)

label_growth_exp = customtkinter.CTkLabel(master=fram_growth, text="Growth exponent:", font=("Roboto", 16)).place(relx=0.05, rely=0.72)

growth_exp = ttk.Entry(master=fram_growth,font =("Aral",10), textvariable = def_growth_exp  , width=15)
growth_exp.place(relx=0.6, rely=0.72)
growth_exp.bind('<FocusIn>', growth_exp_info)
growth_exp.bind('<FocusOut>', info_dis)

#========================================= Buttons==============================================================

run = customtkinter.CTkButton(master=root, text="Run", width = 140 , hover_color="darkgreen"  ,fg_color="green", command=update_parameters )
run.place(x=1335, y=700, anchor=tk.CENTER)

default = customtkinter.CTkButton(master=root, text="Default values",width = 140 ,command=set_default_values )
default.place(x=1335, y=662.5,anchor=tk.CENTER)

About_pro = customtkinter.CTkButton(master=root, text="About programm",width = 140 ,command=About_programm)
About_pro.place(x=1335, y=550,anchor=tk.CENTER)

About_au = customtkinter.CTkButton(master=root, text="About author",width = 140 ,command=About_author)
About_au.place(x=1335, y=587.5,anchor=tk.CENTER)

copyrig = customtkinter.CTkButton(master=root, text="Copyright",width = 140 ,command=Copy_right)
copyrig.place(x=1335, y=625,anchor=tk.CENTER)

info_label = customtkinter.CTkLabel(master=photo_info_fram, text="", font=("Roboto", 14), justify=LEFT)
web_label = customtkinter.CTkLabel(master=photo_info_fram, text="", text_color=('blue'), font=("Roboto", 14), justify=LEFT)


root.mainloop()
