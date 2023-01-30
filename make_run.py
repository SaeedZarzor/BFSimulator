#!/opt/homebrew/bin/Python3
# change the directory above #
######################################

import subprocess
import fileinput
import sys
import time
import psutil
import cv2
import shutil
import os
import re
from os.path import exists
import customtkinter
import tkinter as tk
from tkinter import ttk
from tkinter import  messagebox
from tkinter import filedialog
from pathlib import Path

sys.getfilesystemencoding()
#============================================= functions ===========================================
def cell_video():
    path = 'Folder_Output/Cell_desnity.avi'
    file_exists = exists(path)

    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
        
        Code_name = 'cell_density_video_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',  Code_name])

        
    cap=cv2.VideoCapture(path)
    
    if not cap.isOpened():
        print("Error opening Video File.")

    while(True):
        ret, frame=cap.read()
        if not ret:
            break
        cv2.imshow("Cell density", frame)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

    cv2.waitKey(5000)
    cv2.destroyAllWindows()
    cap.release()
    

def folds_pattern_image():
    path = 'Folder_Output/folds_pattern.png'
    file_exists = exists(path)
    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
           
        Code_name = 'final_folding_pattren_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',   Code_name])

        
    image = cv2.imread(path)
    cv2.imshow ("Folds pattern", image)
    k = cv2.waitKey(0)
    if (k == ord('q')):
        cv2.destroyAllWindows()
        
def velocity_video():
    path = 'Folder_Output/Velocity.avi'
    file_exists = exists(path)
    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
                    
        Code_name = 'velocity_video_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',   Code_name])
        
    cap=cv2.VideoCapture(path)

    if not cap.isOpened():
        print("Error opening Video File.")

    while(True):
        ret, frame=cap.read()
        if not ret:
            break
        cv2.imshow("Velocity", frame)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

    cv2.waitKey(5000)
    cv2.destroyAllWindows()
    cap.release()
    
def stiffness_video():
    path = 'Folder_Output/Stiffness.avi'
    file_exists = exists(path)
    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
        
        Code_name = 'Stiffness_video_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',   Code_name])
        
    cap=cv2.VideoCapture(path)

    if not cap.isOpened():
        print("Error opening Video File.")

    while(True):
        ret, frame=cap.read()
        if not ret:
            break
        cv2.imshow("Stiffness", frame)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

    cv2.waitKey(5000)
    cv2.destroyAllWindows()
    cap.release()
    
def growth_factores_t_video():
    path = 'Folder_Output/growth_factor_t.avi'
    file_exists = exists(path)
    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
                    
        Code_name = 'growth_factors_video_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',   Code_name])
#        export_results.wait()
        
    cap=cv2.VideoCapture(path)

    if not cap.isOpened():
        print("Error opening Video File.")

    while(True):
        ret, frame=cap.read()
        if not ret:
            break
        cv2.imshow("Tangential growth factor", frame)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

    cv2.waitKey(5000)
    cv2.destroyAllWindows()
    cap.release()
    
def growth_factores_r_video():
    path = 'Folder_Output/growth_factor_r.avi'
    file_exists = exists(path)
    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
                
        Code_name = 'growth_factors_video_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',   Code_name])
        
    cap=cv2.VideoCapture(path)

    if not cap.isOpened():
        print("Error opening Video File.")

    while(True):
        ret, frame=cap.read()
        if not ret:
            break
        cv2.imshow("Radial growth factor", frame)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

    cv2.waitKey(5000)
    cv2.destroyAllWindows()
    cap.release()
    
def source_vz_video():
    path = 'Folder_Output/source_vz.avi'
    file_exists = exists(path)
    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
                    
        Code_name = 'source_terms_video_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',   Code_name])
        
    cap=cv2.VideoCapture(path)

    if not cap.isOpened():
        print("Error opening Video File.")

    while(True):
        ret, frame=cap.read()
        if not ret:
            break
        cv2.imshow("RGCs Proliferation", frame)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

    cv2.waitKey(5000)
    cv2.destroyAllWindows()
    cap.release()
    
def source_osvz_video():
    path = 'Folder_Output/source_osvz.avi'
    file_exists = exists(path)
    if not file_exists:
        Files = Path('/Applications')
        for file in Files.iterdir():
            if re.match("ParaView", file.name):
                    Para_view = file
                    
        Code_name = 'source_terms_video_'+str(sys.argv[1])+'.py'
        export_results = subprocess.run([str(Para_view)+'/Contents/bin/pvpython',   Code_name])#        export_results.wait()
        
    cap=cv2.VideoCapture(path)

    if not cap.isOpened():
        print("Error opening Video File.")

    while(True):
        ret, frame=cap.read()
        if not ret:
            break
        cv2.imshow("ORGCs Proliferation", frame)
        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

    cv2.waitKey(5000)
    cv2.destroyAllWindows()
    cap.release()
    
def close_fun():
    close_qes = messagebox.askquestion("","Are you sure you want to continue without saving results?")
    if close_qes =="yes":
        root.destroy()
        directory_folder = "Folder_Output"
        parent_dir_folder = os.getcwd()
        path_folder = os.path.join(parent_dir_folder, directory_folder)
        shutil.rmtree(path_folder)
        sys.exit()
    else:
        save_fun(True)
            
            
def save_fun(in_root):
    if in_root:
        root.destroy()
    subprocess.run(['./save.py'])
    sys.exit()

# =============================================================================================
customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("dark-blue")

#current_task_pid_number = 0
#for process in psutil.process_iter():
#    if process.name() == "Python" :
#        current_task_pid_number = process.pid

cmake_file = 'Makefile'
make_file = 'Brain_growth'



if not exists(cmake_file):
    subprocess.run(['cmake','CMakeLists.txt'])
    
if not exists(make_file):
    subprocess.run(['make'])
    
progress = subprocess.Popen('./progress.py')
subprocess.run(['./Brain_growth', 'Parameters.prm',  sys.argv[1]])

#for process in psutil.process_iter():
#    if re.match("progress", process.name()):
#        os.system("kill /im "+str(process.pid))’

directory_folder = "Folder_Output"
parent_dir_folder = os.getcwd()
path_folder = os.path.join(parent_dir_folder, directory_folder)
if not os.path.exists(path_folder):
    os.mkdir(path_folder)
else:
    shutil.rmtree(path_folder)
    os.mkdir(path_folder)
    
for dirs, subdirs, files in os.walk(parent_dir_folder):
    for file in files:
        if file.endswith('.vtk'):
            filename =os.path.join(parent_dir_folder, dirs, file)
            file_exists =os.path.join(path_folder, file)
            if not os.path.exists(file_exists):
                shutil.move(filename, path_folder)
                
csv_filename =os.path.join(parent_dir_folder,"timeing.csv")
shutil.copy(csv_filename, path_folder)
prm_filename =os.path.join(parent_dir_folder,"Parameters.prm")
shutil.copy(prm_filename, path_folder)
                

post_processing = messagebox.askquestion("", "Done! \n Do you want to see the results?")
if post_processing=="yes":
    root = customtkinter.CTk()
    root.title("Results Output")
    root.geometry("300x500")
    folds_pattern = customtkinter.CTkButton(master=root, text="Folds pattern", width = 250, command=folds_pattern_image)
    folds_pattern.pack(ipadx=10, ipady=5, expand=True)
    cell_video_button = customtkinter.CTkButton(master=root, text="Cell density distribution", width = 250, command=cell_video)
    cell_video_button.pack(ipadx=10, ipady=5, expand=True)
    stiffness_video_button = customtkinter.CTkButton(master=root, text="Stiffness distribution", width = 250, command=stiffness_video)
    stiffness_video_button.pack(ipadx=10, ipady=5, expand=True)
    velocity_video_button = customtkinter.CTkButton(master=root, text="Velocity distribution", width = 250, command=velocity_video)
    velocity_video_button.pack(ipadx=10, ipady=5, expand=True)
    growth_t_video_button = customtkinter.CTkButton(master=root, text="Tangential growth factor distribution", width = 250, command=growth_factores_t_video)
    growth_t_video_button.pack(ipadx=10, ipady=5, expand=True)
    growth_r_video_button = customtkinter.CTkButton(master=root, text="Radial growth factor distribution", width = 250, command=growth_factores_r_video)
    growth_r_video_button.pack(ipadx=10, ipady=5, expand=True)
    source_vz_video_button = customtkinter.CTkButton(master=root, text="RGCs proliferation", width = 250, command=source_vz_video)
    source_vz_video_button.pack(ipadx=10, ipady=5, expand=True)
    source_osvz_video_button = customtkinter.CTkButton(master=root, text="ORGCs proliferation", width = 250, command=source_osvz_video)
    source_osvz_video_button.pack(ipadx=10, ipady=5, expand=True)
    close_button = customtkinter.CTkButton(master=root, text="Close", hover_color="darkred"  ,fg_color="red" ,width = 100, command=close_fun)
    close_button.pack(ipadx=10, ipady=5,side=tk.RIGHT,expand=True, pady=6)
    save_button = customtkinter.CTkButton(master=root, text="Save", hover_color="darkgreen"  ,fg_color="green" ,width = 100, command = lambda:save_fun(True))
    save_button.pack(ipadx=10, ipady=5,side=tk.LEFT,expand=True, pady=6)
    root.mainloop()
        
if post_processing=="no":
    qes = messagebox.askquestion("","Do you want to save the results?")
    if qes =="yes":
        save_fun(False)
    else:
        shutil.rmtree(path_folder)
        sys.exit()
                    




