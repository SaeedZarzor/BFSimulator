#!/usr/bin/python3
# change the directory above #
######################################

import subprocess
import fileinput
import sys
import os
import tkinter as tk
import customtkinter
from PIL import Image, ImageTk, ImageSequence
from PIL.Image import Resampling
import psutil
import time

def play_gif():
    global img
    img = Image.open('Images/Layer-80.gif')
    b1 = customtkinter.CTkLabel(root3, text="In Progress", font=("Roboto", 16, "bold"))
    width, height = img.size;
    new_height  = 250
    new_width = new_height * width / height
    b1.place (x=125, y=125, anchor=tk.CENTER)
    for img in ImageSequence.Iterator(img):
        img = customtkinter.CTkImage(light_image=img, dark_image=img, size=(new_width,new_height))
        b1.configure(image = img)
        root3.update()
        time.sleep(0.03)
    
    root3.after (0,play_gif)

def Stop_function():
    partial_name = "Brain"
    try:
        # Use pgrep to find the PID of the process by a part of its name
        pid = subprocess.check_output(["pgrep", partial_name]).decode().strip()

        # Use the PID to kill the process
        subprocess.run(["kill", pid])
        print(f"Process containing '{partial_name}' in its name (PID: {pid}) killed successfully.")
        
    except subprocess.CalledProcessError:
        print(f"No process containing '{partial_name}' in its name found.")

        
    root3.destroy()


customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("dark-blue")


root3 = customtkinter.CTk()
root3.title("Progress")
screen_width = root3.winfo_screenwidth()
screen_height = root3.winfo_screenheight()
x_position = ((screen_width - 250) // 4 )
y_position = ((screen_height - 300) // 2)-100

root3.geometry(f"250x300+{x_position}+{y_position}")

Stop = customtkinter.CTkButton(master=root3, text="Stop", fg_color="red2", hover_color="darkred", command = Stop_function)
Stop.place(x=125, y=280, anchor=tk.CENTER)




play_gif()

root3.mainloop()
