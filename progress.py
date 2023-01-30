#!/opt/homebrew/bin/Python3
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
    b1 = customtkinter.CTkLabel(root2, text="In Progress", font=("Roboto", 16, "bold"))
    width, height = img.size;
    new_height  = 250
    new_width = new_height * width / height
    b1.place (x=125, y=125, anchor=tk.CENTER)
    for img in ImageSequence.Iterator(img):
        img = customtkinter.CTkImage(light_image=img, dark_image=img, size=(new_width,new_height))
        b1.configure(image = img)
        root2.update()
        time.sleep(0.03)
    
    root2.after (0,play_gif)

def Stop_function():
    for process in psutil.process_iter():
        if process.name() == "Brain_growth":
            os.system(" kill  /im  " + str(process.pid))
    root2.destroy()


customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("dark-blue")

root2 = customtkinter.CTk()
root2.title("Progress")
root2.geometry("250x300")

Stop = customtkinter.CTkButton(master=root2, text="Stop", fg_color="red2", hover_color="darkred", command = Stop_function)
Stop.place(x=125, y=280, anchor=tk.CENTER)




play_gif()

root2.mainloop()
