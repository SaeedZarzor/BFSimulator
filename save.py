#!/opt/homebrew/bin/Python3
# change the directory above #
######################################

from tkinter import  messagebox
from tkinter import filedialog
import shutil
import os
from os.path import exists

my_directroy = filedialog.askdirectory()
directory_folder = "Folder_Output"
parent_dir_folder = os.getcwd()
path_folder = os.path.join(parent_dir_folder, directory_folder)

if my_directroy:
    target_path = os.path.join(my_directroy, directory_folder)
    os.mkdir(target_path)
    for dirs, subdirs, files in os.walk(path_folder):
        for file in files:
            filename = os.path.join(parent_dir_folder, dirs, file)
            shutil.move(filename, target_path)

shutil.rmtree(path_folder)
