# BFSimulator

The project includes two parts. The first and core part is implemented using C++ and the deal.II libraries. It includes the computational model to run the simulations and generate the results. The second part is implemented using Python3 and provides a modern and user-friendly interface to change the model parameters and visualize the results.

![Alt](/Images/example.png)

**Place note:** The current project version works for macOS powered by either intel or apple silicon. In case you would like to run it on other platforms (Linux or Windows), we would like to refer to the other project branches.

In the following, we explain the setup, installation, and configuration of the project step by step. If you do not want to run the second part of the project (user interface), place jump to the "configuration without Python3".

**Place note:** This project is built based on the model introduced in the following paper ["Exploring the role of the outer subventricular zone during cortical folding through a physics-based model"](https://www.biorxiv.org/content/10.1101/2022.09.25.509401v1.abstract)

## Setup and installation
The following packages/libraries are required to be installed before running the project.

* **xcode**

First, install xcode from the app store. After that, you might need to install the command line tools, e.g., via the terminal running the following command:

```
xcode-select --install
```

You can test whether xcode is running using the following command ``` xcode-select --version ```. You should receive the message ```xcode-select version 2396```. Otherwise, xcode was not installed properly.

* **cmake**

To install cmake, you might need to install [Homebrew](https://brew.sh/) first. Open a new terminal and run the following command 

````
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
````

After that, you can install cmake by typing  
````
brew install cmake 
````

You can test whether cmake is running using the following command ``` cmake --version ```.  You should receive the message ```cmake version 3.24.1```. Otherwise, cmake is not installed properly. 

* **compiler and MPI**

````
brew install cmake open-mpi gcc@11
````

* **Python3**

Again, we will use Homebrew to install Python3. Open a terminal and use the command 

````
brew install python@3.10
````

To check the installation success use the command ```which Python3```. You should receive the message ```/opt/homebrew/bin/Python3```

After that, you might need to setup the Tkinter package.  

````
brew install python-tk@3.10
````

Then, you might need to install PIP to easily install and use Python3 packages. First, check if the PIP is already installed using ``` pip --version```. If not, follow the instructions [here](https://www.groovypost.com/howto/install-pip-on-a-mac/#:~:text=To%20install%20PIP%20using%20ensurepip,instructions%20to%20complete%20this%20process.).

* **deal.II**

The next step is to install and setup the [deal.II](https://www.dealii.org/) library by following the steps:

1. open terminal and write ```clang``` to trigger the installation of the command line tools. 
2. download the deal.II library using the terminal command,
````
git clone https://github.com/dealii/candi.git
````
3. navigate the terminal to the download folder ```cd candi```.

4. run the command
````
export OMPI_FC=gfortran-11;export OMPI_CC=clang;export OMPI_CXX=clang++
````
5. install deal.II by running the command

````
 ./candi.sh --packages="dealii"
````
6. Follow the instructions on the screen (you can abort the process by pressing < CTRL > + C)

If you have trouble installing deal.II please see either [here](https://github.com/dealii/candi) or [here](https://github.com/dealii/dealii/wiki/MacOSX).

Note: we recommend installing deal.II version 9.4. Otherwise, you could have trouble with later or older versions.

* **Paraview**

Lastly, you need to download [Paraview](https://www.paraview.org/) and then copy it to the Applications folder.

## Configuration

The next step is to download and configure the BFSimulator project. After navigating the terminal to where you want to download the project, you shall run the command
````
git clone https://github.com/SaeedZarzor/BFSimulator.git
````

Now use ```cd BFSimulator ```to enter the folder. Then you need to install all necessary Python3 packages to run the project using the command
````
pip install -r requirements.txt
````

After that, we want to make the python files executable. In the terminal window, write the command ```which python3``` to get the path to the python3 installed on your machine. Copy the path and open ```BFSimulator.py ```, then paste the path in the first line after ```#!```.
In the terminal use the command ```ls -lh BFSimulator.py```to know the files´s permission. If the response is
````
-rwx------@ 1 saeed  staff    68K Jan 30 10:18 BFSimulator.py
````
it means that the file is already considered executable. Otherwise, run the command ``` chmod 700  BFSimulator.py``` to change the file's permission.

Repeat the previous step again for the files: ``` save.py```, ```make_run.py ```, ```progress.py```


To run the Project use the command
````
./BFSimulator.py 
````
## Configuration without Python3
In this section, we explain how to run the first part of the project independently. If you have already implemented all previous steps and everything works well, you are done. In case you have some trouble with python, however, you might need to proceed with the following section.

First, you might need to install xcode with command line tools, cmake, compiler and MPI, and deal.II (see above). It is also recommended to install Paraview to visualize the results. 

After successful installation, you can download the project by using the following terminal command

````
git clone https://github.com/SaeedZarzor/BFSimulator.git
````

Now use ```cd BFSimulator ```to enter the folder and write the command

````
cmake CMakeLists.txt
````

After cmake generates MakeFile successfully, run the ``` make ``` command.

To start the simulation use the command 

`````
./Brain_growth Parameters.prm 2 
`````

If you want to simulate a 3D case, replace 2 by 3. To change the simulation parameters, open the ``` Parameters.prm ``` file and change them directly.
