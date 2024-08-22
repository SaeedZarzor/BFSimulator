# BFSimulator
This project consists of two main components. The first and primary component is implemented in C++ using the deal.II libraries. It includes the computational model required to run simulations and generate results. The second component is implemented in Python3, providing a modern and user-friendly interface for adjusting model parameters and visualizing the results.

![Alt](/Images/example.png)

**Please note:** The current version of this project is optimized for macOS, compatible with both Intel and Apple Silicon processors. If you wish to run it on other platforms (Linux or Windows), please refer to the corresponding branches of the project.

Below, we provide a step-by-step guide for setting up, installing, and configuring the project. If you do not wish to run the second component (the user interface), you may skip to the section "Configuration without Python3."

**Please note:** This project is based on the model introduced in the paper ["Exploring the role of the outer subventricular zone during cortical folding through a physics-based model"](https://elifesciences.org/articles/82925)

## Setup and installation
The following packages and libraries must be installed before running the project:
* **xcode**

First, install Xcode from the App Store. You may also need to install the command line tools by running the following command in the terminal:

```
xcode-select --install
```

You can verify the installation by running ``` xcode-select --version ```, which should return ```xcode-select version 2396```. If not, the installation was not successful.

* **cmake**

To install cmake, first install [Homebrew](https://brew.sh/). Open a new terminal and run:

````
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
````

Then install cmake:
````
brew install cmake 
````
Verify the installation by running ``` cmake --version ```, which should return ```cmake version 3.24.1```. If not, the installation was not successful.

* **compiler and MPI**
Install these with Homebrew:

````
brew install cmake open-mpi gcc@11
````

* **Python3**

Use Homebrew to install Python3:

````
brew install python@3.10
````

Verify the installation by running ```which Python3```, which should return ```/opt/homebrew/bin/Python3```

Next, set up the Tkinter package:

````
brew install python-tk@3.10
````
Finally, ensure PIP is installed for managing Python3 packages. Check with ``` pip --version```. If PIP is not installed, follow the instructions [here](https://www.groovypost.com/howto/install-pip-on-a-mac/#:~:text=To%20install%20PIP%20using%20ensurepip,instructions%20to%20complete%20this%20process.).

* **deal.II**

Install and set up the [deal.II](https://www.dealii.org/) library as follows:

1. Open the terminal and type ```clang``` to trigger the installation of the command line tools.
2. Download the deal.II library with:
````
git clone https://github.com/dealii/candi.git
````
3. Navigate to the downloaded folder: ```cd candi```.

4. Set the environment variables:
````
export OMPI_FC=gfortran-11;export OMPI_CC=clang;export OMPI_CXX=clang++
````
5. Install deal.II:

````
 ./candi.sh --packages="dealii"
````
6. Follow the on-screen instructions (you can abort the process by pressing <CTRL> + C).

If you have trouble installing deal.II please see either [here](https://github.com/dealii/candi) or [here](https://github.com/dealii/dealii/wiki/MacOSX).

Note: We recommend using deal.II version 9.4 to avoid compatibility issues with other versions.

* **Paraview**

Download [Paraview](https://www.paraview.org/) and copy it to the Applications folder.

## Configuration

Next, download and configure the BFSimulator project. Navigate to your desired download directory and run:
````
git clone https://github.com/SaeedZarzor/BFSimulator.git
````

Enter the project folder:
 ```cd BFSimulator ```
 
Then, install the necessary Python3 packages:
````
pip install -r requirements.txt
````

To make the Python files executable, first find the Python3 path with:
 ```which python3``` 
Copy the path and paste it in the first line of ```BFSimulator.py ``` after ```#!```.
Check file permissions with:

```ls -lh BFSimulator.py```

If the response shows ```-rwx------@ ```, the file is executable. Otherwise, make it executable with ```chmod 700 BFSimulator.py```.

Repeat this step for ``` save.py```, ```make_run.py ```, and ```progress.py```


To run the project, use:
````
./BFSimulator.py 
````
## Configuration without Python3
If you wish to run only the first part of the project, follow these steps:

Ensure Xcode with command line tools, CMake, compiler and MPI, and deal.II are installed (see above). Paraview is also recommended for result visualization.

Download the project:

````
git clone https://github.com/SaeedZarzor/BFSimulator.git
````

Navigate to the folder ```cd BFSimulator ``` Generate the Makefile with:
````
cmake CMakeLists.txt
````

Then run: ``` make ```.

To start the simulation:
`````
./Brain_growth Parameters.prm 2 
`````

For 3D simulations, replace 2 with 3. To modify simulation parameters, edit them directly in the ``` Parameters.prm ``` file.
