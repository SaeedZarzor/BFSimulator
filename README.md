# BFSimulator
This project consists of two main components. The first and primary component is implemented in C++ using the deal.II libraries. It includes the computational model required to run simulations and generate results. The second component is a modern desktop interface built in Python 3 with PySide6 (Qt) — an adaptive dashboard for adjusting model parameters, running the solver, and visualizing the results.

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

Use Homebrew to install Python 3:

````
brew install python@3.10
````

Verify the installation by running ```which python3```, which should return a path such as ```/opt/homebrew/bin/python3```.

The graphical interface is built with **PySide6 (Qt)**; its Python dependencies are installed from `requirements.txt` during configuration (see below) — no separate Tkinter installation is required. Ensure PIP is available for managing Python 3 packages by running ``` pip --version```. If PIP is not installed, follow the instructions [here](https://www.groovypost.com/howto/install-pip-on-a-mac/#:~:text=To%20install%20PIP%20using%20ensurepip,instructions%20to%20complete%20this%20process.).

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
 
Then install the interface's Python dependencies (mainly PySide6):
````
pip install -r requirements.txt
````

To run the interface, use:
````
python3 BFSimulator.py
````

Alternatively, make it executable once with ```chmod +x BFSimulator.py``` and run ```./BFSimulator.py```. The file uses a portable ```#!/usr/bin/env python3``` shebang, so there is no need to edit any path by hand.

### Using the interface

The interface opens as an adaptive **bento dashboard** with the six parameter categories — Geometry, Advection–Diffusion, Mechanical Properties, Discretization, Numerical Solver, and Growth — laid out as cards that reflow to the window size. Focus any field (or click its category) to see that parameter's figure, symbol, unit, recommended range, and validation in the **Parameter Guide** panel; invalid entries are highlighted with an inline message.

The action bar at the bottom provides:

* **Restore Defaults** — load the recommended 2D or 3D preset.
* **Save / Load Parameters** — store the current configuration to a `.prm` file or reload one.
* **Run Simulation** — build (if needed) and run the C++ solver. A progress window shows the solver's live terminal output and a progress bar based on the simulation time.

When the run finishes, you can browse the results — the folding pattern and the cell-density, stiffness, velocity, growth-factor, and proliferation videos (rendered via ParaView) — and save them to a directory of your choice.

> The Python interface is organized as a small set of modules: `BFSimulator.py` (entry point and parameter window), `bf_fields.py` (parameter registry), `bf_widgets.py` (reusable components), `bf_style.py` (theme), `bf_runner.py` (build/run workers), and `bf_results.py` (results browser).
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
