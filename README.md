# BFSimulator

## Setup and installation
This version of the code works on the macOS powered by either intel or apple silicon. The following packages/libraries are required to install before running the code.

* **xcode**

First, install xcode from the app store. After that, you might need to install the command line tools. To do so you might be to open the terminal and run the following command

```
xcode-select --install
```

You can test that xcode running using the following command ``` xcode-select --version ```. The answering message should be like ```xcode-select version 2396```. Otherwise, the xcode is not installed.

* **cmake**

To install cmake you might need to install [Homebrew](https://brew.sh/) first. Open a native terminal and run the following command 

````
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
````

After that, you could install cmake by typing the following command 
````
brew install cmake 
````

You can test that cmake running using the following command ``` cmake --version ```. The answering message should be like ```cmake version 3.24.1```. Otherwise, the cmake is not installed. 

* **compiler and MPI**

````
brew install cmake open-mpi gcc@11
````

* **Python3**

Again, we will use Homebrew to install Python3. Open terminal and use the command 

````
brew install python@3.10
````

To check the installation success use the command ```which Python3```. The answering message should be like ```/opt/homebrew/bin/Python3```

After that, you might need to setup Tkinter package.  

````
brew install python-tk@3.10
````

Then, you might need to install PIP to easily install and use Python3 packages. First, check if the PIP already installed ``` pip --version```, if not, follow the instructions [here](https://www.groovypost.com/howto/install-pip-on-a-mac/#:~:text=To%20install%20PIP%20using%20ensurepip,instructions%20to%20complete%20this%20process.).

* **deal.II**

The next step is to install and setup [deal.II](https://www.dealii.org/) library. You can do so by following the steps:

1. open terminal and write ```clang``` to trigger the installation of the command line tools. 
2. download the deal.II library using the terminal command,
````
git clone https://github.com/dealii/candi.git
````
3. navigate the terminal to the downloads folder ```cd candi```.

4. run the command
````
export OMPI_FC=gfortran-11;export OMPI_CC=clang;export OMPI_CXX=clang++
````
5. install deal.II by runing the command

````
 ./candi.sh --packages="dealii"
````

If you have trouble installing deal.II please see either [here](https://github.com/dealii/candi) or [here](https://github.com/dealii/dealii/wiki/MacOSX).

P.s. we recommended installing deal.II version 9.4. Otherwise, you could have trouble with later or older versions.

* **Paraview**

Last, you need to download [Paraview](https://www.paraview.org/) and then copy it to the Applications folder.

## Configuration

The following step is to configure the BFSimulator project, but first, we should download it. After navigating the terminal to your to where you want to download the project you could run the command
````
git clone https://github.com/SaeedZarzor/BFSimulator.git
