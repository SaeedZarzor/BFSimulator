# BFSimulator

## Setup and installation
This version of the code works on the macOS powered by either intel or apple silicon. The following packages/libraries are required to install before running the code.

* **xcode**

First, install xcode from the app store. After that, you might need to install the command line tools. To do so you might be to open the Terminal and run the following command

```
xcode-select --install
```

You can test that xcode running using the following command ``` xcode-select --version ```. The answering message should be like ```xcode-select version 2396```. Otherwise, the xcode is not installed.

* **cmake**

To install cmake you might need to install [Homebrew](https://brew.sh/) first. Open a native termainl and run the following command 

````
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
````

After that, you could install cmake by typing the following command 
````
brew install cmake 
````

You can test that cmake running using the following command ``` cmake --version ```. The answering message should be like ```cmake version 3.24.1```. Otherwise, the cmake is not installed. 
