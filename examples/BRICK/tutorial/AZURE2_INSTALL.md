# AZURE2 Installation

AZURE2 has a number of dependencies, particularly Minuit2, Qt4, and GSL. Below
are a few notes and tips that ought to make the installation a little easier.

Additionally, there is a [video](https://youtu.be/-5h_gjdwhNs) where the BRICK
author, Daniel Odell, and fellow BAND member, Pablo Giuliani, install AZURE2 and
demonstrate BRICK. It is a little bit long, but we hope it is helpful. 

## Minuit2

Minimally, Minuit2 requires:

* git: `sudo apt install git`
* cmake: `sudo apt install cmake`
* g++: `sudo apt install g++`

On Debian 10.10, this looks like
```
git clone https://github.com/GooFit/Minuit2
cd Minuit
mkdir build
cd build
cmake ..
make
sudo make install
```


## Qt4

AZURE2 requires Qt4. It is not compatible with Qt5 or later. (Obtaining Qt4 is
still fairly accessible. However, doing so on Ubuntu for ARM proved to be
difficult.)

On Debian 10.10, the following commands were successful.

* X11: `sudo apt install libx11-dev`
* qt4: `sudo apt install libqt4-dev`

## GSL

Download and unpack the GSL source.
```
wget https://ftp.wayne.edu/gnu/gsl/gsl-2.7.tar.gz
tar -xvf gsl-2.7.tar.gz
```

Compile and install.
```
cd gsl-2.7
./configure && make
sudo make install
```

## And Finally, AZURE2

* There seems to be a C++ standard issue. Editing line 38 of `AZURE2/include/TargetEffect.h` such that `const` is replaced by `constexpr` seems to fix it.
* If you have installed Minuit2 according to the instruction above, you need to specify its path.
```
cmake -DMINUIT_PATH=$HOME/Minuit2/build/lib -DMINUIT2_INCLUDE_DIR=$HOME/Minuit2/build/src ..
make
make install
```
* If you run into errors related to the C++11 standard, try adding `-DCMAKE_CXX_STANDARD=11` to the `cmake` command.
