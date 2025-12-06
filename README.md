# Installation
##  Clone the repository:

git clone -b release https://github.com/ollyeh/FindingGenerators.git

## Go inside of the cloned repository, create a python environment and activate it

cd FindingGenerators
python3.13 -m venv .venv
source .venv/bin/activate

## Use pip to install the package

### For users:
With graphics:

pip install -v . --config-settings cmake.args="-DGRAPHICS=ON"

Without graphics:

pip install -v . --config-settings cmake.args="-DGRAPHICS=OFF"


### For developers:

pip install -v -e . --config-settings cmake.args="-DGRAPHICS=ON"

if wanted: tweak compile time variables in finding_generators/cpp_src/bindings.cpp

then build from finding_generators/CMakeLists.txt

# Documentation
see testing/full_run.py


