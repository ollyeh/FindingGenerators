# Installation
##  Clone the repository:

git clone https://github.com/ollyeh/FindingGenerators.git

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

SKBUILD_SKIP_BUILD=1 -v -e pip install . --config-settings cmake.args="-DGRAPHICS=ON/OFF"