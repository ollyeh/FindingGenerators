# Installation
git clone https://github.com/ollyeh/FindingGenerators.git
cd FindingGenerators

python3.13 -m venv .venv
source .venv/bin/activate

For regular usage:
pip install . --config-settings cmake.args="-DGRAPHICS=ON/OFF"

For development:
pip install -v -e . --config-settings cmake.args="-DGRAPHICS=ON/OFF"
