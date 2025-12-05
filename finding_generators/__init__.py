from importlib.resources import files
from .finding_generators_binding import set_resource_path

from .main import AtomPermutator, GeneratorFinder

def get_easy3d_install_path():
    path = str(files("finding_generators") / "easy3d" / "resources")
    set_resource_path(path)
    return print(f"From __init__.py: {path}")

get_easy3d_install_path()

__all__ = [
    "AtomPermutator",
    "GeneratorFinder",
]