from importlib.resources import files
from .finding_generators import pipe_resource_path

def get_easy3d_install_path():

    return print(str(files("FindingGenerators") / "easy3d" / "resources"))

get_easy3d_install_path()
