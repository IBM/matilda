from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pybind11 import get_cmake_dir


__version__ = "0.0.1"

_DEBUG = False
_DEBUG_LEVEL = 0

ext_modules = [
    Pybind11Extension(
        "matildacpp",
        [
            "matilda/cpp_src/bindings.cpp",
            "matilda/cpp_src/ModularInt.cpp",
            "matilda/cpp_src/FilteredSimplicialComplex.cpp",
        ],
        define_macros=[("VERSION_INFO", __version__)],
    ),
]

install_requires = [
    "rectangle-packer",
    "numpy",
    "networkx",
    "matplotlib",
    "tqdm",
    "selenium",
    "scipy"
]

setup(
    name="matilda",
    version=__version__,
    author="",
    author_email="",
    url="",
    packages=["matilda"],
    description="installer test",
    long_description="",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
)
