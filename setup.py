import os
import subprocess
import sys
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir="."):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE=Release",
            f"-DCMAKE_Fortran_COMPILER=NOTFOUND"
        ]

        build_args = []
        if sys.platform.startswith("win"):
            # Sanitize PATH only for local builds, since GitHub Actions runners are already clean
            if not os.environ.get("GITHUB_ACTIONS"):
                env_path = os.environ.get("PATH", "")
                clean_paths = [p for p in env_path.split(os.pathsep) if "msys" not in p.lower() and "devkitpro" not in p.lower()]
                os.environ["PATH"] = os.pathsep.join(clean_paths)
                cmake_args += ["-G", "Visual Studio 17 2022", "-A", "x64"]
            else:
                cmake_args += ["-G", "Ninja"]
            
            build_args += ["--config", "Release"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", "-S", ".", "-B", self.build_temp] + cmake_args
        )
        subprocess.check_call(
            ["cmake", "--build", self.build_temp] + build_args
        )

setup(
    name="cefe",
    version="0.2.0",
    author="Ryan Osacra",
    author_email="",
    description="Causal Entropic Field Engine -- QFT/QM simulation on causal diamond lattices",
    long_description=open("README.md", encoding="utf-8").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/DeltaEpiales/CEFE",
    ext_modules=[CMakeExtension("cefe_py.cefe_core")],
    cmdclass={"build_ext": CMakeBuild},
    packages=["cefe_py", "cefe_lab"],
    install_requires=[
        "numpy>=1.20",
        "matplotlib>=3.5",
    ],
    extras_require={
        "gui": [
            "PyQt6>=6.0",
            "pyqtgraph>=0.13",
            "PyOpenGL>=3.1",
        ],
        "dev": [
            "pytest>=7.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "cefe-gui=cefe_py.gui:main",
            "cefe-cli=cefe_py.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False,
    python_requires=">=3.8",
)

