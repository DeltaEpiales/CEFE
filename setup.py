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
            f"-DCMAKE_BUILD_TYPE=Release"
        ]

        build_args = []
        if sys.platform.startswith("win"):
            # Sanitize PATH to remove MSYS2/devkitpro which hijacks CMake
            env_path = os.environ.get("PATH", "")
            clean_paths = [p for p in env_path.split(os.pathsep) if "msys" not in p.lower() and "devkitpro" not in p.lower()]
            os.environ["PATH"] = os.pathsep.join(clean_paths)
            
            cmake_args += ["-G", "Visual Studio 17 2022", "-A", "x64"]
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
    version="0.1.0",
    author="Ryan Osacra",
    description="Causal Entropic Field Engine",
    ext_modules=[CMakeExtension("cefe_py.cefe_core")],
    cmdclass={"build_ext": CMakeBuild},
    packages=["cefe_py"],
    zip_safe=False,
    python_requires=">=3.7",
)
