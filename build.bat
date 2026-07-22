@echo off
echo Cleaning old build files...
rmdir /s /q build
rmdir /s /q cefe.egg-info
del /q cefe_py\*.pyd

echo Building and installing C++ extension in-place...
python setup.py build_ext --inplace

echo Running tests...
python -m pytest tests/

echo Running simulation...
python test_simulation.py

echo Done!
