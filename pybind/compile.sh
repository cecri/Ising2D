g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` Ising2D.cpp -o Ising2D`python3-config --extension-suffix` -L$HOME/local/include
