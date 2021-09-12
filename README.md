
App to calculate kinetics parameters for isothermal and non-isothermal reactions.

Three options to run a program:
1. 
- `g++ main.cpp --std=c++17 -lzip -Wextra -Wall -Werror -pedantic-errors -o KineticsLEMS`
- `./KineticsLEMS <PATH_TO_ZIP>`

2. Using CMake 
- `mkdir build`
- `cd build`
- `cmake ..`
- `cmake --build .`
- `./KineticsLEMS <PATH_TO_ZIP>`

3. Using Docker
- `docker build . -t kinetics-lens:latest`
- `docker run kinetics-lens:latest <PATH_TO_ZIP>`
