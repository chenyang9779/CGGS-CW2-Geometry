# Section 4

To compile the code:
```bash
cd code
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
Alternatively, you can use the `-G` option in `cmake` and create a solution for your chosen C++ IDE (for instance, `-G Xcode`).


To use the code (reproducing results in the report), in the build folder:
```bash
./Section4/Section4 <path/to/mesh.off>
```