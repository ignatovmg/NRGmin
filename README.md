# NRGmin 
Utility for molecular mechanics energy minimization

### Installation ###

1. Install [libmol2](https://bitbucket.org/bu-structure/libmol2/src/master)

2. Build minimization executable 

        mkdir build && cd build
        cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS ..
        make && make test
       
3. The binary is in `build/minimization`