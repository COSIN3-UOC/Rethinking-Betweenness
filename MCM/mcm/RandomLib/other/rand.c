// ---------------------------------------------------------------------
// Example program that uses the RNG class
// To compile this example:
// 1.  Download this file along with rng.C and rng.h
// 2.  Compile rng.C and randtest.C.  For example, with g++, type:
//     g++ randtest.C rng.C -o randtest
// 3.  Run the program "randtest"

#include <iostream>
#include "rng.h"
#include "rng.c"

using namespace std;

int main()
{

 RNG x;  // Not seeded explicitly, so it will be given a random seed.

 for (int i = 0; i < 10000; ++i) {
    cout << x.uniform(0,1) << endl;
 }


  return 0;
}

