#pragma once

#include "MersenneTwister.h"

///////////////////////////////////////////////////////////////////////////////////////

class randomGenerator {

private :

    MTRand mtwist;                   ///< Mersenne twister object

public :

    inline double getRand() {
        return mtwist.rand();
    }
    inline int getRandInt(unsigned n) {
        return mtwist.randInt(n);
    }

   randomGenerator(const int aseed) : mtwist(aseed) {
       if(aseed < 0)                // randomize the seed
           mtwist.seed();
   }
    randomGenerator() : mtwist() {}
};
