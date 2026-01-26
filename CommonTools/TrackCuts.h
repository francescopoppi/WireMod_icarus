#pragma once
#include "DataStructures.h"

namespace Cuts {

inline bool trackLength(const Track& t, float minL=50.) {
    return t.length >= minL;
}

inline bool multiplicity1(const Hit& h) {
    return h.mult == 1;
}

inline bool danglingCable(const Hit& h) {
    return (h.tpc==3 && h.y>=80 && h.z>=300 && h.z<=400);
}

}
