#pragma once
#include <iostream>
#include "Nucleon.hh"

class VCollisionReader {
public:
    virtual NucleonVector GetNucleons() = 0;
};
