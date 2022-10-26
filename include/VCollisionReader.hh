#pragma once
#include <iostream>
#include "Nucleon.hh"
#include "AAMCConstants.hh"

class VCollisionReader {
public:
    VCollisionReader(const VCollisionReader&) = delete;
    VCollisionReader& operator=(const VCollisionReader&) = delete;
    virtual ~VCollisionReader() = 0;

    virtual AAMCCinput operator()() = 0;

protected:
    VCollisionReader() = default;
};

inline VCollisionReader::~VCollisionReader() = default;
