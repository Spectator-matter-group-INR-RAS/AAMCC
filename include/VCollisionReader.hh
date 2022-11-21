#pragma once
#include <memory>

#include "AAMCConstants.hh"

class VCollisionReader {
public:
    VCollisionReader(const VCollisionReader&) = delete;                 //Forbid to copy reader object.
    VCollisionReader& operator=(const VCollisionReader&) = delete;      //It is associated with the file and its state shouldn't be copied.
    virtual ~VCollisionReader() = 0;

    virtual AAMCCinput operator()() = 0;

protected:
    VCollisionReader() = default;                                       //Protected constructor to prevent instancing
};

inline VCollisionReader::~VCollisionReader() = default;
