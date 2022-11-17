#pragma once
#include "AAMCC.hh"

class VReader{
public:
    virtual void operator()(AAMCCrun* run, InitialConditions* InCond) = 0;
    virtual ~VReader() = 0;
};

inline VReader::~VReader() = default;