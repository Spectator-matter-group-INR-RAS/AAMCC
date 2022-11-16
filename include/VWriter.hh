#pragma once
#include "AAMCC.hh"

class VWriter{
public:
    virtual void operator()(AAMCCEvent* ev, AAMCCrun* run, aamcc::NucleonVector* nucleons) = 0;
    virtual ~VWriter() = 0;
};

inline VWriter::~VWriter() = default;