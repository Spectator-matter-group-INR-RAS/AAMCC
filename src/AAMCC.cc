#include "AAMCC.hh"

AAMCC::AAMCC() {

}

AAMCC::~AAMCC() = default;

template<class CollisionGeneratorOut>
G4ReactionProductVector* AAMCC::Do(CollisionGeneratorOut *CollisionOut) {
    __assert(std::is_base_of<VCollisionReader, CollisionGeneratorOut>::value, "type parameter of this method must be derived form VCollisionReader");
    aamcc::NucleonVector nucleons = CollisionOut->operator();
    return nullptr;
}


