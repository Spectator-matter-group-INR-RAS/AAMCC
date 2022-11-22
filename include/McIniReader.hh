//
// Created by apbus_amp_k on 12.11.22.
//

#ifndef GRATE_MCINIREADER_HH
#define GRATE_MCINIREADER_HH

#include "URun.h"
#include "UEvent.h"
#include "UParticle.h"
#include "EventInitialState.h"

#include "VCollisionReader.hh"
#include "AAMCC.hh"

class McIniReader final : public VCollisionReader{
public:
    McIniReader() = delete;
    McIniReader(McIniReader&&) noexcept;
    McIniReader& operator=(McIniReader&&) noexcept;

    McIniReader(const std::unique_ptr<TFile>&);

    ~McIniReader() final = default;

    AAMCCinput operator()() final;
    // temporary methods to provide EventInitialState and UEvenet readings to outside TODO: after implementing normal DTOs should be rewritten to actually return an appropriate DTO
    EventInitialState* getInStateAddress() {return curr_st;}
    UEvent* getEventAdress() {return curr;}

private:
    TTree* ftree;           // Root data structure
    UEvent* curr;           // Event info from file cache
    EventInitialState* curr_st;  //Event initial state
    unsigned short asum;    // Total atomic mass for indexing purposes
    size_t iter;
    size_t treen;
};


#endif //GRATE_MCINIREADER_HH
