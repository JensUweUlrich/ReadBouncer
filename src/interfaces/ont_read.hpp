/*
 * read_holder.h
 *
 *  Created on: October 13, 2021
 *  
 */
#pragma once

#include <map>
#include <list>
#include <string>
#include <memory>
#include <vector>
#include <StopClock.hpp>

#ifndef INTERFACES_READHOLDER_
#define INTERFACES_READHOLDER_

namespace interfaces {

enum ClassifyAction { unblock, stop_receiving };

/// Representation of read data and associated analysis results.
class ONTRead {
public:
   
    uint32_t channelNr{};
    uint32_t readNr{};
    uint32_t readTag{};
    std::string id{};
    std::vector<float> raw_signals{};
    std::string sequence{};
    bool unblock{};

    /// Default constructor.
    ONTRead() {}

    ONTRead(const bool unblock, const ONTRead& rhs);

    // Remove ordinary copy options.
    ONTRead(const ONTRead& rhs);// = delete;

    ONTRead& operator=(const ONTRead& rhs);// = delete;

    /// Move constructor.
    ONTRead(ONTRead&& rhs);

    /// Move assignment operator.
    ONTRead& operator=(ONTRead&& rhs);

    /// Virtual destructor.
    virtual ~ONTRead();

    /// Clears the object.
    void clear();

};

std::ostream& operator<<(std::ostream& strm, const ONTRead& rhs);

typedef std::pair<ONTRead, TimeMeasures> RTPair;

} /* namespace interfaces */

#endif //INTERFACES_READHOLDER_