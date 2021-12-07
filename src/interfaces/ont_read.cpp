/*
 * read_holder.cpp
 *
 *  Created on: October 13, 2021
 *  
 */

#include "ont_read.hpp"


namespace interfaces {

ONTRead::ONTRead(const ONTRead& rhs) :
    channelNr(rhs.channelNr),
    readNr(rhs.readNr),
    id(rhs.id),
    readTag(rhs.readTag),
    sequence(rhs.sequence),
    raw_signals(rhs.raw_signals),
    unblock(rhs.unblock)
{
//    raw_signals.swap(rhs.raw_signals);
}

ONTRead::ONTRead(const bool ubl, const ONTRead& rhs) :
    channelNr(rhs.channelNr),
    readNr(rhs.readNr),
    id(rhs.id),
    readTag(rhs.readTag),
    sequence(rhs.sequence),
    raw_signals(rhs.raw_signals),
    unblock{ubl}
{
}

ONTRead& ONTRead::operator=(const ONTRead& rhs)
{
    id = rhs.id;
    channelNr = rhs.channelNr;
    readNr = rhs.readNr;
    readTag = rhs.readTag;
    sequence = rhs.sequence;
    raw_signals = rhs.raw_signals;
    unblock = rhs.unblock;
    return *this;
}

ONTRead::ONTRead(ONTRead&& rhs) :
        channelNr(rhs.channelNr),
        readNr(rhs.readNr),
        id(rhs.id),
        readTag(rhs.readTag),
        sequence(rhs.sequence),
        unblock(rhs.unblock)
{
    raw_signals.swap(rhs.raw_signals);
    rhs.clear();
}


ONTRead& ONTRead::operator=(ONTRead&& rhs) {
    id = rhs.id;
    channelNr = rhs.channelNr;
    readNr = rhs.readNr;
    readTag = rhs.readTag;
    sequence = rhs.sequence;
    raw_signals.swap(rhs.raw_signals);
    rhs.clear();
    return *this;
}


ONTRead::~ONTRead() = default;


void ONTRead::clear() {
    id.clear();
    channelNr = -1;
    readNr = -1;
    readTag = 0;
    sequence.clear();
    raw_signals.clear();
}

std::ostream& operator<<(std::ostream& strm, const ONTRead& rhs) {
    strm << "ONTRead(id: " << rhs.id
        << ", channelNr: " << rhs.channelNr
        << ", readNr: " << rhs.readNr
        << ", readTag: " << rhs.readTag;
    if (!rhs.raw_signals.empty())
        strm << ", signal size: " << rhs.raw_signals.size() << ", signals: " << std::to_string(rhs.raw_signals[0]) << " " << std::to_string(rhs.raw_signals[1]) << " " << std::to_string(rhs.raw_signals[2])
        << " ... " << std::to_string(rhs.raw_signals[rhs.raw_signals.size() - 3]) << " "
        << std::to_string(rhs.raw_signals[rhs.raw_signals.size() - 2]) << " "
        << std::to_string(rhs.raw_signals[rhs.raw_signals.size() - 1]);
    if (rhs.sequence.size() > 20)
        strm << ", sequence: " << rhs.sequence.substr(0, 10) << "..." << rhs.sequence.substr(rhs.sequence.size() - 11, 10);
    strm << ")";
    return strm;
}


} /* namespace interfaces */
