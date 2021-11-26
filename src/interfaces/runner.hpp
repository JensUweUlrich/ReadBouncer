/*
 * read_holder.h
 *
 *  Created on: October 13, 2021
 *  
 */
#pragma once
#include <atomic>

namespace interfaces {

// TODO: define as singleton class cause only one runner object can exist in the application
    struct Runner
    {
        std::atomic_bool isRunning = true;
    };


} /* namespace interfaces */
