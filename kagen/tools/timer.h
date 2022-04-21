/*******************************************************************************
 * include/tools/timer.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _TIMER_H_
#define _TIMER_H_

#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

namespace kagen {

class Timer {
public:
    Timer() {
        start_ = timestamp();
    }

    void Restart() {
        start_ = timestamp();
    }

    double Elapsed() {
        return timestamp() - start_;
    }

private:
    /** Returns a timestamp ('now') in seconds (incl. a fractional part). */
    inline double timestamp() {
        struct timeval tp;
        gettimeofday(&tp, NULL);
        return double(tp.tv_sec) + tp.tv_usec / 1000000.;
    }

    double start_;
};

} // namespace kagen
#endif
