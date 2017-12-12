/******************************************************************************
 * Timer.h
 *
 * Source of the graph generator
 ******************************************************************************
 * Copyright (C) 2016 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef _TIMER_H_
#define _TIMER_H_

#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

class Timer {
 public:
  Timer() { start_ = timestamp(); }

  void Restart() { start_ = timestamp(); }

  double Elapsed() { return timestamp() - start_; }

 private:
  /** Returns a timestamp ('now') in seconds (incl. a fractional part). */
  inline double timestamp() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return double(tp.tv_sec) + tp.tv_usec / 1000000.;
  }

  double start_;
};

#endif
