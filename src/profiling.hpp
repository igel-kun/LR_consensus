

/** \file profiling.hpp
 * collection of profiling functions (timers, ...)
 */

#ifndef PROFILING_HPP
#define PROFILING_HPP

#include <unistd.h>
#include <sys/times.h>
#include <time.h>

//! get the number of cycles this program has gone through
double get_cycles(){
  return clock();
//  struct tms tms;
//  if ( times( &tms ) != (clock_t)-1 )
//    return (double)tms.tms_utime;
}

//! a simple timer class that can be paused and resumed
struct timer{
  double start_time;  //!< last time the timer was started
  double stop_time;   //!< last time the timer was stopped
  
  //! constructor
  timer(): start_time(0), stop_time(0) {}

  //! start the timer
  void start(){
    start_time = get_cycles();
  }

  //! stop the timer
  void stop(){
    stop_time = get_cycles();
  }

  //! pause the timer (= stop())
  void pause(){
    stop_time = get_cycles();
  }

  //! resume the timer by setting an artificial start_time
  void resume(){
    start_time += get_cycles() - stop_time;
  }

  //! return the number of cycles passed between last start and last stop
  /** the timer has to have been started and stopped for this to make sense */
  double cycles_passed() const {
    return stop_time - start_time;
  }

  //! return the number of seconds passed between last start and last stop (see cycles_passed())
  double seconds_passed() const {
//    return cycles_passed() / (double)sysconf( _SC_CLK_TCK );
    return cycles_passed() / CLOCKS_PER_SEC;
  }
};

#endif

