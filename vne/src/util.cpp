// Copyright 2012 Leonardo Moura
#include <cmath>

#include "../include/util.h"
using std::chrono::high_resolution_clock;

template <typename T>
double truncateBetween(T x, T xmin, T xmax) {
  return max(min(x, xmax), xmin);
}

template <typename T>
Timer<T>::Timer() : start_(high_resolution_clock::now()),
    running_(true)  { }

template <typename T>
Timer<T>::Timer(Timer<T>::R time_limit) : start_(high_resolution_clock::now()),
    running_(true), time_limit_(time_limit) { }

template <typename T>
void Timer<T>::start() { start_ = high_resolution_clock::now(); running_ = true; }

template <typename T>
void Timer<T>::stop() { finish_ = high_resolution_clock::now(); running_ = false; }

template <typename T>
typename Timer<T>::R Timer<T>::total() const {
  if(running_)
    return std::chrono::duration_cast<T>(high_resolution_clock::now() - start_).count();
  else
    return std::chrono::duration_cast<T>(finish_ - start_).count();
}

template <typename T>
bool Timer<T>::reachedTimeLimit() const {
  return total() > time_limit_;
}

template <typename T>
bool Timer<T>::isRunning() const {
  return running_;
}

#define EPSILON 1E-6
bool DBL_EQL(double a, double b) {
    return fabs(a - b) < EPSILON;
}

bool DBL_GE(double a, double b) {
    return a - b > -EPSILON;
}

bool DBL_G (double a, double b) {
    return a - b > EPSILON && fabs(a-b) > EPSILON;
}

bool DBL_L (double a, double b) {
    return a - b < EPSILON && fabs(a-b) > EPSILON;
}

template class Timer<std::chrono::milliseconds>;
template class Timer<std::chrono::seconds>;
