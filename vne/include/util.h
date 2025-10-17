// Copyright 2013 Leonardo Moura
#pragma once
#include <chrono>

template <typename T>
double truncateBetween(T x, T xmin, T xmax);

template<typename T>
class Timer{
  typedef typename T::rep R;
public:
	Timer();
  Timer(R time_limit);
	void start();
	void stop();
  bool isRunning() const;
	R total() const;
  bool reachedTimeLimit() const;
private:
  std::chrono::high_resolution_clock::time_point start_, finish_;
  bool running_;
  R time_limit_;
};

bool DBL_EQL(double, double);
bool DBL_GE (double, double);
bool DBL_G (double, double);
bool DBL_L (double, double);

