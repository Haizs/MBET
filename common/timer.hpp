#pragma once

#include <chrono>
#include <iostream>

class SteadyClockWatch {
 public:
  void Start() { start_ = std::chrono::steady_clock::now(); }

  double Stop() {
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - start_).count();
  }

 private:
  std::chrono::steady_clock::time_point start_;
};

template <auto ID>
class Timer {
  using Watch = SteadyClockWatch;

  class StopWatch : private Watch {
   public:
    StopWatch() { Start(); }

    void Start() { Watch::Start(); }

    double Stop() {
      elapsed_time_ = Watch::Stop();
      sum_time_ += elapsed_time_;
      num_calls_++;
      return elapsed_time_;
    }

    void Report(const std::string& label) {
      std::cout << "[TIME] " << label << ": " << elapsed_time_ << " sec" << std::endl;
    }

   private:
    double elapsed_time_;
  };

  class AutoWatch : private Watch {
   public:
    AutoWatch() { Watch::Start(); }

    ~AutoWatch() {
      sum_time_ += Watch::Stop();
      num_calls_++;
    }
  };

 public:
  static StopWatch StopWatch() { return {}; }
  static AutoWatch AutoWatch() { return {}; }

  static void Report(const std::string& label) {
    std::cout << "[TIME] " << label << ": cnt=" << num_calls_ << ", time=" << sum_time_ << "sec" << std::endl;
  }

  static void Reset() {
    num_calls_ = 0;
    sum_time_ = 0;
  }

 private:
  static size_t num_calls_;
  static double sum_time_;
};

template <auto ID>
size_t Timer<ID>::num_calls_ = 0;

template <auto ID>
double Timer<ID>::sum_time_ = 0;

template <auto ID>
class DebugTimer : public Timer<ID> {
  class DummyWatch {
   public:
    constexpr void Start() {}
    constexpr void Stop() {}
    constexpr void Report(const std::string& label) {}
  };

#ifndef DEBUG
 public:
  static constexpr DummyWatch StopWatch() { return {}; }
  static constexpr DummyWatch AutoWatch() { return {}; }

  static constexpr void Report(const std::string& label) {}
  static constexpr void Reset() {}
#endif
};

using MainTimer = Timer<0>;
using InitTimer = Timer<1>;
using RunTimer = Timer<2>;