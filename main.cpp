#include <iostream>
#include <map>

#include "common/bigraph.hpp"
#include "common/timer.hpp"

#ifdef _MBET
#include "src/MBET.hpp"
#endif

int main(int argc, char **argv) {
  assert(argc >= 2);

  auto file_path = argv[1];
  std::cout << "====================================================" << std::endl;

  auto watch = MainTimer::StopWatch();

  auto bi_graph = BiGraph::LoadFromFile(file_path);

  watch.Stop();
  watch.Report("BiGraph::LoadFromFile");

  if (argc == 2) return 0;

  auto run_times = std::stoi(argv[2]);

  double sum_time = 0;
  std::map<std::string, double> sum_times;

#ifdef _MBET
  for (int run_time = 1; run_time <= run_times; run_time++) {
    auto [cnt, time, times] = MBET(bi_graph).Run();
    std::cout << "[Run-" << run_time << "] cnt=" << cnt << " , time: " << time << std::endl;
    for (auto &[k, v] : times) sum_times[k] += v;
    sum_time += time;
  }
  std::cout << "[Average(" << run_times << ")] avg_time: " << sum_time / run_times << std::endl;
  for (auto &[k, v] : sum_times) {
    std::cout << "[Average(" << run_times << ")] avg_" << k << ": " << v / run_times << std::endl;
  }
#endif

  {
    std::ifstream status("/proc/self/status");
    for (std::string line; std::getline(status, line);) {
      if (line.rfind("VmPeak:", 0) == 0) {
        std::cout << "[MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
      }
      if (line.rfind("VmSize:", 0) == 0) {
        std::cout << "[MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
      }
      if (line.rfind("VmHWM:", 0) == 0) {
        std::cout << "[MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
      }
      if (line.rfind("VmRSS:", 0) == 0) {
        std::cout << "[MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
      }
    }
  }

  return 0;
}
