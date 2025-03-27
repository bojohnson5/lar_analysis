#include <chrono>
#include <iostream>

template <typename T> class ProgressBar {
private:
  T total;
  T current;
  int width;
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

public:
  ProgressBar(T total, int width = 50)
      : total(total), width(width), current(0),
        start_time(std::chrono::high_resolution_clock::now()) {}

  void update(T progress) {
    current = progress;
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(now - start_time).count();

    float ratio = static_cast<float>(current) / static_cast<float>(total);
    int completed = static_cast<int>(ratio * width);

    double estimated_total_time = (elapsed / ratio);
    double remaining_time = estimated_total_time - elapsed;

    // Ensure values make sense
    if (ratio == 0)
      remaining_time = 0;

    std::cout << "\r[";
    for (int i = 0; i < completed; ++i)
      std::cout << "#";
    for (int i = completed; i < width; ++i)
      std::cout << "-";
    std::cout << "] " << static_cast<int>(ratio * 100) << "%";

    std::cout << " | Elapsed: " << format_time(elapsed)
              << " | ETA: " << format_time(remaining_time);

    std::cout.flush();
  }

  void finish() {
    update(total);
    std::cout << std::endl;
  }

private:
  std::string format_time(double seconds) {
    int mins = static_cast<int>(seconds) / 60;
    int secs = static_cast<int>(seconds) % 60;
    return (mins > 0 ? std::to_string(mins) + "m " : "") +
           std::to_string(secs) + "s";
  }
};
