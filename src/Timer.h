#pragma once

#include "Common.h"

#include <chrono>

namespace pbs {

// Simple timer with millisecond resolution.
class Timer {
public:
    Timer() {
        reset();
    }

    void reset() {
        _start = std::chrono::system_clock::now();
    }

    // Returns elapsed milliseconds
    double elapsed() const {
        auto now = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - _start);
        return double(duration.count());
    }

    // Returns elapsed time as a string
    std::string elapsedString(bool precise = false) const {
        return timeString(elapsed(), precise);
    }

private:
    std::chrono::system_clock::time_point _start;
};

} // namespace pbs
