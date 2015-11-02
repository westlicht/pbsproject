#pragma once

#include "Common.h"
#include "Timer.h"

#include <string>
#include <deque>
#include <vector>
#include <numeric>

namespace pbs {

// Simple profiler.
class Profiler {
public:
    struct Item {
        std::string name;
        double avg;
        std::deque<double> history;

        Item(const std::string &name) : name(name), avg(0.0) {}

        void enter() {
            timer.reset();
            active = true;
        }

        void leave() {
            if (active) {
                history.emplace_back(timer.elapsed());
                while (history.size() > 10) { history.pop_front(); }
                avg = std::accumulate(history.begin(), history.end(), 0.0) / history.size();
            }
            active = false;
        }

    private:
        Timer timer;
        bool active = false;
    };

    static void enter(const std::string &name) {
        item(name).enter();
    }
    static void leave(const std::string &name) {
        item(name).leave();
    }

    template<typename Func>
    static void profile(const std::string &name, Func func) {
        enter(name);
        func();
        leave(name);
    }

    static const std::vector<Item> &items() {
        return _items;
    }

    static void dump() {
        for (const auto &item : _items) {
            DBG("%-20s %.1f ms", item.name, item.avg);
        }
    }

private:
    static Item &item(const std::string &name);
    static std::vector<Item> _items;
};

// Profiles the time spent within the current scope.
class ProfileScope {
public:
    ProfileScope(const std::string &name) : _name(name) {
        Profiler::enter(_name);
    }

    ~ProfileScope() {
        Profiler::leave(_name);
    }

private:
    std::string _name;
};

} // namespace pbs
