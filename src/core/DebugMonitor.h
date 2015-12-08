#pragma once

#include <tinyformat.h>

#include <string>
#include <vector>

namespace pbs {

class DebugMonitor {
public:
    struct Item {
        std::string name;
        std::string value;
    };

    static void clear();
    static void addItem(const std::string &name, const std::string &value);
    static const std::vector<Item> &items() { return _items; }

    template<typename... Args>
    static inline void addItem(const std::string &name, const char *fmt, const Args &... args) {
        addItem(name, tfm::format(fmt, args...));
    }

private:
    static std::vector<Item> _items;
};

} // namespace pbs
