#include "Profiler.h"

namespace pbs {

std::vector<Profiler::Item> Profiler::_items;

Profiler::Item &Profiler::item(const std::string &name) {
    auto item = std::find_if(_items.begin(), _items.end(), [&name] (const Item &item) { return item.name == name; });
    if (item == _items.end()) {
        _items.emplace_back(Item(name));
        return _items.back();
    } else {
        return *item;
    }
}


} // namespace pbs
