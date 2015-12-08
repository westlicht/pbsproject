#pragma once

#include "Common.h"

#include <string>
#include <cstring>
#include <ctime>

namespace pbs {
namespace StringUtils {

static std::string lower(const std::string &str) {
    std::string result(str);
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

static std::string upper(const std::string &str) {
    std::string result(str);
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

static std::string formatTime(const std::tm *tm, const std::string &format) {
    char buf[1024];
    std::strftime(buf, sizeof(buf) - 1, format.c_str(), tm);
    buf[sizeof(buf) - 1] = '\0';
    return std::string(buf);
}

// trim from start
static inline std::string ltrim(const std::string &str) {
    std::string ret(str);
    ret.erase(ret.begin(), std::find_if(ret.begin(), ret.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return ret;
}

// trim from end
static inline std::string rtrim(const std::string &str) {
    std::string ret(str);
    ret.erase(std::find_if(ret.rbegin(), ret.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), ret.end());
    return ret;
}

// trim from both ends
static inline std::string trim(const std::string &str) {
    return ltrim(rtrim(str));
}

static bool endsWith(const std::string &str, const char *ending) {
    size_t strLen = str.length();
    size_t endingLen = std::strlen(ending);
    return strLen >= endingLen && 0 == str.compare(strLen - endingLen, endingLen, ending);
}

static bool endsWith(const std::string &str, const std::string &ending) {
    size_t strLen = str.length();
    size_t endingLen = ending.length();
    return strLen >= endingLen && 0 == str.compare(strLen - endingLen, endingLen, ending);
}

static bool startsWith(const std::string &str, const std::string &start) {
    return start.size() <= str.size() && str.compare(0, start.size(), start) == 0;
}

} // namespace StringUtils
} // namespace pbs
