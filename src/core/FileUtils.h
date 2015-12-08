#pragma once

#include "Common.h"
#include "StringUtils.h"

#include <cassert>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#ifdef __WIN32
#include <windows.h>
#else // __WIN32
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#endif // __WIN32

#include <sys/stat.h>

namespace pbs {
namespace FileUtils {

static const char PATH_DELIMITER = '/';

static std::string dirname(const std::string &filename) {
    char buf[1024];
    assert(filename.length() < sizeof(buf) - 1);
#ifdef __WIN32
    _splitpath(filename.c_str(), NULL, buf, NULL, NULL);
    std::string result(buf);
    result.pop_back();
    return result;
#else
    memcpy(buf, filename.c_str(), filename.length());
    buf[filename.length()] = '\0';
    return std::string(::dirname(buf));
#endif
}

static std::string basename(const std::string &filename) {
    char buf[1024];
    assert(filename.length() < sizeof(buf) - 1);
#ifdef __WIN32
    _splitpath(filename.c_str(), NULL, NULL, buf, NULL);
    return std::string(buf);
#else
    memcpy(buf, filename.c_str(), filename.length());
    buf[filename.length()] = '\0';
    return std::string(::basename(buf));
#endif
}

static std::string realpath(const std::string &filename) {
#ifdef __WIN32
    char buf[1024];
    return std::string(_fullpath(buf, filename.c_str(), sizeof(buf)));
#else
    char *path = ::realpath(filename.c_str(), NULL);
    if (!path)
        return filename;
    std::string result = std::string(path);
    free(path);
    return result;
#endif
}

static std::pair<std::string, std::string> splitFilename(const std::string &path) {
    unsigned found = path.find_last_of("/\\");
    return std::pair<std::string, std::string>(path.substr(0, found), path.substr(found+1));
}

static std::pair<std::string, std::string> splitExtension(const std::string &path) {
    unsigned found = path.find_last_of('.');
    return std::pair<std::string, std::string>(path.substr(0, found), path.substr(found+1));
}

static std::string extractFilename(const std::string &path) {
    return splitFilename(path).second;
}

static std::string extractDirectory(const std::string &path) {
    return splitFilename(path).first;
}

static std::string extractExtension(const std::string &path) {
    return splitExtension(path).second;
}

static bool hasExtension(const std::string &path, const std::string &ext) {
    if (ext.empty())
        return false;
    const std::string fileExt = StringUtils::upper(extractExtension(path));
    const std::string checkExt = StringUtils::upper(ext[0] == '.' ? ext.substr(1) : ext);
    return fileExt == checkExt;
}

static std::string replaceExtension(const std::string &path, const std::string &ext) {
    return splitExtension(path).first + '.' + ext;
}

static std::string appendTrailingDelimiter(const std::string &filename) {
    if (filename.length() == 0 || filename.back() != PATH_DELIMITER)
        return filename + PATH_DELIMITER;
    return filename;
}

static std::string join(const std::string &dir, const std::string &filename) {
    return appendTrailingDelimiter(dir) + filename;
}

static bool fileExists(const std::string &filename) {
    std::ifstream is(filename, std::ifstream::in | std::ifstream::binary);
    return is.good();
}

static std::string readFile(const std::string &filename) {
    std::ifstream t(filename);
    std::string str;

    t.seekg(0, std::ios::end);   
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(t)),
                std::istreambuf_iterator<char>());
    
    return str;
}

static bool dirExists(const std::string &dirname) {
#ifdef __WIN32
    DWORD fattr = GetFileAttributesA(realpath(dirname).c_str());
    return 
        (fattr != INVALID_FILE_ATTRIBUTES) &&
        (fattr & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat sb;
    return 
        (!stat(dirname.c_str(), &sb)) &&
        (S_ISDIR(sb.st_mode));
#endif
}

static bool createDir(const std::string &dirname) {
#ifdef __WIN32
    return CreateDirectory(dirname.c_str(), NULL);
#else
    return mkdir(dirname.c_str(), 0755) != -1;
#endif
}

static bool deleteFile(const std::string &filename) {
#ifdef __WIN32
    return DeleteFile(filename.c_str());
#else
    return !unlink(filename.c_str());
#endif
}

static bool deleteDir(const std::string &dirname) {
#ifdef __WIN32
    return RemoveDirectory(dirname.c_str());
#else
    return !rmdir(dirname.c_str());
#endif    
}

static bool changeCurrentDir(const std::string &dir) {
#ifdef __WIN32
    return SetCurrentDirectoryA(dir.c_str()) != 0;
#else
    return chdir(dir.c_str()) == 0;
#endif
    return false;
}

static std::string getCurrentDir() {
    char buf[1024];
#ifdef __WIN32
    DWORD size = GetCurrentDirectory(sizeof(buf), buf);
    if (size >= sizeof(buf)) {
        std::unique_ptr<char[]> tmpBuf(new char[size + 1]);
        size = GetCurrentDirectory(size + 1, tmpBuf.get());
        if (size) {
            return appendTrailingDelimiter(std::string(tmpBuf.get(), size));
        }
    } else if (size != 0) {
        return appendTrailingDelimiter(std::string(buf, size));
    }
#else
    if (getcwd(buf, sizeof(buf)))
        return appendTrailingDelimiter(std::string(buf));
#endif
    return std::string();
}
    
} // namespace FileUtils
} // namespace pbs
