#pragma once

#include "Common.h"

#include <istream>
#include <ostream>
#include <vector>

namespace pbs {
namespace Serialize {

template<typename T>
static void write(std::ostream &os, const T &value)
{
    os.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template<typename T>
static void read(std::istream &is, T &value)
{
    is.read(reinterpret_cast<char *>(&value), sizeof(T));
}

template<typename T>
static void writeVector(std::ostream &os, const std::vector<T> &vec)
{
    std::size_t size = vec.size();
    write(os, size);
    os.write(reinterpret_cast<const char *>(vec.data()), size * sizeof(T));
}

template<typename T>
static void readVector(std::istream &is, std::vector<T> &vec)
{
    std::size_t size;
    read(is, size);
    vec.resize(size);
    is.read(reinterpret_cast<char *>(vec.data()), size * sizeof(T));
}

template<typename Matrix>
static void writeMatrix(std::ostream &os, const Matrix &m)
{
    size_t rows = m.rows(), cols = m.cols();
    write(os, rows);
    write(os, cols);
    os.write(reinterpret_cast<const char *>(m.data()), rows * cols * sizeof(typename Matrix::Scalar));
}

template<typename Matrix>
static void readMatrix(std::istream &is, Matrix &m)
{
    size_t rows, cols;
    read(is, rows);
    read(is, cols);
    m.resize(rows, cols);
    is.read(reinterpret_cast<char *>(m.data()), rows * cols * sizeof(typename Matrix::Scalar));
}

} // namespace Serialize
} // namespace pbs

