#pragma once

// 3D morton code.
class Morton3D {
public:
    static inline uint32_t interleave10bit(uint32_t x) {
        x &= 0x3ff;
        x = (x | x << 16) & 0x30000ff;
        x = (x | x << 8) & 0x300f00f;
        x = (x | x << 4) & 0x30c30c3;
        x = (x | x << 2) & 0x9249249;
        return x;
    }

    static inline uint64_t interleave21bit(uint64_t x) {
        x &= 0x1fffff;
        x = (x | x << 32) & 0x1f00000000ffff;
        x = (x | x << 16) & 0x1f0000ff0000ff;
        x = (x | x << 8) & 0x100f00f00f00f00f;
        x = (x | x << 4) & 0x10c30c30c30c30c3;
        x = (x | x << 2) & 0x1249249249249249;
        return x;
    }

    static inline uint32_t morton10bit(uint32_t x, uint32_t y, uint32_t z) {
        return interleave10bit(x) | (interleave10bit(y) << 1) | (interleave10bit(z) << 1);
    }

    static inline uint64_t morton21bit(uint64_t x, uint64_t y, uint64_t z) {
        return interleave21bit(x) | (interleave21bit(y) << 1) | (interleave21bit(z) << 1);
    }

};
