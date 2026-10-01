#include "png.h"

#include <algorithm>
#include <fstream>

namespace mo {

static uint32_t crc32(const uint8_t* d, size_t n, uint32_t crc = 0) {
    static uint32_t table[256];
    static bool init = false;
    if (!init) {
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            table[i] = c;
        }
        init = true;
    }
    crc = ~crc;
    for (size_t i = 0; i < n; ++i) crc = table[(crc ^ d[i]) & 0xFF] ^ (crc >> 8);
    return ~crc;
}

static void put32(std::vector<uint8_t>& v, uint32_t x) {
    for (int s = 24; s >= 0; s -= 8) v.push_back(uint8_t(x >> s));
}

static void chunk(std::vector<uint8_t>& out, const char* type, const std::vector<uint8_t>& data) {
    put32(out, uint32_t(data.size()));
    std::vector<uint8_t> td(type, type + 4);
    td.insert(td.end(), data.begin(), data.end());
    out.insert(out.end(), td.begin(), td.end());
    put32(out, crc32(td.data(), td.size()));
}

namespace {
struct BitWriter {
    std::vector<uint8_t>& out;
    uint32_t acc = 0;
    int n = 0;
    explicit BitWriter(std::vector<uint8_t>& o) : out(o) {}
    void bits(uint32_t v, int cnt) {  // LSB d'abord
        acc |= v << n;
        n += cnt;
        while (n >= 8) { out.push_back(uint8_t(acc & 0xFF)); acc >>= 8; n -= 8; }
    }
    void huff(uint32_t code, int len) {  // code de Huffman : MSB d'abord
        uint32_t r = 0;
        for (int i = 0; i < len; ++i) r |= ((code >> i) & 1u) << (len - 1 - i);
        bits(r, len);
    }
    void flush() { if (n > 0) { out.push_back(uint8_t(acc & 0xFF)); acc = 0; n = 0; } }
};

void putLiteral(BitWriter& w, int lit) {
    if (lit < 144) w.huff(0x30 + lit, 8);
    else if (lit < 256) w.huff(0x190 + lit - 144, 9);
    else if (lit < 280) w.huff(lit - 256, 7);
    else w.huff(0xC0 + lit - 280, 8);
}

// Deflate à codes de Huffman fixes + LZ77 (chaînes de hachage).
std::vector<uint8_t> deflateFixed(const std::vector<uint8_t>& d) {
    static const int lenBase[29] = {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258};
    static const int lenExtra[29] = {0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0};
    static const int distBase[30] = {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
    static const int distExtra[30] = {0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13};
    std::vector<uint8_t> out = {0x78, 0x01};
    BitWriter w(out);
    w.bits(1, 1);  // dernier bloc
    w.bits(1, 2);  // Huffman fixe
    const int HB = 15, HS = 1 << HB, WIN = 32768, MAXCHAIN = 24;
    std::vector<int> head(HS, -1), prev(d.size(), -1);
    auto hash = [&](size_t i) { return ((d[i] << 10) ^ (d[i + 1] << 5) ^ d[i + 2]) & (HS - 1); };
    size_t i = 0;
    while (i < d.size()) {
        int bestLen = 0, bestDist = 0;
        if (i + 2 < d.size()) {
            int h = hash(i);
            int cand = head[h], chain = 0;
            while (cand >= 0 && int(i) - cand <= WIN && chain++ < MAXCHAIN) {
                int l = 0, maxl = int(std::min<size_t>(258, d.size() - i));
                while (l < maxl && d[cand + l] == d[i + l]) ++l;
                if (l > bestLen) { bestLen = l; bestDist = int(i) - cand; if (l == maxl) break; }
                cand = prev[cand];
            }
        }
        size_t adv = 1;
        if (bestLen >= 3) {
            int lc = 28;
            while (lenBase[lc] > bestLen) --lc;
            putLiteral(w, 257 + lc);
            if (lenExtra[lc]) w.bits(uint32_t(bestLen - lenBase[lc]), lenExtra[lc]);
            int dc = 29;
            while (distBase[dc] > bestDist) --dc;
            w.huff(uint32_t(dc), 5);
            if (distExtra[dc]) w.bits(uint32_t(bestDist - distBase[dc]), distExtra[dc]);
            adv = size_t(bestLen);
        } else {
            putLiteral(w, d[i]);
        }
        for (size_t k = 0; k < adv; ++k, ++i)
            if (i + 2 < d.size()) { int h = hash(i); prev[i] = head[h]; head[h] = int(i); }
    }
    putLiteral(w, 256);
    w.flush();
    return out;
}
}  // namespace

std::vector<uint8_t> encodePNG(int w, int h, const std::vector<uint8_t>& rgb) {
    std::vector<uint8_t> raw;
    raw.reserve(size_t(h) * (w * 3 + 1));
    for (int y = 0; y < h; ++y) {
        raw.push_back(1);  // filtre « Sub »
        const uint8_t* row = rgb.data() + size_t(y) * w * 3;
        for (int x = 0; x < w * 3; ++x) raw.push_back(uint8_t(row[x] - (x >= 3 ? row[x - 3] : 0)));
    }
    std::vector<uint8_t> z = deflateFixed(raw);
    uint32_t a = 1, b = 0;
    for (uint8_t c : raw) { a = (a + c) % 65521; b = (b + a) % 65521; }
    put32(z, (b << 16) | a);
    std::vector<uint8_t> out = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    std::vector<uint8_t> ihdr;
    put32(ihdr, w);
    put32(ihdr, h);
    ihdr.insert(ihdr.end(), {8, 2, 0, 0, 0});
    chunk(out, "IHDR", ihdr);
    chunk(out, "IDAT", z);
    chunk(out, "IEND", {});
    return out;
}

bool writePNG(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    auto b = encodePNG(w, h, rgb);
    f.write(reinterpret_cast<const char*>(b.data()), std::streamsize(b.size()));
    return true;
}

std::string base64(const std::vector<uint8_t>& d) {
    static const char* T = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string s;
    s.reserve((d.size() + 2) / 3 * 4);
    for (size_t i = 0; i < d.size(); i += 3) {
        uint32_t v = d[i] << 16;
        if (i + 1 < d.size()) v |= d[i + 1] << 8;
        if (i + 2 < d.size()) v |= d[i + 2];
        s += T[(v >> 18) & 63];
        s += T[(v >> 12) & 63];
        s += i + 1 < d.size() ? T[(v >> 6) & 63] : '=';
        s += i + 2 < d.size() ? T[v & 63] : '=';
    }
    return s;
}

}  // namespace mo
