// Encodeur PNG minimal (RVB 8 bits, deflate « stocké ») et base64.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace mo {
std::vector<uint8_t> encodePNG(int w, int h, const std::vector<uint8_t>& rgb);
bool writePNG(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb);
std::string base64(const std::vector<uint8_t>& data);
}  // namespace mo
