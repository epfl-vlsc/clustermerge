
#include "debug.h"

std::string PrintNormalizedProtein(const char* seq, size_t len) {
  std::vector<char> scratch;
  scratch.resize(len);
  memcpy(&scratch[0], seq, len);
  for (size_t i = 0; i < len; i++) {
    scratch[i] = scratch[i] + 'A';
  }
  return std::string(&scratch[0], len);
}
