
#include <zlib.h>
#include <fstream>
#include "src/agd/buffer.h"

class BufferedCompressor {
 public:
  BufferedCompressor(uint32_t buf_size);
  ~BufferedCompressor();

  // initialize the compressor
  // must call before Write()
  int Init(const std::string& filename);

  // write bytes
  // compressor will buffer bytes to buf_size, then compress and flush to file
  int Write(const char* data, size_t sz);

 private:
  // data into buffer_, buffer_ is compressed into compressed_buffer_,
  // compressed_buffer_ dumped to file
  std::ofstream outfile_;
  agd::Buffer buffer_;
  agd::Buffer compressed_buffer_;

  z_stream stream_ = {0};

  int compress_and_flush(const char* segment, size_t segment_size);

  void ensure_extend_capacity(size_t capacity);

  int finish();
};