
#include "src/common/buffered_compressor.h"
#include "absl/strings/str_cat.h"
#include <iostream>

using namespace std;

static const int ENABLE_ZLIB_GZIP = 32;
static const int ENABLE_ZLIB_GZIP_COMPRESS = 16;
static const int window_bits = 15;
static const size_t extend_length = 1024 * 1024 * 8; // 2 Mb
static const size_t reserve_factor = 3;


BufferedCompressor::BufferedCompressor(uint32_t buf_size) {
  buffer_.reserve(buf_size);
  compressed_buffer_.reserve(buf_size);
}

BufferedCompressor::~BufferedCompressor() {
  int res = finish();
  if (res) {
    cout << "BufferedCompressor finish failed" << endl;
  }
  outfile_.close();
}

int BufferedCompressor::Init(const std::string& filename) {
  string compressed_name = absl::StrCat(filename, ".gz");
  cout << "BufferedCompressor opening file: " << compressed_name << std::endl;
  outfile_.open(compressed_name, std::ios::out);

  if (!outfile_.good()) {
    cout << "failed to open file " << filename << endl;
    return -1;
  }

  stream_ = {0};
  // Not sure if this is necessary
  stream_.zalloc = Z_NULL;
  stream_.zfree = Z_NULL;
  stream_.opaque = Z_NULL;

  int status = deflateInit2(&stream_, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                            window_bits | ENABLE_ZLIB_GZIP_COMPRESS,
                            9,  // higher memory, better speed
                            Z_DEFAULT_STRATEGY);

  if (status != Z_OK) {
    std::cout << "deflateInit() didn't return Z_OK. Return " << status
              << " with 2nd param " << Z_DEFAULT_COMPRESSION << std::endl;
  }

  return 0;
}

int BufferedCompressor::compress_and_flush(const char* segment,
                                           size_t segment_size) {

  // compress `buffer_` into `compressed_buffer_` and flush to file 

  compressed_buffer_.reset();
  stream_.avail_out = 0;
  stream_.next_in = const_cast<unsigned char*>(
      reinterpret_cast<const unsigned char*>(segment));
  stream_.avail_in = segment_size;

  size_t bytes = stream_.total_out;

  int status;
  while (stream_.avail_in != 0) {
    if (stream_.avail_out == 0)
      ensure_extend_capacity((stream_.avail_in / 2) +
                             512);  // in case of round off at the end
    status = deflate(&stream_, Z_NO_FLUSH);
    if (status != Z_OK) {
      cout << "deflate(Z_NO_FLUSH) return status " << status << std::endl;
      return -1;
    }
    compressed_buffer_.resize(stream_.total_out - bytes);
  }

  // according to the documentation, this is the assumption when avail_out > 0
  // and all input has been consumed
  if (stream_.avail_in != 0) {
    std::cout << "Compressor: stream.avail in > 0! Got " << stream_.avail_in
              << std::endl;
    return -1;
  }

  cout << "compressed and flushed " << stream_.total_out - bytes << " bytes, avail out is " << stream_.avail_out << endl;

  // the chunk is compressed, dump it to file
  outfile_.write(compressed_buffer_.data(), compressed_buffer_.size());

  return 0;
}

int BufferedCompressor::finish() {
  cout << "buffered compressor finishing ..." << endl;
  if (buffer_.size() > 0) {
    if (compress_and_flush(buffer_.data(), buffer_.size())) {
      return -1;
    }
  }

  compressed_buffer_.reset();
  size_t bytes = stream_.total_out;

  int status;
  stream_.next_in = nullptr;
  stream_.avail_in = 0;
  do {
    ensure_extend_capacity(32);
    status = deflate(&stream_, Z_FINISH);
    if (status != Z_STREAM_END) {
      cout << "deflate(Z_FINISH) return status " << status << endl;
      return -1;
    }
    compressed_buffer_.resize(stream_.total_out - bytes);
  } while (status == Z_OK);

  // flush all remaining output to the output buffer
  status = deflateEnd(&stream_);
  if (status != Z_OK) {
    cout << "deflateEnd() didn't receive Z_OK. Got: " << status << endl;
    return -1;
  }

  compressed_buffer_.resize(stream_.total_out - bytes);
  cout << "finish compression still produced " << stream_.total_out - bytes << endl;
  
  outfile_.write(compressed_buffer_.data(), compressed_buffer_.size());
  return 0;
}

void BufferedCompressor::ensure_extend_capacity(size_t capacity) {
  size_t avail_capacity = compressed_buffer_.capacity() - compressed_buffer_.size();
  if (avail_capacity < capacity) {
    compressed_buffer_.extend_allocation(capacity - avail_capacity);
  }

  size_t sz = compressed_buffer_.size();
  stream_.avail_out = compressed_buffer_.capacity() - sz;
  stream_.next_out = const_cast<unsigned char*>(
      reinterpret_cast<const unsigned char*>(&compressed_buffer_[sz]));
}

int BufferedCompressor::Write(const char* data, size_t sz) {
  if (sz > buffer_.capacity() - buffer_.size()) {
    // it will not fit, compress the buffer and flush to file
    // res is 0 on success
    int res = compress_and_flush(buffer_.data(), buffer_.size());
    if (res) {
      cout << "unable to compress and flush " << endl;
      return res;
    }
    buffer_.reset();
  }

  if (buffer_.AppendBuffer(data, sz) != agd::Status::OK()) {
    return -1;
  } else {
    return 0;
  }
}