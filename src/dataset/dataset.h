
#pragma once

#include "src/agd/status.h"

using agd::Status;

// dataset interface
// data, sz for protein seqs, meta, meta_sz for metadata
class Dataset {
  public: 

    /*virtual static Status Create(const std::string file_path,
                              std::unique_ptr<FastaDataset>& dataset) = 0;*/
    virtual ~Dataset() {}

    uint32_t Size() const { return total_records_; }

    const std::string& Name() const { return genome_; }

    virtual Status GetNextRecord(const char** data, size_t* sz,
                              const char** meta = nullptr, size_t* meta_sz = nullptr) = 0;

    virtual Status GetRecordAt(const size_t index, const char** data, size_t* sz,
                            const char** meta = nullptr, size_t* meta_sz = nullptr) const = 0;

    virtual void Reset() { current_record_ = 0; }

  protected: 
    std::string genome_;
    size_t total_records_;
    size_t current_record_;
};