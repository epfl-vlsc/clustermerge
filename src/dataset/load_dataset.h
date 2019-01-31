
#pragma once

#include "agd_protein_dataset.h"
#include "args.h"
#include "fasta_dataset.h"

// utilities to facilitate loading datasets

agd::Status LoadDataset(const std::string& ext, const std::string& filename,
                        std::unique_ptr<Dataset>& dataset);

agd::Status LoadDatasetsPositional(
    args::PositionalList<std::string>& datasets_opts,
    std::vector<std::unique_ptr<Dataset>>* datasets);

agd::Status LoadDatasetsJSON(const std::string& dataset_file,
                             std::vector<std::unique_ptr<Dataset>>* datasets);
