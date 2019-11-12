
#include "cluster.h"
#include <iostream>

agd::Status Cluster::AlignReps(const Cluster& other,
                               ProteinAligner::Alignment* alignment,
                               ProteinAligner* aligner) {
  const auto& this_rep = seqs_.front();
  const auto& other_rep = other.seqs_.front();

  return aligner->AlignSingle(all_seqs_->at(this_rep).Seq().data(), all_seqs_->at(other_rep).Seq().data(),
                              all_seqs_->at(this_rep).Seq().size(), all_seqs_->at(other_rep).Seq().size(),
                              *alignment);
}

bool Cluster::PassesThreshold(const Cluster& other, ProteinAligner* aligner) {
  const auto& this_rep = seqs_.front();
  const auto& other_rep = other.seqs_.front();

  return aligner->PassesThreshold(all_seqs_->at(this_rep).Seq().data(), all_seqs_->at(other_rep).Seq().data(),
                              all_seqs_->at(this_rep).Seq().size(), all_seqs_->at(other_rep).Seq().size());
}

void Cluster::AddSequence(uint32_t seq) {
  // make sure we aren't adding duplicate
  bool found = false;
  for (const auto& s : seqs_) {
    if (s == seq){
      found = true;
      break;
    }
  }
  if (!found) {
    seqs_.push_back(seq);
  }
}

void Cluster::Merge(Cluster* other, ProteinAligner* aligner) {
  const auto& other_seqs = other->Sequences();
  //seqs_.push_back(other_seqs.front());  // the rep matches, or we wouldnt be here
  bool found = false;
  for (const auto& s : seqs_) {
    if (s == other_seqs.front()) {
      found = true;
      break;
    }
  }
  if (!found) {
    seqs_.push_back(other_seqs.front());  // the rep matches, or we wouldnt be here
  } 

  bool first = true;  // to skip first
  for (const auto& seq : other_seqs) {
    if (first) {
      first = false;
      continue;
    }
    const auto& rep = seqs_.front();
    bool found = false;
    for (const auto& s : seqs_) {
      if (s == seq) {
        found = true;
        break;
      }
    }
    if (!found) {
      if (aligner->PassesThreshold(all_seqs_->at(rep).Seq().data(), all_seqs_->at(seq).Seq().data(),
                                   all_seqs_->at(rep).Seq().size(), all_seqs_->at(seq).Seq().size())) {
        seqs_.push_back(seq);
      }
    }
  }

  other->MergeOther(this, aligner);
}

void Cluster::MergeOther(Cluster* other, ProteinAligner* aligner) {
  const auto& other_seqs = other->Sequences();
  bool found = false;
  for (const auto& s : seqs_) {
    if (s == other_seqs.front()) {
      found = true;
      break;
    }
  }
  if (!found) {
    seqs_.push_back(other_seqs.front());  // the rep matches, or we wouldnt be here
  } 

  bool first = true;  // to skip first
  for (const auto& seq : other_seqs) {
    if (first) {
      first = false;
      continue;
    }
    const auto& rep = seqs_.front();
    bool found = false;
    for (const auto& s : seqs_) {
      if (s == seq) {
        found = true;
        break;
      }
    }
    if (!found) {
      if (aligner->PassesThreshold(all_seqs_->at(rep).Seq().data(), all_seqs_->at(seq).Seq().data(),
                                   all_seqs_->at(rep).Seq().size(), all_seqs_->at(seq).Seq().size())) {
        seqs_.push_back(seq);
      }
    }
  }
}
