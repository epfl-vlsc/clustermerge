
#include "cluster.h"
#include <iostream>

agd::Status Cluster::AlignReps(const Cluster& other,
                               ProteinAligner::Alignment* alignment,
                               ProteinAligner* aligner) {
  const auto& this_rep = seqs_.front();
  const auto& other_rep = other.seqs_.front();

  return aligner->AlignSingle(this_rep.Seq().data(), other_rep.Seq().data(),
                              this_rep.Seq().size(), other_rep.Seq().size(),
                              *alignment);
}

bool Cluster::PassesThreshold(const Cluster& other, ProteinAligner* aligner) {
  const auto& this_rep = seqs_.front();
  const auto& other_rep = other.seqs_.front();

  return aligner->PassesThreshold(this_rep.Seq().data(), other_rep.Seq().data(),
                                  this_rep.Seq().size(),
                                  other_rep.Seq().size());
}

void Cluster::AddSequence(const Sequence& seq) {
  // make sure we aren't adding duplicate
  bool found = false;
  for (const auto& s : seqs_) {
    if (s.ID() == seq.ID()) {
      found = true;
      break;
    }
  }
  if (!found) {
    residue_total_ += seq.Seq().size();
    if (seq.Seq().size() > longest_) {
      longest_ = seq.Seq().size();
    }
    seqs_.push_back(seq);
  }
}

void Cluster::Merge(Cluster* other, ProteinAligner* aligner) {
  const auto& other_seqs = other->Sequences();
  seqs_.push_back(other_seqs.front());  // the rep matches, or we wouldnt be here

  bool first = true;  // to skip first
  for (const auto& seq : other_seqs) {
    if (first) {
      first = false;
      continue;
    }
    const auto& rep = seqs_.front();
    bool found = false;
    for (const auto& s : seqs_) {
      if (s.ID() == seq.ID()) {
        found = true;
        break;
      }
    }
    if (!found) {
      if (aligner->PassesThreshold(rep.Seq().data(), seq.Seq().data(),
                                   rep.Seq().size(), seq.Seq().size())) {
        residue_total_ += seq.Seq().size();
        if (seq.Seq().size() > longest_) {
          longest_ = seq.Seq().size();
        }
        seqs_.push_back(seq);
      }
    }
  }

  other->MergeOther(this, aligner);
}

void Cluster::MergeOther(Cluster* other, ProteinAligner* aligner) {
  const auto& other_seqs = other->Sequences();
  seqs_.push_back(other_seqs.front());  // the rep matches, or we wouldnt be here

  bool first = true;  // to skip first
  for (const auto& seq : other_seqs) {
    if (first) {
      first = false;
      continue;
    }
    const auto& rep = seqs_.front();
    bool found = false;
    for (const auto& s : seqs_) {
      if (s.ID() == seq.ID()) {
        found = true;
        break;
      }
    }
    if (!found) {
      if (aligner->PassesThreshold(rep.Seq().data(), seq.Seq().data(),
                                   rep.Seq().size(), seq.Seq().size())) {
        residue_total_ += seq.Seq().size();
        if (seq.Seq().size() > longest_) {
          longest_ = seq.Seq().size();
        }
        seqs_.push_back(seq);
      }
    }
  }
}
