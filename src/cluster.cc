
#include "cluster.h"

agd::Status Cluster::AlignReps(const Cluster& other,
                               ProteinAligner::Alignment* alignment,
                               ProteinAligner* aligner) {
  const auto& this_rep = seqs_[0];
  const auto& other_rep = other.seqs_[0];
  return aligner->AlignSingle(this_rep.Seq().data(), other_rep.Seq().data(),
                              this_rep.Seq().size(), other_rep.Seq().size(),
                              *alignment);
}

bool Cluster::PassesThreshold(const Cluster& other, ProteinAligner* aligner) {
  const auto& this_rep = seqs_[0];
  const auto& other_rep = other.seqs_[0];

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
    residue_total += seq.Seq().size();
    seqs_.push_back(std::move(seq));
  }
}

void Cluster::Merge(Cluster* other, ProteinAligner* aligner) {
  const auto& other_seqs = other->Sequences();
  seqs_.push_back(other_seqs[0]);  // the rep matches, or we wouldnt be here

  bool first = true;  // to skip first
  for (const auto& seq : other_seqs) {
    if (first) {
      first = false;
      continue;
    }
    const auto& rep = seqs_[0];
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
        residue_total += seq.Seq().size();
        seqs_.push_back(seq);
      }
    }
  }

  other->MergeOther(this, aligner);
}

void Cluster::MergeOther(Cluster* other, ProteinAligner* aligner) {
  const auto& other_seqs = other->Sequences();
  seqs_.push_back(other_seqs[0]);  // the rep matches, or we wouldnt be here

  bool first = true;  // to skip first
  for (const auto& seq : other_seqs) {
    if (first) {
      first = false;
      continue;
    }
    const auto& rep = seqs_[0];
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
        residue_total += seq.Seq().size();
        seqs_.push_back(seq);
      }
    }
  }
}
