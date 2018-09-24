
#include "cluster.h"
#include <iostream>

agd::Status Cluster::AlignReps(const Cluster& other,
                               ProteinAligner::Alignment* alignment,
                               ProteinAligner* aligner) {
  const auto& this_rep = seqs_[0];
  const auto& other_rep = other.seqs_[0];
  return aligner->AlignSingle(this_rep.Seq().data(), other_rep.Seq().data(),
                             this_rep.Seq().size(), other_rep.Seq().size(),
                             *alignment);
}

bool Cluster::PassesThreshold(const Cluster& other, ProteinAligner* aligner, SequenceIDMap* id_map) {
  const auto& this_rep = seqs_[0];
  const auto& other_rep = other.seqs_[0];

  bool passed = false;
  SequencePair sp;
  if (this_rep.ID() < other_rep.ID()) {
    sp = std::make_pair(this_rep.ID(), other_rep.ID());
  } else {
    sp = std::make_pair(other_rep.ID(), this_rep.ID());
  }

  bool found = id_map->Exists(sp, &passed);
  if (found) {
    return passed;
  } else {
    passed = aligner->PassesThreshold(this_rep.Seq().data(), other_rep.Seq().data(),
                                  this_rep.Seq().size(),
                                  other_rep.Seq().size());
    id_map->Insert(sp, passed);
    return passed;
  }
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
    seqs_.push_back(std::move(seq));
  }
}

void Cluster::Merge(const Cluster& other, ProteinAligner* aligner, SequenceIDMap* id_map) {

  const auto& other_seqs = other.Sequences();
  seqs_.push_back(other_seqs[0]); // the rep matches, or we wouldnt be here

  bool first = true; // to skip first
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
      bool passed = false;
      SequencePair sp;
      if (seq.ID() < rep.ID()) {
        sp = std::make_pair(seq.ID(), rep.ID());
      } else {
        sp = std::make_pair(rep.ID(), seq.ID());
      }

      bool found = id_map->Exists(sp, &passed);
      if (found && passed) {
        std::cout << "avoided a recompute!!\n";
        seqs_.push_back(seq);
      } else {
        passed = aligner->PassesThreshold(rep.Seq().data(), seq.Seq().data(),
                                      rep.Seq().size(),
                                      seq.Seq().size());
        id_map->Insert(sp, passed);
        if (passed) {
          seqs_.push_back(seq);
        }
      }
    }
  }
}
