
#pragma once

#include <queue>
#include "src/common/cluster_set.h"
#include "src/common/alignment_environment.h"
#include "src/common/aligner.h"
#include "src/common/params.h"

void MergeBatch(std::deque<ClusterSet>& sets_to_merge, ProteinAligner* aligner);
