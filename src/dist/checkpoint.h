
// checkpoint file format and IO
// Stuart Byma, EPFL

#include "src/agd/status.h"
#include "src/common/concurrent_queue.h"
#include "src/comms/requests.h"

// file is a blob
// just a dump of serialized clusterSets in the sets to merge queue

agd::Status WriteCheckpointFile(const std::string& path,
    const ConcurrentQueue<MarshalledClusterSet>& queue);

agd::Status LoadCheckpoinFile(const std::string& path,
    const ConcurrentQueue<MarshalledClusterSet>& queue);

bool CheckpointFileExists();