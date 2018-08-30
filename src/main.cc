
#include <iostream>
#include "agd/agd_dataset.h"

using std::cout;
using std::string;
using std::unique_ptr;

string PrintNormalizedProtein(const char* seq, size_t len) {
  std::vector<char> scratch;
  scratch.resize(len);
  memcpy(&scratch[0], seq, len);
  for (size_t i = 0; i < len; i++) {
    scratch[i] = scratch[i] + 'A';
  }
  return string(&scratch[0], len);
}

int main(int argc, char** argv) {
  string json_path(argv[1]);
  unique_ptr<agd::AGDDataset> dataset;

  cout << "JSON path is " << json_path << "\n";
  agd::Status s = agd::AGDDataset::Create(json_path, dataset);

  if (!s.ok()) {
    cout << s.ToString() << "\n";
    return 0;
  }

  agd::AGDDataset::ColumnIterator iter;
  s = dataset->Column("prot", &iter);
  
  if (!s.ok()) {
    cout << s.ToString() << "\n";
    return 0;
  }

  const char* data;
  size_t size;
  s = iter.GetNextRecord(&data, &size);
  while (s.ok()) {
    cout << PrintNormalizedProtein(data, size) << "\n\n";
    s = iter.GetNextRecord(&data, &size);
  }

  return (0);
}
