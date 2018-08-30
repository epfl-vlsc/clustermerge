#pragma once

namespace agd {
namespace errors {
enum Code {
  OK,
  OUT_OF_RANGE,
  RESOURCE_EXHAUSTED,
  UNKNOWN,
  INTERNAL,
  NOT_FOUND,
  INVALID_ARGUMENT,
  UNAVAILABLE
};

}  // namespace errors
}  // namespace agd
