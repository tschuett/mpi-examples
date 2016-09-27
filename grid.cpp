#include "grid.hpp"

#include <cassert>

using namespace std;

Grid::Grid(std::array<uint32_t, 3> _dims, uint32_t _rank, uint32_t size)
    : dims(_dims), rank(_rank) {
  assert(dims[0] * dims[1] * dims[2] == size);

  coords = rank2coords(rank);
}

array<uint32_t, 3> Grid::rank2coords(uint32_t rank) const {
  array<uint32_t, 3> res;
  res[2] = rank / (dims[0] * dims[1]);
  res[1] = (rank - res[2] * dims[0] * dims[1]) / dims[0];
  res[0] = rank - res[1] * dims[0] - res[2] * dims[0] * dims[1];

  return res;
}

uint32_t Grid::coords2rank(array<uint32_t, 3> coords) const {
  return coords[0] + dims[0] * coords[1] + dims[0] * dims[1] * coords[2];
}

uint32_t Grid::xplus() const {
  return coords2rank({(coords[0] + 1) % dims[0], coords[1], coords[2]});
}

uint32_t Grid::xminus() const {
  return coords2rank({(coords[0] - 1) % dims[0], coords[1], coords[2]});
}

uint32_t Grid::yplus() const {
  return coords2rank({coords[0], (coords[1] + 1) % dims[1], coords[2]});
}

uint32_t Grid::yminus() const {
  return coords2rank({coords[0], (coords[1] - 1) % dims[1], coords[2]});
}

uint32_t Grid::zplus() const {
  return coords2rank({coords[0], coords[1], (coords[2] + 1) % dims[2]});
}

uint32_t Grid::zminus() const {
  return coords2rank({coords[0], coords[1], (coords[2] - 1) % dims[2]});
}

Grid get_3d_grid(uint32_t size, uint32_t rank) {
  switch (size) {
  case 5440:
    return Grid({16, 17, 20}, rank, size);
  case 5372:
    return Grid({4, 17, 79}, rank, size);
  case 5056:
    return Grid({8, 8, 79}, rank, size);
  case 5120:
    return Grid({16, 16, 20}, rank, size);
  case 272:
    return Grid({4, 4, 17}, rank, size);
  case 256:
    return Grid({4, 8, 8}, rank, size);
  case 16:
    return Grid({4, 5, 8}, rank, size);
  case 136:
    return Grid({2, 4, 17}, rank, size);
  case 128:
    return Grid({4, 4, 8}, rank, size);
  case 80:
    return Grid({4, 4, 5}, rank, size);
  case 68:
    return Grid({2, 2, 17}, rank, size);
  case 64:
    return Grid({4, 4, 4}, rank, size);
  case 36:
    return Grid({2, 3, 6}, rank, size);
  case 24:
    return Grid({2, 2, 6}, rank, size);
  case 12:
    return Grid({2, 2, 3}, rank, size);
  default:
    assert(false);
  }
}
