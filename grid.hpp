#pragma once

#include <array>
#include <cassert>
#include <cstdint>

class Grid {
  std::array<uint32_t, 3> dims;
  std::array<uint32_t, 3> coords;
  uint32_t rank;

  std::array<uint32_t, 3> rank2coords(uint32_t rank) const;
  uint32_t coords2rank(std::array<uint32_t, 3> coords) const;

public:
  Grid(std::array<uint32_t, 3> dims, uint32_t rank, uint32_t size);

  uint32_t xplus() const;
  uint32_t xminus() const;

  uint32_t yplus() const;
  uint32_t yminus() const;

  uint32_t zplus() const;
  uint32_t zminus() const;
};

extern Grid get_3d_grid(uint32_t size, uint32_t rank);
