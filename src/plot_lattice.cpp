#include <array>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "honeycomb_lattice.h"
#include "squarelattice.h"

namespace {

constexpr double kDefaultScaleSquare = 120.0;
constexpr double kDefaultScaleHoneycomb = 160.0;

struct GeometryRecord {
  std::string name;
  std::function<void(size_t, size_t)> generator;
};

std::string BuildOutputName(const std::string &geometry, size_t Ly, size_t Lx) {
  return geometry + "_Ly" + std::to_string(Ly) + "_Lx" + std::to_string(Lx) + ".svg";
}

void EnsureOutputDir(const std::string &path) {
  std::filesystem::path output(path);
  if (output.has_parent_path()) {
    std::filesystem::create_directories(output.parent_path());
  }
}

std::function<std::string(size_t, size_t)> HoneycombColorizer(
    const HoneycombLattice &lattice) {
  return [&lattice](size_t a, size_t b) {
    const auto dir = lattice.ClassifyNearestNeighbor(a, b);
    return lattice.DirectionColor(dir);
  };
}

void DumpSquare(const std::string &filename, const LatticeBase &lattice,
                const std::vector<LatticeBase::Link> &winding) {
  lattice.DumpSVG(filename, lattice.nearest_neighbor_links, winding,
                  kDefaultScaleSquare);
  std::cout << "[SVG] Wrote " << filename << '\n';
}

void DumpHoneycomb(const std::string &filename, const HoneycombLattice &lattice) {
  auto color_fn = HoneycombColorizer(lattice);
  lattice.DumpSVG(filename, lattice.nearest_neighbor_links, {},
                  kDefaultScaleHoneycomb, color_fn);
  std::cout << "[SVG] Wrote " << filename << '\n';
}

void HandleSquareLattices(size_t Ly, size_t Lx) {
  {
    SquareCylinder lattice(Ly, Lx);
    const auto filename = BuildOutputName("SquareCylinder", Ly, Lx);
    EnsureOutputDir(filename);
    DumpSquare(filename, lattice, lattice.wind_nn_links);
  }

  {
    SquareOBC lattice(Ly, Lx);
    const auto filename = BuildOutputName("SquareOBC", Ly, Lx);
    EnsureOutputDir(filename);
    DumpSquare(filename, lattice, {});
  }

  {
    SquareRotatedCylinder lattice(Ly, Lx);
    const auto filename = BuildOutputName("SquareRotated", Ly, Lx);
    EnsureOutputDir(filename);
    DumpSquare(filename, lattice, lattice.wind_nn_links);
  }

  {
    SquareTorus lattice(Ly, Lx);
    const auto filename = BuildOutputName("SquareTorus", Ly, Lx);
    EnsureOutputDir(filename);
    DumpSquare(filename, lattice, lattice.wind_nn_links);
  }

  {
    SquareLadder lattice(Lx);
    const auto filename = BuildOutputName("SquareLadder", lattice.Ly, lattice.Lx);
    EnsureOutputDir(filename);
    DumpSquare(filename, lattice, {});
  }

  {
    Chain lattice(Lx);
    const auto filename = BuildOutputName("SquareChain", lattice.Ly, lattice.Lx);
    EnsureOutputDir(filename);
    DumpSquare(filename, lattice, lattice.wind_nn_links);
  }
}

void HandleHoneycombLattices(size_t Ly, size_t Lx) {
  {
    HoneycombXCCylinder lattice(Ly, Lx);
    const auto filename = BuildOutputName("HoneycombXC", Ly, Lx);
    EnsureOutputDir(filename);
    DumpHoneycomb(filename, lattice);
  }

  {
    HoneycombYCCylinder lattice(Ly, Lx);
    const auto filename = BuildOutputName("HoneycombYC", Ly, Lx);
    EnsureOutputDir(filename);
    DumpHoneycomb(filename, lattice);
  }
}

std::pair<size_t, size_t> ParseDimensions(int argc, char *argv[]) {
  if (argc == 1) {
    return {4, 8};
  }
  if (argc != 3) {
    throw std::runtime_error(
        "Usage: ./plot_lattice [Ly Lx]\nExample: ./plot_lattice 4 8");
  }
  const size_t Ly = static_cast<size_t>(std::stoul(argv[1]));
  const size_t Lx = static_cast<size_t>(std::stoul(argv[2]));
  if (Ly == 0 || Lx == 0) {
    throw std::runtime_error("Ly and Lx must be positive integers.");
  }
  return {Ly, Lx};
}

} // namespace

int main(int argc, char *argv[]) {
  try {
    const auto [Ly, Lx] = ParseDimensions(argc, argv);
    HandleSquareLattices(Ly, Lx);
    HandleHoneycombLattices(Ly, Lx);
    std::cout << "[SVG] Generation completed.\n";
  } catch (const std::exception &ex) {
    std::cerr << "[error] " << ex.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
