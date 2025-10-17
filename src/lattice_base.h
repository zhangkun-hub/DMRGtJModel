// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TJMODEL_SRC_LATTICE_BASE_H
#define TJMODEL_SRC_LATTICE_BASE_H

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

class LatticeBase {
 public:
  using Link = std::tuple<size_t, size_t>;

  LatticeBase() = default;
  LatticeBase(size_t Ly, size_t Lx, size_t sites_per_unit = 1)
      : Ly(Ly), Lx(Lx), N(Lx * Ly * sites_per_unit) {}
  virtual ~LatticeBase() = default;

  inline void Print() const {
    std::cout << "System size (Lx, Ly) = ( " << Lx << ", " << Ly << ")" << "\n";
  }

  virtual std::pair<double, double> SiteCoordinate(size_t site) const = 0;

  inline void DumpSVG(const std::string &filename,
                      const std::vector<Link> &links,
                      const std::vector<Link> &winding_links,
                      double scale = 60.0,
                      const std::function<std::string(size_t, size_t)> &color_fn =
                          [](size_t, size_t) { return std::string("#444444"); }) const {
    std::ofstream f(filename);
    if (!f.is_open()) {
      std::cerr << "[DumpSVG] Failed to open file: " << filename << std::endl;
      return;
    }
    f << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 1200 1200'>\n";

    auto emit_point = [&](size_t idx, const std::pair<double, double> &coord) {
      const double cx = coord.first * scale + 100.0;
      const double cy = 1000.0 - coord.second * scale;
      f << "  <circle cx='" << cx << "' cy='" << cy
        << "' r='6' fill='#1f77b4'/>\n";
      f << "  <text x='" << cx + 8 << "' y='" << cy + 4
        << "' font-size='16' fill='#333333'>" << idx << "</text>\n";
    };

    auto emit_line = [&](const Link &link, const std::string &style) {
      const auto [src, dst] = link;
      const auto src_coord = SiteCoordinate(src);
      const auto dst_coord = SiteCoordinate(dst);
      const double x1 = src_coord.first * scale + 100.0;
      const double y1 = 1000.0 - src_coord.second * scale;
      const double x2 = dst_coord.first * scale + 100.0;
      const double y2 = 1000.0 - dst_coord.second * scale;
      f << "  <line x1='" << x1 << "' y1='" << y1
        << "' x2='" << x2 << "' y2='" << y2
        << "' stroke-width='2' " << style << "/>\n";
    };

    auto emit_colored_line = [&](const Link &link, const std::string &color,
                                 bool dashed) {
      std::string style = "stroke='" + color + "'";
      if (dashed) {
        style += " stroke-dasharray='8,4'";
      }
      emit_line(link, style);
    };

    auto emit_winding_curve = [&](const Link &link, const std::string &color) {
      const auto [src, dst] = link;
      const auto src_coord = SiteCoordinate(src);
      const auto dst_coord = SiteCoordinate(dst);
      const double x1 = src_coord.first * scale + 100.0;
      const double y1 = 1000.0 - src_coord.second * scale;
      const double x2 = dst_coord.first * scale + 100.0;
      const double y2 = 1000.0 - dst_coord.second * scale;

      const double dx = std::abs(src_coord.first - dst_coord.first);
      const double dy = std::abs(src_coord.second - dst_coord.second);

      double ctrl_x = (x1 + x2) / 2.0;
      double ctrl_y = (y1 + y2) / 2.0;

      const double bend = 120.0;
      if (dy > dx) {
        ctrl_x += (x1 <= x2 ? 1 : -1) * bend;
      } else if (dx > dy) {
        ctrl_y += (y1 <= y2 ? -1 : 1) * bend;
      } else {
        ctrl_x += bend;
        ctrl_y += bend;
      }

      f << "  <path d='M " << x1 << " " << y1
        << " Q " << ctrl_x << " " << ctrl_y << " " << x2 << " " << y2
        << "' fill='none' stroke='" << color
        << "' stroke-width='2' stroke-dasharray='8,4'/>\n";
    };

    for (const auto &link : links) {
      const auto color = color_fn(std::get<0>(link), std::get<1>(link));
      emit_colored_line(link, color, false);
    }
    for (const auto &link : winding_links) {
      const auto color = color_fn(std::get<0>(link), std::get<1>(link));
      emit_winding_curve(link, color.empty() ? std::string("#d62728") : color);
    }
    for (size_t site = 0; site < N; ++site) {
      emit_point(site, SiteCoordinate(site));
    }
    f << "</svg>\n";
  }

  size_t Ly{};
  size_t Lx{};
  size_t N{};
  std::vector<Link> nearest_neighbor_links;
  std::vector<Link> next_nearest_neighbor_links;
  std::vector<Link> third_nearest_neighbor_links;
  std::vector<Link> direct_nn_links;
  std::vector<Link> wind_nn_links;
  std::vector<Link> direct_nnn_links;
  std::vector<Link> wind_nnn_links_1;
  std::vector<Link> wind_nnn_links_2;
};

#endif // TJMODEL_SRC_LATTICE_BASE_H


