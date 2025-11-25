/*
 * File Name: my_measure.h
 * Description: Custom measurement functions for t-J model
 *
 */

#ifndef TJMODEL_SRC_MY_MEASURE_H
#define TJMODEL_SRC_MY_MEASURE_H

#include "qlmps/one_dim_tn/mps/finite_mps/finite_mps.h"
#include "qlmps/one_dim_tn/mps/finite_mps/finite_mps_measu_memory.h"
#include "qlten/qlten.h"

#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <map>
#include <utility>

namespace qlmps {
using namespace qlten;

//helper
// FindLeftBoundary has been removed as it is obsolete in current UltraDMRG context.

/**
 * Measure bond correlations (two-point functions) for a given list of bonds.
 * Optimized for measuring many pairs by sweeping through the MPS.
 * 
 * @param mps The FiniteMPS object (tensors stored on disk).
 * @param mps_path Path to MPS tensors.
 * @param bonds List of bonds (pairs of sites) to measure. Assumes site1 < site2.
 * @param op1 Operator at site1.
 * @param op2 Operator at site2.
 * @param inst_op Optional insertion operator (e.g. Jordan-Wigner string) applied between site1 and site2.
 * @return MeasuRes containing the results.
*/
template<typename TenElemT, typename QNT>
MeasuRes<TenElemT> MeasureBonds(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::string &mps_path,
    const std::vector<std::pair<size_t, size_t>> &bonds,
    const QLTensor<TenElemT, QNT> &op1,
    const QLTensor<TenElemT, QNT> &op2,
    const QLTensor<TenElemT, QNT> &inst_op = QLTensor<TenElemT, QNT>()
) {
  using TenT = QLTensor<TenElemT, QNT>;
  
  // Group bonds by site1
  std::map<size_t, std::vector<size_t>> bonds_map;
  for (const auto &bond : bonds) {
    size_t s1 = std::min(bond.first, bond.second);
    size_t s2 = std::max(bond.first, bond.second);
    bonds_map[s1].push_back(s2);
  }

  // Sort site2 for each site1 to process in order
  for (auto &entry : bonds_map) {
    std::sort(entry.second.begin(), entry.second.end());
  }

  MeasuRes<TenElemT> results;
  results.reserve(bonds.size());

  // Prepare MPS state
  // We assume the canonical center is at 0 (UltraDMRG convention).
  // We start scanning from the first site that has bonds.
  
  if (bonds_map.empty()) return results;

  size_t start_site = bonds_map.begin()->first;
  
  // Helper to move center and manage memory
  // We need to be at site1 (Left Canonical up to site1-1)

  assert(mps.empty());
  size_t current_center = 0;

  // Move center to start_site
  // We assume the center is at 0 initially.
  // We can just move forward using LeftCanonicalize.

  auto id_op_set = mps.GetSitesInfo().id_ops;
  bool has_inst = (inst_op != TenT());

  for (auto it = bonds_map.begin(); it != bonds_map.end(); ++it) {
    size_t site1 = it->first;
    const auto &targets = it->second;

    // Move center to site1 if needed
    while (current_center < site1) {
        if (mps(current_center) == nullptr) {
            mps.LoadTen(current_center, GenMPSTenName(mps_path, current_center));
        }
        mps.LoadTen(current_center + 1, GenMPSTenName(mps_path, current_center + 1));
        mps.LeftCanonicalizeTen(current_center);
        mps.dealloc(current_center);
        current_center++;
    }
    // Now current_center == site1 and mps[site1] is loaded and is the canonical center
    
    // Load site1 if it wasn't moved to (e.g. site1 == 0)
    if(mps(site1) == nullptr) {
        mps.LoadTen(site1, GenMPSTenName(mps_path, site1));
    }

    // Contract at site1
    // temp = mps[site1] * op1
    TenT temp_ten0;
    Contract(&mps[site1], &op1, {{1}, {0}}, &temp_ten0);
    
    // ptemp = temp * dag(mps[site1])
    // This creates the environment after site1
    TenT *ptemp_ten = new TenT();
    TenT mps_ten_dag = Dag(mps[site1]);
    Contract(&temp_ten0, &mps_ten_dag, {{0, 2}, {0, 1}}, ptemp_ten);
    
    // We keep mps[site1] loaded because it is the current center. 
    // We will move the center from here in the next iteration of site1.
    
    size_t eated_site = site1;
    
    for (size_t site2 : targets) {
        // Move from eated_site to site2
        // Apply inst_op or identity
        while (eated_site < site2 - 1) {
            size_t eating_site = eated_site + 1;
            mps.LoadTen(eating_site, GenMPSTenName(mps_path, eating_site));
            
            const auto &op_mid = has_inst ? inst_op : id_op_set[eating_site];
            
            CtrctMidTen(mps, eating_site, op_mid, id_op_set[eating_site], ptemp_ten);
            
            mps.dealloc(eating_site);
            eated_site = eating_site;
        }
        
        // Now eated_site == site2 - 1
        // Load site2
        mps.LoadTen(site2, GenMPSTenName(mps_path, site2));
        
        // Contract site2
        // ptemp * mps[site2] * op2 * dag(mps[site2])
        
        TenT temp_ten2, temp_ten3, res_ten;
        Contract(&mps[site2], ptemp_ten, {{0}, {0}}, &temp_ten2);
        
        Contract(&temp_ten2, &op2, {{0}, {0}}, &temp_ten3);
        
        TenT mps_ten_dag2 = Dag(mps[site2]);
        Contract(&temp_ten3, &mps_ten_dag2,
                 {{0, 1, 2}, {2, 0, 1}},
                 &res_ten);
        
        results.push_back(MeasuResElem<TenElemT>({site1, site2}, res_ten()));
        
        // We don't need to explicitly dealloc site2 here if we want to rely on OS cache,
        // but we do dealloc it in standard usage unless it's the center.
        // However, for correctness of 'eated_site' logic in next iteration, 
        // site2 is just treated as the end of THIS measurement.
        // For the next target (if any), 'eated_site' logic will re-load site2 (as intermediate) if needed.
        mps.dealloc(site2);
    } // end targets loop
    
    delete ptemp_ten;
    
    // Before moving to next site1, we need to make sure site1 is handled?
    // We are at `site1`. We contracted it.
    // Loop `it` moves to next `site1_next`.
    // At start of loop, we `while(current_center < site1_next)`.
    // This will `LeftCanonicalize` site1 and dealloc it.
    // So we are good.

  } // end site1 loop
  
  // Cleanup final center
  if(mps(current_center) != nullptr) {
      mps.dealloc(current_center);
  }

  return results;
}

} // namespace qlmps

#endif // TJMODEL_SRC_MY_MEASURE_H
