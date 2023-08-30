#ifndef HOMOLOGYCOMPUTER_HPP

#define HOMOLOGYCOMPUTER_HPP
#include "SparseMatrix.hpp"
#include "FilteredSimplicialComplex.hpp"
#include <set>
#include <map>
#include <unordered_map>
#include <tuple>
#include "globals.hpp"
#include <limits>

namespace matilda
{

    template <typename Coefficient, typename Row>
    class PersistentHomologyComputer
    {
    public:
        using SparseVector = std::map<int_fast64_t, Coefficient>;
        using Simplex = std::vector<int_fast64_t>;
        std::unordered_map<int_fast64_t, std::map<int_fast64_t, std::array<float, 2>>> bars;
        std::map<Simplex, int_fast64_t> simplices_lookup; /**Simplex added at step 5 in filtration order maps to 5
                                                          * In general the following equality should hold:
                                                          * simplices_lookup[fsc.simplices[fsc.simplices_indices[i]]] == i
                                                          */
        std::unordered_map<int_fast64_t, std::map<int_fast64_t, Row>> persistent_cycles;
        bool persistent_cycles_computed = false;
        /**Helper variables for homology computation*/
        std::unordered_map<uint64_t, std::unordered_map<uint64_t, std::vector<int_fast64_t>>> proto_bars;
        SparseMatrix<Coefficient, Row> store_matrix;
        SparseMatrix<Coefficient, Row> cycle_bounders;
        SparseMatrix<Coefficient, Row> boundary_matrix;
        std::set<int_fast64_t> unmatched;
        std::unordered_map<int_fast64_t, int_fast64_t> pivots_lookup; /** Holds indices of `store_matrix` rows having a rightmost
                                                                       * nonzero entry at a given column. E.g. `pivots_lookup[col] == r`
                                                                       * means that the `r`-th row of `store_matrix` has its rightmost
                                                                       *  nonzero entry at column `col`.
                                                                       */
        Row boundary(const Simplex &s)
        {
            Row result;
            Simplex face;
            for (int_fast64_t i = 0; i < s.size(); ++i)
            {
                face = s;
                face.erase(face.begin() + i);
                if (i % 2 == 0)
                {
                    result.emplace(simplices_lookup[face], Coefficient(1));
                }
                else
                {
                    result.emplace(simplices_lookup[face], Coefficient(-1));
                }
            }
            return result;
        }

        int_fast64_t compute_persistent_homology(const FilteredSimplicialComplex &fsc,
                                                 int_fast64_t top_degree,
                                                 bool verbose,
                                                 bool partial = false,
                                                 int_fast64_t partial_start = 0,
                                                 int_fast64_t partial_end = -1,
                                                 int_fast64_t modulus_param = 3)
        {
            if (!partial || partial_start == 0)
            {
                set_modulus(modulus_param);
                std::copy(fsc.simplices_indices.begin(), fsc.simplices_indices.end(), std::inserter(unmatched, unmatched.begin()));
                for (int_fast64_t i = 0; i < fsc.simplices.size(); ++i)
                {
                    simplices_lookup[fsc.simplices[fsc.simplices_indices[i]]] = i;
                }
                partial_start = 0;
            }
            if (!partial || partial_end == -1)
            {
                partial_end = fsc.simplices.size();
            }
            std::vector<int_fast64_t> current_simplex;
            Coefficient coeff;
            int_fast64_t current_pivot, current_pivot_lookup;
            std::unordered_map<int_fast64_t, int_fast64_t>::iterator next_pivot_it = pivots_lookup.end();
            typename std::unordered_map<int_fast64_t, Row>::iterator current_row_it;
            typename std::unordered_map<int_fast64_t, Row>::iterator current_pivot_lookup_row_it;
            bool current_row_deleted = false;
            for (int_fast64_t current_row = partial_start, s_index; current_row < partial_end; ++current_row)
            {
                s_index = fsc.simplices_indices[current_row];
                if (fsc.simplices[s_index].size() > 1)
                {
                    current_row_it = store_matrix.rows.emplace(std::piecewise_construct,
                                                               std::forward_as_tuple(current_row),
                                                               std::forward_as_tuple(boundary(fsc.simplices[s_index])))
                                         .first;

                    // naive attempt
                    boundary_matrix.rows.emplace(std::piecewise_construct,
                                                 std::forward_as_tuple(current_row),
                                                 std::forward_as_tuple(boundary(fsc.simplices[s_index])));

                    current_row_deleted = false;
                    next_pivot_it = pivots_lookup.find(max_iterator(current_row_it->second)->first);
                }
                else
                {
                    current_row_deleted = true;
                }
                cycle_bounders.rows.emplace(std::piecewise_construct,
                                            std::forward_as_tuple(current_row),
                                            std::forward_as_tuple(current_row, Coefficient(1)));
                while (!current_row_deleted and next_pivot_it != pivots_lookup.end())
                {
                    current_pivot = max_iterator(current_row_it->second)->first;
                    current_pivot_lookup = next_pivot_it->second;
                    current_pivot_lookup_row_it = store_matrix.rows.find(current_pivot_lookup);
                    coeff = -(current_row_it->second[current_pivot] / current_pivot_lookup_row_it->second[current_pivot]);
                    current_row_deleted = store_matrix.addscale_from_to(current_pivot_lookup_row_it, current_row_it, coeff);
                    cycle_bounders.addscale_from_to(current_pivot_lookup, current_row, coeff);
                    if (!current_row_deleted)
                    {
                        next_pivot_it = pivots_lookup.find(max_iterator(current_row_it->second)->first);
                    }
                }
                if (!current_row_deleted)
                {
                    int_fast64_t start = max_iterator(current_row_it->second)->first;
                    pivots_lookup[start] = current_row;
                    proto_bars[fsc.simplices[s_index].size() - 2][start] = std::vector<int_fast64_t>{start, current_row};
                    unmatched.erase(start);
                    unmatched.erase(current_row);
                }
            }
            if (partial_end < fsc.simplices.size() - 1)
            {
                return 1;
            }
            /** Add infinite proto bars*/
            for (const auto i : unmatched)
            {
                proto_bars[fsc.simplices[fsc.simplices_indices[i]].size() - 1][i] = std::vector<int_fast64_t>{i, -1};
            }
            /** End of Add infinite proto bars*/

            /** Compute generators */
            for (const auto &bars_at_dimension : proto_bars)
            {
                for (const auto &bar : bars_at_dimension.second)
                {
                    if (bar.second[1] != -1 && fsc.appears_at[fsc.simplices_indices[bar.second[0]]] == fsc.appears_at[fsc.simplices_indices[bar.second[1]]])
                    {
                        continue;
                    }
                    Row result;
                    for (const auto &cycle_vector : cycle_bounders.rows[bar.second[0]])
                    {
                        if (fsc.simplices[fsc.simplices_indices[bar.second[0]]].size() == 1)
                        {
                            result.emplace(bar.second[0], Coefficient(1));
                        }
                        else
                        {
                            result.emplace(cycle_vector.first, cycle_vector.second);
                        }
                    }
                    for (auto it = result.begin(); it != result.end();)
                    {
                        if (it->second == 0)
                        {
                            it = result.erase(it);
                        }
                        else
                        {
                            ++it;
                        }
                    }
                    persistent_cycles[bars_at_dimension.first][bar.first] = result;
                }
            }
            persistent_cycles_computed = true;
            /** End of Compute generators    */

            /** Compute bars with filtration values*/
            for (const auto &bars_at_dimension : proto_bars)
            {
                for (const auto &bar : bars_at_dimension.second)
                {
                    if (bar.second[1] != -1)
                    {
                        if (fsc.appears_at[fsc.simplices_indices[bar.second[0]]] < fsc.appears_at[fsc.simplices_indices[bar.second[1]]])
                        {
                            bars[bars_at_dimension.first][bar.first][0] = fsc.appears_at[fsc.simplices_indices[bar.second[0]]];
                            bars[bars_at_dimension.first][bar.first][1] = fsc.appears_at[fsc.simplices_indices[bar.second[1]]];
                        }
                    }
                    else
                    {
                        bars[bars_at_dimension.first][bar.first][0] = fsc.appears_at[fsc.simplices_indices[bar.second[0]]];
                        bars[bars_at_dimension.first][bar.first][1] = std::numeric_limits<float>::infinity();
                    }
                }
            }
            /** End of Compute bars with filtration values*/
            return 0;
        }
        int_fast64_t compute_persistent_homology_no_representatives(const FilteredSimplicialComplex &fsc,
                                                                    int_fast64_t top_degree,
                                                                    bool verbose,
                                                                    bool partial = false,
                                                                    int_fast64_t partial_start = 0,
                                                                    int_fast64_t partial_end = -1,
                                                                    int_fast64_t modulus_param = 3)
        {
            if (!partial || partial_start == 0)
            {
                set_modulus(modulus_param);
                std::copy(fsc.simplices_indices.begin(), fsc.simplices_indices.end(), std::inserter(unmatched, unmatched.begin()));
                for (int_fast64_t i = 0; i < fsc.simplices.size(); ++i)
                {
                    simplices_lookup[fsc.simplices[fsc.simplices_indices[i]]] = i;
                }
                partial_start = 0;
            }
            if (!partial || partial_end == -1)
            {
                partial_end = fsc.simplices.size();
            }
            std::vector<int_fast64_t> current_simplex;
            Coefficient coeff;
            int_fast64_t current_pivot, current_pivot_lookup;
            std::unordered_map<int_fast64_t, int_fast64_t>::iterator next_pivot_it = pivots_lookup.end();
            typename std::unordered_map<int_fast64_t, Row>::iterator current_row_it;
            typename std::unordered_map<int_fast64_t, Row>::iterator current_pivot_lookup_row_it;
            bool current_row_deleted = false;
            for (int_fast64_t current_row = partial_start, s_index; current_row < partial_end; ++current_row)
            {
                s_index = fsc.simplices_indices[current_row];
                if (fsc.simplices[s_index].size() > 1)
                {
                    current_row_it = store_matrix.rows.emplace(std::piecewise_construct,
                                                               std::forward_as_tuple(current_row),
                                                               std::forward_as_tuple(boundary(fsc.simplices[s_index])))
                                         .first;
                    current_row_deleted = false;
                    next_pivot_it = pivots_lookup.find(max_iterator(current_row_it->second)->first);
                }
                else
                {
                    current_row_deleted = true;
                }
                while (!current_row_deleted and next_pivot_it != pivots_lookup.end())
                {
                    current_pivot = max_iterator(current_row_it->second)->first;
                    current_pivot_lookup = next_pivot_it->second;
                    current_pivot_lookup_row_it = store_matrix.rows.find(current_pivot_lookup);
                    coeff = -(current_row_it->second[current_pivot] / current_pivot_lookup_row_it->second[current_pivot]);
                    current_row_deleted = store_matrix.addscale_from_to(current_pivot_lookup_row_it, current_row_it, coeff);
                    if (!current_row_deleted)
                    {
                        next_pivot_it = pivots_lookup.find(max_iterator(current_row_it->second)->first);
                    }
                }
                if (!current_row_deleted)
                {
                    int_fast64_t start = max_iterator(current_row_it->second)->first;
                    pivots_lookup[start] = current_row;
                    proto_bars[fsc.simplices[s_index].size() - 2][start] = std::vector<int_fast64_t>{start, current_row};
                    unmatched.erase(start);
                    unmatched.erase(current_row);
                }
            }
            if (partial_end < fsc.simplices.size() - 1)
            {
                return 1;
            }
            /** Add infinite proto bars*/
            for (const auto i : unmatched)
            {
                proto_bars[fsc.simplices[fsc.simplices_indices[i]].size() - 1][i] = std::vector<int_fast64_t>{i, -1};
            }
            /** End of Add infinite proto bars*/

            /** Compute bars with filtration values*/
            for (const auto &bars_at_dimension : proto_bars)
            {
                for (const auto &bar : bars_at_dimension.second)
                {
                    if (bar.second[1] != -1)
                    {
                        if (fsc.appears_at[fsc.simplices_indices[bar.second[0]]] < fsc.appears_at[fsc.simplices_indices[bar.second[1]]])
                        {
                            bars[bars_at_dimension.first][bar.first][0] = fsc.appears_at[fsc.simplices_indices[bar.second[0]]];
                            bars[bars_at_dimension.first][bar.first][1] = fsc.appears_at[fsc.simplices_indices[bar.second[1]]];
                        }
                    }
                    else
                    {
                        bars[bars_at_dimension.first][bar.first][0] = fsc.appears_at[fsc.simplices_indices[bar.second[0]]];
                        bars[bars_at_dimension.first][bar.first][1] = std::numeric_limits<float>::infinity();
                    }
                }
            }
            /** End of Compute bars with filtration values*/
            return 0;
        }
    };

}
#endif /* HOMOLOGYCOMPUTER_HPP*/
