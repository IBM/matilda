#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP
#include <tuple>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <climits>
#include <cassert>
#include <map>
#include <unordered_map>
#include "RowVector.hpp"

namespace matilda
{

    template <typename Row>
    typename std::unordered_map<int_fast64_t, Row>::iterator max_iterator(std::unordered_map<int_fast64_t, Row> &x)
    {
        if (x.empty())
        {
            return x.end();
        }
        auto it = x.begin();
        typename std::unordered_map<int_fast64_t, Row>::iterator result;
        int_fast64_t max = -1;
        for (; it != x.end(); ++it)
        {
            if (it->first > max)
            {
                result = it;
                max = it->first;
            }
        }
        return result;
    }

    template <typename Coefficient>
    typename std::map<int_fast64_t,
                      Coefficient>::iterator
    max_iterator(typename std::map<int_fast64_t,
                                   Coefficient> &x)
    {

        return --x.end();
    }

    /**
    Class to model sparse column matrices in column compressed form.
    */

    template <typename Coefficient, typename Row>
    class SparseMatrix
    {
    public:
        int_fast64_t shape[2];
        SparseMatrix() : shape{0, 0} {}
        std::unordered_map<int_fast64_t, Row> rows;
        std::map<int_fast64_t, int_fast64_t> pivots;
        void update_shape()
        {
            if (rows.empty())
            {
                shape[0] = shape[1] = 0;
                return;
            }
            shape[0] = max_iterator(rows)->first + 1;
            int_fast64_t max_column = -1;
            int_fast64_t current_max_column;
            for (const auto &it : rows)
            {
                current_max_column = max_iterator(it.second)->first + 1;
                if (max_column < current_max_column)
                {
                    max_column = current_max_column;
                }
            }
            shape[1] = max_column;
        }
        void print()
        {
            if (!rows.empty())
            {
                update_shape();
                for (int i = 0; i < shape[0]; ++i)
                {
                    for (int j = 0; j < shape[1]; ++j)
                    {
                        const auto row_i = rows.find(i);
                        if (row_i != rows.end())
                        {
                            const auto entry = row_i->second.find(j);
                            if (entry != row_i->second.end())
                            {
                                std::cout << entry->second.value;
                            }
                            else
                            {
                                std::cout << "0";
                            }
                        }
                        else
                        {
                            std::cout << "0";
                        }
                    }
                    std::cout << "\n";
                }
            }
        }
        /** IMPORTANT: */
        /**vvv Assumes non-zero coefficient , also we assume that destination row always exists vvv*/
        bool addscale_from_to(int_fast64_t from_index, int_fast64_t to_index, Coefficient coeff)
        {
            bool result = false;
            auto dest_it = rows.find(to_index);
            Row *destination = &(dest_it->second);
            for (auto &it : rows[from_index])
            {
                auto it_in_destination = destination->find(it.first);
                if (it_in_destination != destination->end())
                {
                    it_in_destination->second.mult_add(coeff.value, it.second.value);
                    if (it_in_destination->second == 0)
                    {
                        destination->erase(it_in_destination);
                    }
                }
                else
                {
                    destination->emplace(it.first, coeff * it.second);
                }
            }
            if (destination->empty())
            {
                rows.erase(dest_it);
                result = true;
            }
            return result;
        }

        /** IMPORTANT: */
        /**vvv Assumes non-zero coefficient , also we assume that destination row always exists vvv*/
        bool addscale_from_to(typename std::unordered_map<int_fast64_t, Row>::iterator &from_index,
                              typename std::unordered_map<int_fast64_t, Row>::iterator &to_index,
                              Coefficient coeff)
        {
            bool result = false;
            for (auto &it : from_index->second)
            {
                auto it_in_destination = to_index->second.find(it.first);
                if (it_in_destination != to_index->second.end())
                {
                    it_in_destination->second.mult_add(coeff.value, it.second.value);
                    if (it_in_destination->second == 0)
                    {
                        to_index->second.erase(it_in_destination);
                    }
                }
                else
                {
                    to_index->second.emplace(it.first, coeff * it.second);
                }
            }
            if (to_index->second.empty())
            {
                rows.erase(to_index);
                result = true;
            }
            return result;
        }
    };
}
#endif /* SPARSEMATRIX_HPP */
