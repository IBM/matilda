#ifndef ROWVECTOR_HPP
#define ROWVECTOR_HPP
#include <tuple>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <climits>
#include <cassert>
#include <map>
#include <unordered_map>

namespace matilda{
template <typename Coefficient>
class Pair
{
public:
    Pair() {}
    Pair(int_fast64_t a, Coefficient b): first(a), second(b) {}
    int_fast64_t first;
    Coefficient second;
};


template <typename Coefficient>
class RowVector
{
public:
    std::vector<std::pair<int_fast64_t,Coefficient>> entries;
    RowVector() {}
    RowVector(int_fast64_t a, Coefficient b): entries{{a,b}} {}
    typedef typename std::vector<std::pair<int_fast64_t,Coefficient>>::iterator iterator;
    typedef typename std::vector<std::pair<int_fast64_t,Coefficient>>::const_iterator const_iterator;
    iterator begin()
    {
        return entries.begin();
    }
    iterator end()
    {
        return entries.end();

    }
    const_iterator cbegin() const
    {
        return entries.cbegin();
    }
    const_iterator cend() const
    {
        return entries.cend();

    }
    iterator find(int_fast64_t x) /** Assumes non empty entries vector */
    {
        return std::find_if(entries.begin(), entries.end(), [&]
                         (std::pair<int_fast64_t,Coefficient> &a){ return a.first == x;} );

    } // Experiments with binary search + sorting when emplacing were much slower than with linear search

    iterator erase(int_fast64_t i)
    {
        iterator it;
        for( it = entries.begin(); it!=entries.end(); ++it)
        {
            if (it->first == i)
            {
                return entries.erase(it);
            }
        }
        return entries.end();
    }
    iterator erase(iterator it)
    {
        return entries.erase(it);
    }

    void emplace(int_fast64_t key, Coefficient value)
    {
        entries.emplace_back(key,value);
    }

    size_t size()
    {
        return entries.size();
    }
    
    bool empty()
    {
        return entries.empty();
    }

    bool empty() const
    {
        return entries.empty();
    }
    Coefficient operator[](int_fast64_t key)
    {
        auto it = find(key);
        return it->second;
    }
};

template <typename Coefficient>
typename RowVector<Coefficient>::const_iterator max_iterator(const RowVector<Coefficient> &x)
/**
 * @brief Returns an iterator pointing to the element with maximum column position in x. In short,
 *         the non-zero element furthest to the right in x.
 * 
 */
{
    if (x.empty())
    {
        return x.cend();
    }
    auto it = x.cbegin();
    typename RowVector<Coefficient>::const_iterator result;
    int_fast64_t max=-1;
    for( ; it!=x.cend(); ++it)
    {
        if (it->first > max)
        {
            result = it;
            max = it->first;
        }
    }
    return result;
}



}
#endif // ROWVECTOR_HPP
