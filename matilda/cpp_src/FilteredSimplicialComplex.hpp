#ifndef FILTEREDSIMPLICIALCOMPLEX_HPP
#define FILTEREDSIMPLICIALCOMPLEX_HPP

#include <stdint.h>
#include <vector>
#include <climits>
#include <set>
#include "Matrix.hpp"

namespace matilda
{

class FilteredSimplicialComplex
{
public:
    using Simplex = std::vector<int_fast64_t> ;
    std::vector<Simplex> simplices;
    std::vector<int_fast64_t> simplices_indices;
    int_fast64_t dimension;
    std::vector<float> appears_at;
    FilteredSimplicialComplex() : dimension(0),diameter_helper(0.) {}
    void construct_vietoris_from_metric(const Matrix & matrix,
                                        int_fast64_t dimension,
                                        float diameter);
    /** Helpers for VR complexes, for now in FSC class but might change to special VR class in future */
    Matrix matrix_helper;
    float diameter_helper;
    std::vector<std::vector<int_fast64_t>> neighbors_helper;
    void sort_indices();
    void add_cofaces(const std::vector<int_fast64_t> &simplex,
                     const std::vector<int_fast64_t> & neighbors);
    void compute_filter_function();
    void compute_all_neighbors();
    float max_filter(const std::vector<int_fast64_t> &simplex);
    std::vector<int_fast64_t> compute_right_neighbors(int_fast64_t i);
    /** End of helpers for VR complexes */
};
}

#endif /* FILTEREDSIMPLICIALCOMPLEX_HPP */
