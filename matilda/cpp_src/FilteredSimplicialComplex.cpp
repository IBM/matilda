#include "FilteredSimplicialComplex.hpp"
#include "set"
#include <numeric>
#include <iostream>
#include <algorithm>


namespace matilda{
void FilteredSimplicialComplex::add_cofaces(const std::vector<int_fast64_t> &simplex,
                                            const std::vector<int_fast64_t> & neighbors)
{
    simplices.push_back(simplex);
    if (simplex.size()<dimension+1)
    {
        for(auto i: neighbors)
        {   
            std::vector<int_fast64_t> new_simplex=simplex;
            std::vector<int_fast64_t> new_neighbors;
            std::set_intersection(neighbors.begin(),
                                    neighbors.end(),
                                    neighbors_helper[i].begin(), 
                                    neighbors_helper[i].end(),
                                    std::back_inserter(new_neighbors) ) ;
            new_simplex.push_back(i);
            add_cofaces(new_simplex, new_neighbors);
        }
    }
};

void FilteredSimplicialComplex::sort_indices()
{
    simplices_indices.resize(simplices.size());
    std::iota(simplices_indices.begin(), simplices_indices.end(),0);
    std::stable_sort(simplices_indices.begin(),
                     simplices_indices.end(),
                     [this](int_fast64_t i, int_fast64_t j)
                     {return appears_at[i]<appears_at[j];});
}

float FilteredSimplicialComplex::max_filter(const std::vector<int_fast64_t> &simplex)
{

        if (simplex.size() == 1)
        {
            return 0;
        }
        if (simplex.size() == 2)
        {
            return matrix_helper(simplex[0],simplex[1]);
        }
        float result = 0;
        for (int_fast64_t i=0; i<simplex.size(); ++i)
        {
            std::vector<int_fast64_t> new_simplex=simplex;
            new_simplex.erase(new_simplex.begin()+i);
            if (result<max_filter(new_simplex))
            {
                result = max_filter(new_simplex);
            }
        }
        return result;
}

void FilteredSimplicialComplex::compute_filter_function()
{
    for (int_fast64_t i=0; i<simplices.size(); ++i)
    {
        if (simplices[i].size() == 1)
        {
            appears_at.push_back(0);
        }
        else if (simplices[i].size() == 2)
        {
            appears_at.push_back(matrix_helper(simplices[i][0],simplices[i][1]));
        }
        else 
        {
            appears_at.push_back( max_filter(simplices[i]));
        }
    }
}

std::vector<int_fast64_t> FilteredSimplicialComplex::compute_right_neighbors( int_fast64_t i)
{
    std::vector<int_fast64_t>  result;
    for(int_fast64_t j=i+1; j<matrix_helper.shape[0]; ++j)
    {
        if(matrix_helper(i,j)<diameter_helper)
        {
            result.push_back(j);
        }
    }
    return result;
}


void FilteredSimplicialComplex::compute_all_neighbors()
{
    for (int_fast64_t i = 0; i<matrix_helper.shape[0]; ++i)
    {
        std::vector<int_fast64_t>  current;
        for(int_fast64_t j=0; j<i; ++j)
        {
            if(matrix_helper(i,j)<diameter_helper)
            {
                current.push_back(j);
            }
        }
        neighbors_helper.push_back(current);
    }
}


void \
FilteredSimplicialComplex::\
construct_vietoris_from_metric(const Matrix &matrix,  
                                            int_fast64_t dimension,
                                            float diameter) 
{
    int_fast64_t matrix_size[2];
    matrix_size[0] = matrix.shape[0];
    matrix_size[1] = matrix.shape[1];
    this->dimension = dimension;
    this->matrix_helper = matrix;
    this->diameter_helper = diameter;
    compute_all_neighbors();
    for(int_fast64_t i=0; i<matrix.shape[0];++i )
    {
        add_cofaces(Simplex{i},neighbors_helper[i]);
    }
    for (int_fast64_t i=0; i<simplices.size(); ++i)
    {
        std::reverse(simplices[i].begin(),simplices[i].end());
    }
    compute_filter_function();
    sort_indices();
}
}
