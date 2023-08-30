#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <numeric>
#include <iostream>
/* Matrix is meant only to be read, not written to (its destructor
 does NOT free memory nor its constructor reserves memory) */
class Matrix
{
public:

    float *entries;
    int_fast64_t shape[2];
    Matrix(){};
    Matrix(float *entries, int_fast64_t n_rows, int_fast64_t n_columns)
    {
        this->entries = entries;
        shape[0] = n_rows;
        shape[1] = n_columns;
    }
    inline
    const float operator() (int_fast64_t i, int_fast64_t j) const
    {
        return entries[shape[1]*i + j];
    }
    float & operator() (int_fast64_t i, int_fast64_t j)
    {
        return entries[shape[1]*i + j];
    }
    ~Matrix() 
    {
        return;
    }
    void print()
    {
        std::cout.precision(3);
        std::cout<<std::fixed;
        for (int i=0; i < shape[0]; ++i)
        {
            for (int j=0; j < shape[1]; ++j)
            {
                std::cout << this->operator()(i, j)<<" ";
            }
            std::cout << "\n";
        }
    }
};

#endif /* MATRIX_HPP */
