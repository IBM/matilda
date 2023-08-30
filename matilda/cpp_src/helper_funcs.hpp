#ifndef HELPER_FUNCS_HPP
#define HELPER_FUNCS_HPP
#include <stdint.h>


int_fast64_t binomial_coeff(int_fast64_t n, int_fast64_t k)
{
    int_fast64_t result=1;
    if (k>0)
    {   
        return (n * binomial_coeff(n-1,k-1))/k;
    }
    if (k==0)
    {
        return 1;
    }
    return result;
}

int_fast64_t are_epsilon_close(const Matrix &m, int_fast64_t i, int_fast64_t j, float r)
{
    return m(i,j) < r ? 0:1;
}

#endif /* HELPER_FUNCS_HPP */
