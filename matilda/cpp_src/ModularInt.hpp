#ifndef MODULARINT_HPP
#define MODULARINT_HPP
#include <iostream>
#include <cstdlib>
#include <map>
#include "globals.hpp"
/** Finite characteristic coefficient class. For use as Coefficient in templates. **/


namespace matilda
{

inline
int_fast64_t mod_inverse(int_fast64_t a) // Assumes a non zero
{
    if (modulus == 2)
    {
        return 1;
    }
    a = (a%modulus + modulus)%modulus;
    int_fast64_t old_r = modulus;
    int_fast64_t r = a;
    int_fast64_t old_t=0;
    int_fast64_t t = 1;
    int_fast64_t q, temp_1, temp_2;
    while (r!=0)
    {
        q = old_r/r;
        temp_1 = old_t;
        temp_2 = old_r;
        old_t = t;
        t = temp_1 - q*t;
        old_r = r;
        r = temp_2 - q*r;
    }
    return old_t;
}


class ModularInt{
public:
    
    int_fast64_t value;
    ModularInt(){};
    explicit ModularInt(int_fast64_t x):value(modulus_lookup[(modulus_lookup[x+MOD_LOOKUP_SIZE/2]
                                               + modulus) + modulus + MOD_LOOKUP_SIZE/2])  {};
    ModularInt inverse() const;
    ModularInt operator-()
    {
        return ModularInt(-this->value);
    }
    ModularInt operator*=(const ModularInt & b)
    {
        this->value= modulus_lookup[(this->value * b.value) + MOD_LOOKUP_SIZE/2];
        return *this;
    }
    friend bool operator==(ModularInt  a, const ModularInt & b)
    {
        return a.value == b.value;
    }
    friend bool operator==(ModularInt  a, const int_fast64_t b)
    {
        return a.value == b;
    }
    
    friend bool operator!=(ModularInt  a, const ModularInt & b)
    {
        return a.value != b.value;
    }
    friend ModularInt operator*(ModularInt  a, const ModularInt & b)
    {
        a *= b;
        return a;
    }
    ModularInt operator+=(const ModularInt & b)
    {
        this->value= modulus_lookup[(this->value + b.value)+MOD_LOOKUP_SIZE/2];
        return *this;
    }
    friend ModularInt operator+(ModularInt  a, const ModularInt & b)
    {
        a += b;
        return a;
    }
    ModularInt operator/=(const ModularInt & b) /// Assumes b!=0
    {
        if (modulus == 2)
        {
            this->value = 1;
            return *this;
        }
        this ->value = modulus_lookup[modulus_lookup[this->value * modulus_inverse_lookup[b.value + MOD_INV_LOOKUP_SIZE/2] + MOD_LOOKUP_SIZE/2 ] + modulus + MOD_LOOKUP_SIZE/2];
        return *this;
    }
    friend ModularInt operator/(ModularInt  a, const ModularInt & b)
    {
        a /= b;
        return a;
    }
    void mult_add(int_fast64_t b, int_fast64_t c )
    {
        this->value+= b*c;
        this->value = modulus_lookup[this->value + MOD_LOOKUP_SIZE/2];
    }
    
    friend std::ostream& operator<<(std::ostream& out, const ModularInt & x)
    {
        out<<x.value;
        return out;
    }
};
}
#endif /* MODULARINT_HPP */
