#ifndef REAL_HPP
#define REAL_HPP
#include <iostream>
#include <cstdlib>
#include <map>
#include <math.h>
#include "globals.hpp"
/** Characteristic zero coefficient class. For use as Coefficient in templates. **/


namespace matilda
{

class Real{
public:
    
    float value;
    Real(){};
    explicit Real(float x):value(x)  {};
    Real inverse() const;
    Real operator-()
    {
        return Real(-this->value);
    }
    Real operator*=(const Real & b)
    {
        this->value*=b.value;
        return *this;
    }
    friend bool operator==(Real  a, const Real & b)
    {
        return a.value == b.value;
    }
    friend bool operator==(Real  a, const int_fast64_t b)
    {
        return ::fabsf(a.value - float(b)) < 1e-8;
    }
    
    friend bool operator!=(Real  a, const Real & b)
    {
        return ::fabsf(a.value - b.value)>1e-8;
    }
    friend Real operator*(Real  a, const Real & b)
    {
        a *= b;
        return a;
    }
    Real operator+=(const Real & b)
    {
        this->value+=b.value;
        return *this;
    }
    friend Real operator+(Real  a, const Real & b)
    {
        a += b;
        return a;
    }
    Real operator/=(const Real & b) /// Assumes b!=0
    {
        this->value/=b.value;
        return *this;
    }
    friend Real operator/(Real  a, const Real & b)
    {
        a /= b;
        return a;
    }
    void mult_add(float b, float c )
    {
        this->value+= b*c;
    }
    
    friend std::ostream& operator<<(std::ostream& out, const Real & x)
    {
        out<<x.value;
        return out;
    }
};
}
#endif /* REAL_HPP */
