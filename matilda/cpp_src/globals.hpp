#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#define MOD_LOOKUP_SIZE 1024
#define MOD_INV_LOOKUP_SIZE 1024

#include <cstdlib>
#include <map>

namespace matilda
{
extern int_fast64_t modulus;
extern int_fast64_t modulus_lookup[MOD_LOOKUP_SIZE];
extern int_fast64_t modulus_inverse_lookup[MOD_INV_LOOKUP_SIZE];
void set_modulus(int_fast64_t new_modulus);
}

#endif /* GLOBALS_HPP */