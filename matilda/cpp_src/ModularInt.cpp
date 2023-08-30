#include "ModularInt.hpp"

namespace matilda
{
int_fast64_t modulus=3;
int_fast64_t modulus_lookup[MOD_LOOKUP_SIZE];
int_fast64_t modulus_inverse_lookup[MOD_INV_LOOKUP_SIZE];


void set_modulus(int_fast64_t new_modulus)
{
    if (new_modulus == 0)
    {
        return;
    }
    modulus = new_modulus;
    for(int_fast64_t i = -MOD_LOOKUP_SIZE/2; i<MOD_LOOKUP_SIZE/2; ++i)
    {
        modulus_lookup[i + MOD_LOOKUP_SIZE/2] = (i%modulus + modulus)%modulus;
    }
    for(int_fast64_t i = -MOD_INV_LOOKUP_SIZE/2; i<MOD_INV_LOOKUP_SIZE/2; ++i)
    {
        if (i!=0)
        {
            modulus_inverse_lookup[i + MOD_INV_LOOKUP_SIZE/2] = (mod_inverse(i)%modulus + modulus)%modulus;
        }
    }
}
}
