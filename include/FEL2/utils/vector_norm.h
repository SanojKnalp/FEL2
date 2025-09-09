#ifndef MAX_NORM_VEC
#define  MAX_NORM_VEC

#include <iostream>
#include <cmath>
#include <stdexcept>


template <typename Vector_type>
class Norms {
public:
    double compute_max_norm(const Vector_type& x);
};

template <typename Vector_type>
double Norms<Vector_type>::compute_max_norm(const Vector_type &x)
{
    if (x.empty()) { throw std::runtime_error("compute_max_norm: Cannot compute maximum norm of empty vector."); }

    double highest_value = std::abs(x[0]);
    for (std::size_t i = 1; i < x.size(); i++) {
        highest_value = std::max(highest_value, std::abs(x[i]));
    }
    return highest_value;
}
#endif 
