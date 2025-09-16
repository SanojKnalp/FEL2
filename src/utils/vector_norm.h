#ifndef MAX_NORM_VEC
#define  MAX_NORM_VEC

#include <iostream>
#include <cmath>
#include <stdexcept>

/**
 * @brief A utility class providing norm computations for vector-like types.
 *
 * This class template implements methods for computing different norms
 * of vectors. Currently, only the maximum norm is implemented.
 *
 * @tparam Vector_type A container type supporting:
 *   - operator[] for element access
 *   - size() method
 *   - empty() method
 */
template <typename Vector_type>
class Norms {
public:
    /**
     * @brief Compute the maximum norm of a vector.
     *
     * @param x The input vector.
     * @return The maximum absolute value among the vector's elements.
     *
     * @throws std::runtime_error if the input vector is empty.
     */
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