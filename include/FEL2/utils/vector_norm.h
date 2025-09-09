#ifndef MAX_NORM_VEC
#define  MAX_NORM_VEC

#include <iostream>
#include <cmath>


template <typename Vector_type>
class Norms {
public:
    double compute_max_norm(const Vector_type& x);
};


template <typename Vector_type>
double Norms<Vector_type>::compute_max_norm(const Vector_type &x)
{
    int size = x.size();
    double highest_value = std::abs(x[0]) ;
    for(int i = 0; i < size; i++){
        if(abs(x[i])>highest_value){
            highest_value = std::abs(x[i]);
        }
    }

    return highest_value;
}
#endif 
