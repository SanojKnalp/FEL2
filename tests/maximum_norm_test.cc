#include <vector>
#include <iostream>
#include  <complex>
#include "../src/utils/vector_norm.h"


int main(){
    std::vector<int> v1 = {10, 20, 35, 30};
    std::vector<float> v2 = {10.0, 35.0, 20.0, 8.0};
    std::vector<std::complex<double>> v3 = {{2.0, 1.0}, {0.0, 4.0}, {1.0, 5.0}, {3.0, 0.0}};
    Norms<std::vector<int>> v1_max;
    double norm = v1_max.compute_max_norm(v1);
    std::cout << "Norm1: " << norm << std::endl;
    Norms<std::vector<float>> v2_max;
    double norm2 = v2_max.compute_max_norm(v2);
    std::cout << "Norm2: " << norm2 << std::endl;
    Norms<std::vector<std::complex<double>>> v3_max;
    double norm3 = v3_max.compute_max_norm(v3);
    std::cout << "Norm3: " << norm3 << std::endl;

}
