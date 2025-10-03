#include "../src/LinearAlgebra/sparse/FEL2_linalg_sparsematrix.hpp"
#include <armadillo>

int main() {
    using namespace LinearAlgebra;
    SparseMatrix<double> A(5,5);
    arma::vec x = arma::randu<arma::vec>(5);
    arma::vec b = arma::randu<arma::vec>(5);
    arma::vec r;

    A.residual(r, x, b);
}