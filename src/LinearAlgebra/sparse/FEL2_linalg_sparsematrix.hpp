#ifndef FEL2_LINALG_SPARSEMATRIX_HPP
#define FEL2_LINALG_SPARSEMATRIX_HPP

#include <armadillo>

namespace LinearAlgebra
{
    template <typename T>
    class SparseMatrix {
        public:
            arma::SpMat<T> mat;
            SparseMatrix() = default;
            SparseMatrix(size_t rows, size_t cols) : mat(rows, cols){};





            //Compute the residual
            // Residual: r = b - A * x
            void residual(arma::Col<T>& r, const arma::Col<T>& x, const arma::Col<T>& b) const {
                if (mat.n_cols != x.n_rows || mat.n_rows != b.n_rows) {
                    throw std::invalid_argument("Dimension mismatch in residual()");
                }
                r = b - mat * x;
            }

            //Compute the Frobenius Norm
            T frobenius_norm() const {
                return arma::norm(mat,"fro");
            }

            void mmult(SparseMatrix<T>& C, const SparseMatrix<T>& B) const {
                if (mat.n_cols != B.mat.n_rows) {
                    throw std::invalid_argument("Dimension mismatch in mmult()");
                }
                C.mat = mat * B.mat;
            }

            void vmult(arma::Col<T>& c, const arma::Col<T>& x) const {
                if (mat.n_cols != B.mat.n_rows) {
                    throw std::invalid_argument("Dimension mismatch in mmult()");
                }
                c = mat * x;
            }

            void smult(SparseMatrix<T> &C, T scalar) const {
                if (mat.n_cols != B.mat.n_rows) {
                    throw std::invalid_argument("Dimension mismatch in mmult()");
                }
                C = scalar*mat;
            }
             






        };
}

#endif