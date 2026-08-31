#pragma once
#ifndef COME_LINALG_FULLMATRIX_HPP
#define COME_LINALG_FULLMATRIX_HPP

#include <armadillo>
#include <stdexcept>
#include "COME_Linalg_vector.hpp"

namespace LinearAlgebra {

    template <typename T>
    class FullMatrix {
    public:
        arma::Mat<T> mat;

        // Constructors
        FullMatrix() = default;
        FullMatrix(size_t rows, size_t cols) : mat(rows, cols) {}
        FullMatrix(const arma::SpMat<T>& A) : mat(A) {}

        // Shape
        size_t rows() const { return mat.n_rows; }
        size_t cols() const { return mat.n_cols; }

        // Element access
        const T& operator()(size_t i, size_t j) const {
            if (i >= rows() || j >= cols())
                throw std::out_of_range("FullMatrix index out of bounds");
            return mat(i, j);
        }
        //Write access
        auto operator()(size_t i, size_t j)-> decltype(mat(i, j)) {
            if (i >= rows() || j >= cols())
                throw std::out_of_range("FullMatrix index out of bounds");
            return mat(i, j);
        }


        // Set an element (inserts/updates value). Use assignment into SpMat.
        // Note: mat(i,j) = value will insert the value for Full matrices.
        void set(size_t i, size_t j, const T& value) {
            if (i >= rows() || j >= cols())
                throw std::out_of_range("FullMatrix index out of bounds");
            mat(i, j) = value;
        }

        // reinits the matrix with the proper size
        void reinit(size_t size) {
            mat.zeros(size, size);
        }

        // Frobenius norm
        T norm() const {
            return arma::norm(mat, "fro");
        }

        // Matrix-vector multiplication: y = A*x
        Vector<T> operator*(const Vector<T>& x) const {
            if (cols() != x.size())
                throw std::invalid_argument("Dimension mismatch: A*x");

            Vector<T> out(rows());
            out.col = mat * x.col;
            return out;
        }

        // Matrix-matrix multiplication: C = A*B
        FullMatrix<T> operator*(const FullMatrix<T>& B) const {
            if (cols() != B.rows())
                throw std::invalid_argument("Dimension mismatch: A*B");

            FullMatrix<T> C(rows(), B.cols());
            C.mat = mat * B.mat;
            return C;
        }

        // Scalar multiplication
        FullMatrix<T> operator*(T scalar) const {
            FullMatrix<T> out(rows(), cols());
            out.mat = scalar * mat;
            return out;
        }

        //overload operator 
        FullMatrix<T>& operator*=(T scalar) {
            mat *= scalar;
            return *this;
        }

        // Addition
        FullMatrix<T> operator+(const FullMatrix<T>& B) const {
            if (rows() != B.rows() || cols() != B.cols())
                throw std::invalid_argument("Dimension mismatch: A+B");

            FullMatrix<T> out(rows(), cols());
            out.mat = mat + B.mat;
            return out;
        }

        // Subtraction
        FullMatrix<T> operator-(const FullMatrix<T>& B) const {
            if (rows() != B.rows() || cols() != B.cols())
                throw std::invalid_argument("Dimension mismatch: A-B");

            FullMatrix<T> out(rows(), cols());
            out.mat = mat - B.mat;
            return out;
        }
        //Extract the diagonal Matrix
        // No pass by reference cause we create something new
        FullMatrix<T> diagonalMatrix() const
        {
            FullMatrix<T> C;
            C.mat = arma::diagmat(mat);
            return C;
        }

        //function that extract the diagonal and inverts the diagonal. Not super happy, but ok.
        FullMatrix<T> inverseDiagonalMatrix() const
        {
            FullMatrix<T> C(rows(), cols());

            for (size_t i = 0; i < rows(); ++i) {
                T d = mat(i, i);
                if (d == T(0))
                    throw std::runtime_error("Zero on diagonal in Jacobi");
                C.mat(i, i) = T(1) / d;
            }
            return C;
        }

        //Eliminates row and column from matrix. Important for BCs.
        FullMatrix<T> eliminateRowAndColumn(unsigned int index)
        {
            this->mat.shed_col(index);
            this->mat.shed_row(index);
            return *this;
        }
        //analytical formula for fast 2x2 inverse. SHOULD ONLY BE USED FOR 2x2 and really only if one is 100% sure it's a 2x2 matrix.
        FullMatrix<T>& inverse2x2()
        {

            auto detA = arma::det(mat);
            if (detA == T(0))
            {
                throw std::runtime_error("inverse2x2(): singular matrix");
            }

            T a = mat(0, 0);
            T b = mat(0, 1);
            T c = mat(1, 0);
            T d = mat(1, 1);

            mat(0, 0) = d;
            mat(0, 1) = -b;
            mat(1, 0) = -c;
            mat(1, 1) = a;

            mat /= detA;
            return *this;
        }

        // analytical formula for fast 3x3 inverse.SHOULD ONLY BE USED FOR 3x3 and really only if one is 100 % sure it's a 3x3 matrix.
        FullMatrix<T>& inverse3x3()
        {
            auto detA = arma::det(mat);
            if (detA == T(0))
            {
                throw std::runtime_error("inverse3x3(): singular matrix");
            }
            T a = mat(0, 0);
            T b = mat(0, 1);
            T c = mat(0, 2);
            T d = mat(1, 0);
            T e = mat(1, 1);
            T f = mat(1, 2);
            T g = mat(2, 0);
            T h = mat(2, 1);
            T i = mat(2, 2);

            mat(0, 0) = e * i - f * h;
            mat(0, 1) = c * h - b * i;
            mat(0, 2) = b * f - c * e;
            mat(1, 0) = f * g - d * i;
            mat(1, 1) = a * i - c * g;
            mat(1, 2) = c * d - a * f;
            mat(2, 0) = d * h - e * g;
            mat(2, 1) = b * g - a * h;
            mat(2, 2) = a * e - b * d;

            mat /= detA;
            return *this;
        }

        T det()
        {
            return arma::det(mat);

        }

        FullMatrix<T> inverse1x1_new()
        {
            if (mat(0, 0) == 0)
            {
                throw std::runtime_error("inverse1x1(): singular matrix");
            }
            FullMatrix<T> C(1, 1);
            C(0, 0) = 1 / mat(0, 0);
            return C;
        }

        FullMatrix<T> inverse2x2_new()
        {
            auto detA = arma::det(mat);
            if (detA == T(0))
            {
                throw std::runtime_error("inverse2x2(): singular matrix");
            }
            FullMatrix<T> C(2, 2);
            T a = mat(0, 0);
            T b = mat(0, 1);
            T c = mat(1, 0);
            T d = mat(1, 1);

            C(0, 0) = d;
            C(0, 1) = -b;
            C(1, 0) = -c;
            C(1, 1) = a;

            C /= detA;
            return C;
        }

        FullMatrix<T> inverse3x3_new()
        {
            auto detA = arma::det(mat);
            if (detA == T(0))
            {
                throw std::runtime_error("inverse3x3(): singular matrix");
            }
            FullMatrix<T> C(3, 3);
            T a = mat(0, 0);
            T b = mat(0, 1);
            T c = mat(0, 2);
            T d = mat(1, 0);
            T e = mat(1, 1);
            T f = mat(1, 2);
            T g = mat(2, 0);
            T h = mat(2, 1);
            T i = mat(2, 2);

            C(0, 0) = e * i - f * h;
            C(0, 1) = c * h - b * i;
            C(0, 2) = b * f - c * e;
            C(1, 0) = f * g - d * i;
            C(1, 1) = a * i - c * g;
            C(1, 2) = c * d - a * f;
            C(2, 0) = d * h - e * g;
            C(2, 1) = b * g - a * h;
            C(2, 2) = a * e - b * d;

            C /= detA;
            return C;

        }

        FullMatrix<T> small_dim_inverse()
        {
            if constexpr (this->cols() == 1)
            {
                return this->inverse1x1_new();
            }
            else if constexpr (this->cols() == 2)
            {
                return this->inverse2x2_new();
            }
            else if constexpr (this->cols() == 3)
            {
                return this->inverse3x3_new();
            }
        }

        FullMatrix<T>& transpose()
        {
            return this->mat.t();
        }


    };

    template <typename T>
    FullMatrix<T> operator*(T s, const FullMatrix<T>& A) {
        return A * s;
    }



} // namespace LinearAlgebra

#endif // COME_LINALG_FULLMATRIX_HPP
