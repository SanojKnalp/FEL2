#include "../catch.hpp"
#include "LinearAlgebra/COME_Linalg_dense.hpp"
#include "LinearAlgebra/COME_Linalg_vector.hpp"
#include <random>
#include <cmath>

// Helper: fill a Vector<double> with random values
void fill_random_dense(LinearAlgebra::Vector<double>& v) {
    static std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (size_t i = 0; i < v.size(); ++i) {
        v(i) = dist(rng);
    }
}

TEST_CASE("Full matrix residual test", "[Full][residual]") {
    using namespace LinearAlgebra;

    FullMatrix<double> A(5, 5);
    Vector<double> x(5);
    Vector<double> b(5);

    fill_random_dense(x);
    fill_random_dense(b);

    Vector<double> r = b - (A * x);

    // For a zero-initialized FullMatrix, b - A*x should equal b
    for (size_t i = 0; i < r.size(); ++i) {
        REQUIRE(r(i) == Approx(b(i)).epsilon(1e-12));
    }
}

TEST_CASE("Full matrix Frobenius norm", "[Full][norm]") {
    using namespace LinearAlgebra;

    FullMatrix<double> A(5, 5);

    A.set(0, 0, 5.0);
    A.set(1, 1, 1.0);
    A.set(2, 2, 2.0);
    A.set(3, 3, 3.0);
    A.set(4, 4, 4.0);

    double frob = A.norm();

    REQUIRE(frob == Approx(std::sqrt(55.0)).epsilon(1e-12));
}

TEST_CASE("Full matrix scalar multiplication", "[Full][scalarmult]") {
    using namespace LinearAlgebra;

    FullMatrix<double> A(5, 5);

    A(0, 0) = 5.0;
    A(1, 1) = 1.0;
    A(2, 2) = 2.0;
    A(3, 3) = 3.0;
    A(4, 4) = 4.0;


    double factor = 2.0;
    FullMatrix<double> B = A * factor;

    REQUIRE(B(0, 0) == Approx(10.0));
    REQUIRE(B(1, 1) == Approx(2.0));
    REQUIRE(B(2, 2) == Approx(4.0));
    REQUIRE(B(3, 3) == Approx(6.0));
    REQUIRE(B(4, 4) == Approx(8.0));
}

TEST_CASE("Full matrix matrix-vector product")
{
    using namespace LinearAlgebra;

    FullMatrix<float> A(2, 2);
    Vector<float> b(2);
    A(0, 0) = 1.0;
    A(1, 1) = 2.0;
    b(1) = 1.0;
    b(0) = 1.0;
    Vector<float> c = A * b;
    REQUIRE(c(0) == Approx(1.0));
    REQUIRE(c(1) == Approx(2.0));

}

TEST_CASE("Full matrix matrix-matrix product")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(5, 5);
    FullMatrix<double> B(5, 5);

    for (int i = 0; i < A.rows(); i++)
    {
        A(i, i) = 1.0;
        B(i, i) = i * 1.0;
    }
    FullMatrix<double> C = A * B;
    REQUIRE(C(0, 0) == Approx(0.0));
    REQUIRE(C(1, 1) == Approx(1.0));
    REQUIRE(C(2, 2) == Approx(2.0));
    REQUIRE(C(3, 3) == Approx(3.0));
    REQUIRE(C(4, 4) == Approx(4.0));
}

TEST_CASE("Full matrix matrix addition")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(5, 5);
    FullMatrix<double> B(5, 5);

    for (int i = 0; i < A.rows(); i++)
    {
        A(i, i) = 1.0;
        B(i, i) = i * 1.0;
    }
    FullMatrix<double> C = A + B;
    REQUIRE(C(0, 0) == Approx(1.0));
    REQUIRE(C(1, 1) == Approx(2.0));
    REQUIRE(C(2, 2) == Approx(3.0));
    REQUIRE(C(3, 3) == Approx(4.0));
    REQUIRE(C(4, 4) == Approx(5.0));
}

TEST_CASE("Full matrix matrix subtraction")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(5, 5);
    FullMatrix<double> B(5, 5);

    for (int i = 0; i < A.rows(); i++)
    {
        A(i, i) = 1.0;
        B(i, i) = i * 1.0;
    }
    FullMatrix<double> C = A - B;
    REQUIRE(C(0, 0) == Approx(1.0));
    REQUIRE(C(1, 1) == Approx(0.0));
    REQUIRE(C(2, 2) == Approx(-1.0));
    REQUIRE(C(3, 3) == Approx(-2.0));
    REQUIRE(C(4, 4) == Approx(-3.0));
}

TEST_CASE("Full matrix extract diagonal")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(5, 5);
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            A(i, j) = i * j * 1.0;
        }
    }
    FullMatrix<double> C = A.diagonalMatrix();

    REQUIRE(C(0, 0) == Approx(0.0));
    REQUIRE(C(1, 1) == Approx(1.0));
    REQUIRE(C(2, 2) == Approx(4.0));
    REQUIRE(C(3, 3) == Approx(9.0));
    REQUIRE(C(4, 4) == Approx(16.0));
    REQUIRE(C(1, 2) == Approx(0.0));
    REQUIRE(C(4, 1) == Approx(0.0));
}

TEST_CASE("Full matrix invert diagonal")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(4, 4);
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            A(i, j) = i * j * 1.0 + 1.0;
        }
    }
    FullMatrix<double> C = A.inverseDiagonalMatrix();

    REQUIRE(C(0, 0) == Approx(1.0));
    REQUIRE(C(1, 1) == Approx(0.5));
    REQUIRE(C(2, 2) == Approx(0.2));
    REQUIRE(C(3, 3) == Approx(0.1));
    REQUIRE(C(1, 2) == Approx(0.0));
    REQUIRE(C(3, 1) == Approx(0.0));
}

TEST_CASE("Full matrix eliminateRowAndColumn")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(4, 4);
    A.eliminateRowAndColumn(1);
    REQUIRE(A.cols() == 3);
    REQUIRE(A.rows() == 3);


}

TEST_CASE("Full matrix 2x2 inversion")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(2, 2);
    A(0, 0) = 2.0;
    A(1, 1) = 2.0;
    A(0, 1) = -1.0;
    A.inverse2x2();
    REQUIRE(A(0,0) == 0.5);
    REQUIRE(A(1, 1) == 0.5);
    REQUIRE(A(0, 1) == 0.25);
    REQUIRE(A(1, 0) == 0.0);


}

TEST_CASE("Full matrix 3*3 inversion")
{
    using namespace LinearAlgebra;
    FullMatrix<double> A(3, 3);
    A(0, 0) = A(1,1) = A(2,2) = 2.0;
    A(0, 1) = A(1,2) =  - 1.0;
    A.inverse3x3();
    REQUIRE(A(0, 0) == 0.5);
    REQUIRE(A(1, 1) == 0.5);
    REQUIRE(A(2, 2) == 0.5);
    REQUIRE(A(0, 1) == 0.25);
    REQUIRE(A(1, 2) == 0.25);
    REQUIRE(A(0, 2) == 0.125);
    REQUIRE(A(1, 0) == 0.0);
}

TEST_CASE("Full Matrix reinit")
{
    using namespace LinearAlgebra;
    FullMatrix<float> A;
    REQUIRE(A.cols() == 0);
    A.reinit(5);
    REQUIRE(A.cols() == 5);
}

TEST_CASE("determinant")
{
    using namespace LinearAlgebra;
    FullMatrix<float> A(3,3);
    A(0, 0) = 1.0;
    A(1, 1) = 2.0;
    A(2, 2) = 3.0;

    REQUIRE(A.det() == Approx(6.0));

}