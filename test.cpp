#include "Vector.h"
#include "gtest/gtest.h"

using Vector = tmx::Vector<double, 3>;
using Matrix = tmx::Matrix<double, 3, 3>;

TEST(VectorTest, plus) {
    Vector a{ 1., 2., 3. };
    Vector b{ -1., 1., -1. };
    Vector r = a + b;

    ASSERT_DOUBLE_EQ(0., r[0]);
    ASSERT_DOUBLE_EQ(3., r[1]);
    ASSERT_DOUBLE_EQ(2., r[2]);
}

TEST(VectorTest, minus) {
    Vector a{ 1., 2., 3. };
    Vector b{ -1., 1., -1. };
    Vector r = a - b;

    ASSERT_DOUBLE_EQ(2., r[0]);
    ASSERT_DOUBLE_EQ(1., r[1]);
    ASSERT_DOUBLE_EQ(4., r[2]);
}

TEST(VectorTest, unary_minus) {
    Vector a{ 1., 2., 3. };

    auto r = -a;

    ASSERT_DOUBLE_EQ(-1., r[0]);
    ASSERT_DOUBLE_EQ(-2., r[1]);
    ASSERT_DOUBLE_EQ(-3., r[2]);
}

TEST(VectorTest, mul) {
    Vector a{ 1., 2., -3. };
    Vector r = a * 2;

    ASSERT_DOUBLE_EQ(2., r[0]);
    ASSERT_DOUBLE_EQ(4., r[1]);
    ASSERT_DOUBLE_EQ(-6., r[2]);
}

TEST(VectorTest, div) {
    Vector a{ 1., 2., -3. };
    Vector r = a / 2.;

    ASSERT_DOUBLE_EQ(0.5, r[0]);
    ASSERT_DOUBLE_EQ(1., r[1]);
    ASSERT_DOUBLE_EQ(-1.5, r[2]);
}

TEST(VectorTest, dot) {
    Vector a{ 1., 2., -3. };
    Vector b{ -1., 1., -1. };

    double r = a * b;
    ASSERT_DOUBLE_EQ(4., r);
}

TEST(VectorTest, cross) {
    Vector a{ 1., 2., -3. };
    Vector b{ -1., 1., -1. };

    Vector r = cross(a, b);

    ASSERT_DOUBLE_EQ(1., r[0]);
    ASSERT_DOUBLE_EQ(4., r[1]);
    ASSERT_DOUBLE_EQ(3, r[2]);
}

TEST(VectorTest, len) {
    Vector a{ 1., 2., -3. };

    ASSERT_DOUBLE_EQ(3.74165738677394, len(a));
    ASSERT_DOUBLE_EQ(14., sqlen(a));
}


class MatrixTest : public ::testing::Test {
public:
    Matrix i {{
        { 1., 0., 0. },
        { 0., 1., 0. },
        { 0., 0., 1. },
    }};

    Matrix m1 {{
        {  2., -1., 3. },
        { -2.,  0., 1. },
        { -1., -3., 2. },
    }};

    Matrix m2 {{
        { -1., -2.,  2. },
        { -3., -1.,  3. },
        {  2., -1., -3. },
    }};
};

#define ASSERT_MATRIX_EQ(m1, m2)\
    ASSERT_EQ(m1.sizeX(), (m2).sizeX());\
    ASSERT_EQ(m1.sizeY(), (m2).sizeY());\
    ASSERT_DOUBLE_EQ(m1(0, 0), (m2)(0, 0));\
    ASSERT_DOUBLE_EQ(m1(0, 1), (m2)(0, 1));\
if (m1.sizeY() > 2 && (m2).sizeY() > 2)\
    ASSERT_DOUBLE_EQ(m1(0, 2), (m2)(0, 2));\
    ASSERT_DOUBLE_EQ(m1(1, 0), (m2)(1, 0));\
    ASSERT_DOUBLE_EQ(m1(1, 1), (m2)(1, 1));\
if (m1.sizeY() > 2 && (m2).sizeY() > 2) \
    ASSERT_DOUBLE_EQ(m1(1, 2), (m2)(1, 2));\
if (m1.sizeX() > 2 && (m2).sizeX() > 2) \
    ASSERT_DOUBLE_EQ(m1(2, 0), (m2)(2, 0));\
if (m1.sizeX() > 2 && (m2).sizeX() > 2) \
    ASSERT_DOUBLE_EQ(m1(2, 1), (m2)(2, 1));\
if (m1.sizeX() > 2 && (m2).sizeX() > 2 && m1.sizeY() > 2 && (m2).sizeY() > 2) \
    ASSERT_DOUBLE_EQ(m1(2, 2), (m2)(2, 2))


TEST_F(MatrixTest, plus) {
    Matrix check {{
        {  1., -3.,  5. },
        { -5., -1.,  4. },
        {  1., -4., -1. }
    }};

    ASSERT_MATRIX_EQ(check, m1 + m2);
}

TEST_F(MatrixTest, minus) {
    Matrix check {{
        {  3.,  1.,  1. },
        {  1.,  1., -2. },
        { -3., -2.,  5. }
    }};

    ASSERT_MATRIX_EQ(check, m1 - m2);
}

TEST_F(MatrixTest, unary_minus) {
    Matrix check {{
        {  1., 2.,  -2. },
        {  3., 1.,  -3. },
        { -2., 1.,   3. },
    }};

    ASSERT_MATRIX_EQ(check, -m2);

}

TEST_F(MatrixTest, div) {
    Matrix check {{
        {   0.5, -0.25, 0.75 },
        {  -0.5,    0., 0.25 },
        { -0.25, -0.75, 0.5 }
    }};

    ASSERT_MATRIX_EQ(check, m1 / 4);
}

TEST_F(MatrixTest, mul_scalar) {
    Matrix check {{
        { 4., -2., 6. },
        { -4.,  0., 2. },
        { -2., -6., 4. },
    }};

    ASSERT_MATRIX_EQ(check, m1 * 2);
}


TEST_F(MatrixTest, mul_matrix) {
    Matrix check {{
        {  7., -6.,  -8. },
        {  4.,  3.,  -7. },
        { 14.,  3., -17. }
    }};

    ASSERT_MATRIX_EQ(check, m1 * m2);
}


TEST_F(MatrixTest, mul_matrix_non_square) {
    tmx::Matrix<double, 3, 2> mx1 {{
        {  2., -1. },
        { -2.,  0. },
        { -1., -3. }
    }};

    tmx::Matrix<double, 2, 3> mx2{ {
        { 3., -1., -3.},
        { 1.,  2.,  3. },
    }};

    Matrix check1 {{
        {  5., -4., -9. },
        { -6.,  2.,  6. },
        { -6., -5., -6. }
    }};

    tmx::Matrix<double, 2, 2> check2 {{
        { 11.,   6. },
        { -5., -10. }
    }};

    ASSERT_MATRIX_EQ(check1, mx1 * mx2);
    ASSERT_MATRIX_EQ(check2, mx2 * mx1);
}

TEST(MatrixVectorTest, mul_right) {
    Matrix m1 {{
        { 2., -1., 3. },
        { -2.,  0., 1. },
        { -1., -3., 2. },
    }};
    Vector a { 1., 2., -3. };

    Vector check{ -9., -5., -13. };

    ASSERT_DOUBLE_EQ(check[0], (m1 * a)[0]);
    ASSERT_DOUBLE_EQ(check[1], (m1 * a)[1]);
    ASSERT_DOUBLE_EQ(check[2], (m1 * a)[2]);
}

TEST(MatrixVectorTest, mul_left) {
    Matrix m1 {{
        { 2., -1., 3. },
        { -2.,  0., 1. },
        { -1., -3., 2. },
    }};
    Vector a { 1., 2., -3. };

    Vector check{ 1., 8., -1. };

    ASSERT_DOUBLE_EQ(check[0], (a * m1)[0]);
    ASSERT_DOUBLE_EQ(check[1], (a * m1)[1]);
    ASSERT_DOUBLE_EQ(check[2], (a * m1)[2]);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
 }
