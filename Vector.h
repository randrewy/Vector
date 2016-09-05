#ifndef TMX_MATRIX_VECTOR
#define TMX_MATRIX_VECTOR

#include <iostream>
#include <vector>
#include <cmath>

#ifdef __clang__
    #ifndef MSC_VER
        #define CONSTEXPR constexpr
    #endif
#else
#define CONSTEXPR
#endif

namespace tmx {
namespace traits {

template<typename... T>
struct all_same : std::false_type { };

template<typename T>
struct all_same<T, T> : std::true_type { };

template<typename T, typename... Ts>
struct all_same<T, T, Ts...> : all_same<T, Ts...> { };

template<typename Type, typename... T>
constexpr bool AllSame() {
    return all_same<Type, T...>::value;
}

template<typename Type, typename... T>
constexpr bool AllSameDecay() {
    return all_same<std::decay_t<Type>, std::decay_t<T>...>::value;
}

} // namespace traits

constexpr size_t MAX_STACK_SIZE = 17;
enum class Init { raw };
constexpr Init raw = Init::raw;


template<typename Type, size_t N, typename storage_type = std::conditional_t<(N < MAX_STACK_SIZE), Type[N], std::vector<Type>>>
struct Storage;

template<typename Type, size_t N>
struct Storage<Type, N, std::enable_if_t<(N < MAX_STACK_SIZE), Type[N]>> {
protected:
    Type data[N];
public:
   template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N || sizeof...(Args) == 0>>
   constexpr Storage(const Args&... args) noexcept : data {args... } {}
   Storage(const Init&) noexcept {}
};

template<typename Type, size_t N>
struct Storage<Type, N, std::vector<Type>> {
protected:
    std::vector<Type> data;
public:
    constexpr Storage() : data(N) {}
    template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N>>
    constexpr Storage(const Args&... args) : data {args... } {}
};


template<typename Type, size_t N>
using Initializer = const Type(&)[N];

template<typename Type, size_t N_i, size_t N_j>
using Initializer2D = const Type(&)[N_i][N_j];

template<typename Type, size_t N_x, size_t N_y>
struct Storage2D : public Storage<Type, N_x * N_y> {
    using Storage<Type, N_x * N_y>::data;

    CONSTEXPR void init(Initializer2D<Type, N_x, N_y> list) {
        for (size_t i = 0; i < N_x; ++i) {
            for (size_t j = 0; j < N_y; ++j) {
                data[i*N_y + j] = list[i][j];
            }
        }
    }

    using Storage<Type, N_x * N_y>::Storage; // Allow Constexpr init without inner braces

    template<typename T, size_t N_i, size_t N_j> // template guarantee narrowing is not possible
    CONSTEXPR Storage2D(Initializer2D<T, N_i, N_j>& list) {
        static_assert(N_i == N_x, "Row number missmatch");
        static_assert(N_j == N_y, "Col number missmatch");
        init(list);
    }
};


template<class T, size_t N>
class Vector : public Storage<T, N> {
    using Storage<T, N>::data;
public:
    constexpr size_t size() const noexcept { return N; }

    CONSTEXPR T& operator[](size_t i) noexcept { return data[i]; }
    CONSTEXPR T operator[](size_t i) const noexcept { return data[i]; }

    CONSTEXPR inline Vector& operator += (const Vector& v) noexcept;
    CONSTEXPR inline Vector& operator -= (const Vector& v) noexcept;

    template<class Scalar>
    CONSTEXPR inline Vector& operator *= (Scalar s) noexcept;

    template<class Scalar>
    CONSTEXPR inline Vector& operator /= (Scalar s) noexcept;

    constexpr Vector() noexcept(noexcept(Storage<T, N>())) {}

    // disable Vector&/&& here, so copy and move are default (though template instantiated copy c-tor should not
    // prevent move op generating, it is not true for some compilers)
    template<typename... Args, typename = std::enable_if_t<!traits::AllSameDecay<Vector<T, N>, Args...>()> >
    constexpr Vector(const Args&... args) noexcept(noexcept(Storage<T, N>()))
        : Storage<T, N>(args...) {
        static_assert(sizeof...(Args) == N, "Vector<_, N> constructor should be called with exactly N arguments");
    }

    Vector(const Init& init_tag) noexcept : Storage<T, N>(init_tag) {}
};


template<class T, size_t N>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator += (const Vector<T, N>& v) noexcept {
    for (size_t i = 0; i < N; ++i) {
        data[i] += v.data[i];
    }
    return *this;
}

template<class T, size_t N>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator -= (const Vector<T, N>& v) noexcept {
    for (size_t i = 0; i < N; ++i) {
        data[i] -= v.data[i];
    }
    return *this;
}

template<class T, size_t N>
template<class Scalar>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator *= (Scalar s) noexcept {
    for (size_t i = 0; i < N; ++i) {
        data[i] *= s;
    }
    return *this;
}

template<class T, size_t N>
template<class Scalar>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator /= (Scalar s) noexcept {
    for (size_t i = 0; i < N; ++i) {
        data[i] /= s;
    }
    return *this;
}


template<class T, size_t N>
CONSTEXPR inline Vector<T, N>& operator - (const Vector<T, N>& v) noexcept {
    Vector<T, N> tmp;
    for (size_t i = 0; i < N; ++i) {
        tmp[i] = -v[i];
    }
    return tmp;
}

template<class T, size_t N>
constexpr inline Vector<T, N> operator - (Vector<T, N> lhs, const Vector<T, N>& rhs) noexcept {
    return lhs -= rhs;
}

template<class T, size_t N>
constexpr inline Vector<T, N> operator + (Vector<T, N> lhs, const Vector<T, N>& rhs) noexcept {
    return lhs += rhs;
}


template<class T, size_t N>
CONSTEXPR inline T operator * (const Vector<T, N>& lhs, const Vector<T, N>& rhs) noexcept {
    T res = 0;
    for (size_t i = 0; i < N; ++i) {
        res += lhs[i] * rhs[i];
    }
    return res;
}

template<class T, size_t N, class S>
CONSTEXPR inline Vector<T, N> operator * (Vector<T, N> lhs, S d) noexcept {
    return lhs *= d;
}

template<class T, size_t N, class S>
CONSTEXPR inline Vector<T, N> operator * (S d, Vector<T, N> lhs) noexcept {
    return lhs * d;
}

template<class T, size_t N, class S>
CONSTEXPR inline Vector<T, N> operator / (Vector<T, N> lhs, S d) noexcept {
    return lhs *= 1./d;
}

template<class T, size_t N>
constexpr inline T sqlen(const Vector<T, N>& vec) noexcept {
    return vec * vec;
}

template<class T, size_t N>
constexpr inline T len(const Vector<T, N>& vec) noexcept {
    return sqrt(sqlen(vec));
}

template<class T>
constexpr Vector<T, 3> cross (const Vector<T, 3>& a, const Vector<T, 3>& b) {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}

template<class T, size_t N>
constexpr Vector<T, 3> cross (Initializer<T, N>& a, Initializer<T, N>& b) {
    static_assert(N == 3, "3D cross product should have 3D vector arguments");
    return cross(Vector<T, 3>{a[0], a[1], a[2]}, Vector<T, 3>{b[0], b[1], b[2]});
}

template<class T, size_t N, template<typename, size_t> class... Args>
constexpr inline std::enable_if_t<traits::AllSame<Vector<T, N>, Args<T, N>...>(),
Vector<T, N>> general_cross(const Args<T, N>&... args) noexcept {
    static_assert(sizeof...(Args) == N - 1, "N-D space cross-product must have (N-1) arguments");
    static_assert(sizeof...(Args) == 2, "general_cross is not implemented for N != 3");
    return cross(args...);
}

template<typename Type, size_t N_x, size_t N_y>
struct Matrix : public Storage2D<Type, N_x, N_y> {
    using Storage2D<Type, N_x, N_y>::data;
public:

    constexpr size_t sizeX() const { return N_x; }
    constexpr size_t sizeY() const { return N_y; }
    constexpr size_t size() const { return N_x * N_y; }

    CONSTEXPR Type& operator () (size_t i, size_t j) { return data[i*N_y + j]; }
    CONSTEXPR Type operator () (size_t i, size_t j) const { return data[i*N_y + j]; }


    CONSTEXPR inline Matrix& operator += (const Matrix& v) noexcept;
    CONSTEXPR inline Matrix& operator -= (const Matrix& v) noexcept;

    template<class Scalar>
    CONSTEXPR inline Matrix& operator *= (Scalar s) noexcept;

    template<class Scalar>
    CONSTEXPR inline Matrix& operator /= (Scalar s) noexcept;


    using Storage2D<Type, N_x, N_y>::Storage2D;
    constexpr Matrix() noexcept : Storage2D<Type, N_x, N_y>() {}
};

template<class T, size_t N_x, size_t N_y>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator += (const Matrix<T, N_x, N_y>& m) noexcept {
    for (size_t i = 0; i < size(); ++i) {
        data[i] += m.data[i];
    }
    return *this;
}

template<class T, size_t N_x, size_t N_y>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator -= (const Matrix<T, N_x, N_y>& m) noexcept {
    for (size_t i = 0; i < size(); ++i) {
        data[i] -= m.data[i];
    }
    return *this;
}

template<class T, size_t N_x, size_t N_y>
template<class Scalar>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator *= (Scalar s) noexcept {
    for (size_t i = 0; i < size(); ++i) {
        data[i] *= s;
    }
    return *this;
}

template<class T, size_t N_x, size_t N_y>
template<class Scalar>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator /= (Scalar s) noexcept {
    for (size_t i = 0; i < size(); ++i) {
        data[i] /= s;
    }
    return *this;
}


template<class T, size_t N_i, size_t N_j>
constexpr inline Matrix<T, N_i, N_j> operator + (Matrix<T, N_i, N_j> m1, const Matrix<T, N_i, N_j>& m2) noexcept {
    return m1 += m2;
}

template<class T, size_t N_i, size_t N_j>
constexpr inline Matrix<T, N_i, N_j> operator - (Matrix<T, N_i, N_j> m1, const Matrix<T, N_i, N_j>& m2) noexcept {
    return m1 -= m2;
}

template<class T, size_t N_i, size_t N_j>
CONSTEXPR inline Matrix<T, N_i, N_j> operator - (const Matrix<T, N_i, N_j>& m) noexcept {
    Matrix<T, N_i, N_j> tmp;
    for (size_t i = 0; i < N_i; ++i) {
        for (size_t j = 0; j < N_j; ++j) {
            tmp(i, j) = -m(i, j);
        }
    }
    return tmp;
}

template<class T, size_t N_i, size_t N_j, class Scalar>
constexpr inline Matrix<T, N_i, N_j> operator * (Matrix<T, N_i, N_j> m, const Scalar& s) noexcept {
    return m *= s;
}

template<class T, size_t N_i, size_t N_j, class Scalar>
constexpr inline Matrix<T, N_i, N_j> operator * (const Scalar& s, const Matrix<T, N_i, N_j>& m) noexcept {
    return m * s;
}

template<class T, size_t N_i, size_t N_j, class Scalar>
constexpr inline Matrix<T, N_i, N_j> operator / (Matrix<T, N_i, N_j> m, const Scalar& s) noexcept {
    return m /= s;
}

// naive n^3
template<class T, size_t N_i, size_t N_j, size_t N_k>
CONSTEXPR inline Matrix<T, N_i, N_k> operator * (const Matrix<T, N_i, N_j>& m1, const Matrix<T, N_j, N_k>& m2) noexcept {
    Matrix<T, N_i, N_k> result;
    for (size_t i = 0; i < N_i; ++i) {
        for (size_t k = 0; k < N_k; ++k) {
            T t{};
            for (size_t j = 0; j < N_j; ++j) {
                t += m1(i, j) * m2(j, k);
            }
            result(i, k) = t;
        }
    }
    return result;
}

template<class T, size_t N_i, size_t N_j>
CONSTEXPR inline Vector<T, N_i> operator * (const Matrix<T, N_i, N_j>& m, const Vector<T, N_j>& v) noexcept {
    Vector<T, N_i> result;
    for (size_t i = 0; i < N_i; ++i) {
        for (size_t j = 0; j < N_j; ++j) {
            result[i] += m(i, j) * v[j];
        }
    }
    return result;
}

template<class T, size_t N_i, size_t N_j>
CONSTEXPR inline Vector<T, N_j> operator * (const Vector<T, N_i>& v, const Matrix<T, N_i, N_j>& m) noexcept {
    Vector<T, N_j> result;
    for (size_t i = 0; i < N_i; ++i) {
        for (size_t j = 0; j < N_j; ++j) {
            result[j] += m(i, j) * v[i];
        }
    }
    return result;
}


template<class T, size_t N>
std::ostream& operator << (std::ostream& os, const Vector<T, N>& v) {
    os << "(";
    for (size_t i = 0; i < N - 1; ++i) {
       os << v[i] << ", ";
    }
    if (0 < N) {
       os << v[N - 1];
    }
    return os << ")";
}

template<class T, size_t N_i, size_t N_j>
std::ostream& operator << (std::ostream& os, const Matrix<T, N_i, N_j>& m) {
    os << "{";
    for (size_t i = 0; i < N_i; ++i) {
        os << "(";
        for (size_t j = 0; j < N_j - 1; ++j) {
            os << m(i, j) << ", ";
        }
        if (0 < N_j) {
            os << m(i, N_j - 1) << ")";
        }
        if (i != N_i - 1) {
            os << ", ";
        }
    }
    return os << "}";
}

} // namespace tmx


#endif // !TMX_MATRIX_VECTOR
