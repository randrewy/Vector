#ifndef TMX_MATRIX_VECTOR
#define TMX_MATRIX_VECTOR

#include <iostream>
#include <vector>
#include <cmath>

#if    (defined(__clang_major__) && (__clang_major__ * 10 + __clang_minor__ >= 32)) \
    || (defined (_MSC_VER) && _MSC_VER >= 1700 && 0)/* MSC still has some troubles*/ \
    || (defined (__GNUC__) && __GNUC__ >= 5)
    #define CONSTEXPR constexpr
#else
#define CONSTEXPR
#endif // __clang__

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

enum class RawInitTag { raw };
constexpr RawInitTag raw = RawInitTag::raw;

// QtCreator goes crazy if you use '<' in "N < MAX_STACK_SIZE
template<typename Type, size_t N, typename = std::conditional_t<(N <= MAX_STACK_SIZE - 1), Type[N], std::vector<Type>>>
struct Storage;


// is it better than anonymous struct? (which is deprecated since c++11)
// though any good compiler will eliminate all branches from switch
#ifdef TMX_USE_NAMED_MEMBERS
template<typename Type>
struct alignas(2 * sizeof(Type)) Storage<Type, 2, Type[2]> {
    Type x;
    union { Type y; Type __constexpt_fail; };

    constexpr Storage(const Type& x = {}, const Type& y = {}) : x(x), y(y) {}
    Storage(const RawInitTag&) {}

    constexpr Type& operator[](size_t i) {
        switch (i) {
        case 0: return x;
        case 1: return y;
        default: return __constexpt_fail;
        }
    }

    constexpr Type operator[](size_t i) const {
        switch (i) {
        case 0: return x;
        case 1: return y;
        default: return __constexpt_fail;
        }
    }

    const Type* plainData() const { return &x; }
};

template<typename Type>
struct alignas(4 * sizeof(Type)) Storage<Type, 3, Type[3]> {
    Type x;
    Type y;
    union { Type z; Type __constexpt_fail; };

    constexpr Storage(const Type& x = {}, const Type& y = {}, const Type& z = {}) : x(x), y(y), z(z) {}
    Storage(const RawInitTag&) {}

    constexpr Type& operator[](size_t i) {
        switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: return __constexpt_fail;
        }
    }

    constexpr Type operator[](size_t i) const {
        switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: return __constexpt_fail;
        }
    }

    const Type* plainData() const { return &x; }
};
#endif // TMX_USE_NAMED_MEMBERS

template<typename Type, size_t N>
struct alignas(sizeof(Type) >= 16 ? 4 * sizeof(Type) : 8 * sizeof(Type)) Storage<Type, N, Type[N]> {
protected:
    Type data[N];
public:
    template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N || sizeof...(Args) == 0>>
    constexpr Storage(const Args&... args) : data {args... } {}
   
    Storage(const RawInitTag&) {}

    constexpr Type& operator[](size_t i) { return data[i]; }
    constexpr Type  operator[](size_t i) const { return data[i]; }

    const Type* plainData() const { return data; }
};

template<typename Type, size_t N>
struct Storage<Type, N, std::vector<Type>> {
protected:
    std::vector<Type> data;
public:
    constexpr Storage() : data(N) {}
    template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N>>
    constexpr Storage(const Args&... args) : data {args... } {}

    constexpr Type& operator[](size_t i) { return data[i]; }
    constexpr Type  operator[](size_t i) const { return data[i]; }

    const Type* plainData() const { return data.data(); }
};


template<typename Type, size_t N>
using Initializer = const Type(&)[N];

template<typename Type, size_t N_i, size_t N_j>
using Initializer2D = const Type(&)[N_i][N_j];

template<typename Type, size_t N_x, size_t N_y>
struct Storage2D : public Storage<Type, N_x * N_y> {
protected:
    constexpr Type& get(size_t i, size_t j) { return Storage<Type, N_x * N_y>::operator[](i*N_y + j); }
    constexpr Type  get(size_t i, size_t j) const { return Storage<Type, N_x * N_y>::operator[](i*N_y + j); }
    constexpr Type& get(size_t i) { return Storage<Type, N_x * N_y>::operator[](i); }
    constexpr Type  get(size_t i) const { return Storage<Type, N_x * N_y>::operator[](i); }

    constexpr void init(Initializer2D<Type, N_x, N_y> list) {
        for (size_t i = 0; i < N_x; ++i) {
            for (size_t j = 0; j < N_y; ++j) {
                get(i, j) = list[i][j];
            }
        }
    }
public:
    constexpr Storage2D() : Storage<Type, N_x * N_y>() {}

    template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N_x * N_y>>
    constexpr Storage2D(const Args&... args) : Storage<Type, N_x * N_y>( args... ) {} // Allow Constexpr init without inner braces

    template<typename T, size_t N_i, size_t N_j> // template guarantee narrowing is not possible
    constexpr Storage2D(Initializer2D<T, N_i, N_j>& list) {
        static_assert(N_i == N_x, "Row number missmatch");
        static_assert(N_j == N_y, "Col number missmatch");
        init(list);
    }
};


template<class T, size_t N>
class Vector : public Storage<T, N> {
public:
    using Storage<T, N>::operator[];

    constexpr size_t size() const { return N; }

    CONSTEXPR Vector& operator += (const Vector& v);
    CONSTEXPR Vector& operator -= (const Vector& v);

    template<class Scalar>
    CONSTEXPR Vector& operator *= (Scalar s);

    template<class Scalar>
    CONSTEXPR Vector& operator /= (Scalar s);

    constexpr Vector() {}

    // disable Vector&/&& here, so copy and move are default (though template instantiated copy c-tor should not
    // prevent move op generating, it is not true for some compilers)
    template<typename... Args, typename = std::enable_if_t<!traits::AllSameDecay<Vector<T, N>, Args...>()> >
    constexpr Vector(const Args&... args)
        : Storage<T, N>(args...) {
        static_assert(sizeof...(Args) == N, "Vector<_, N> constructor should be called with exactly N arguments");
    }

    constexpr Vector(const Vector& vec) = default;
    constexpr Vector(Vector&& vec) = default;

    Vector(const RawInitTag& init_tag) : Storage<T, N>(init_tag) {}

    const T* data() const { return Storage<T, N>::plainData(); }
};


template<class T, size_t N>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator += (const Vector<T, N>& v) {
    for (size_t i = 0; i < N; ++i) {
        operator[](i) += v[i];
    }
    return *this;
}

template<class T, size_t N>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator -= (const Vector<T, N>& v) {
    for (size_t i = 0; i < N; ++i) {
        operator[](i) -= v[i];
    }
    return *this;
}

template<class T, size_t N>
template<class Scalar>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator *= (Scalar s) {
    for (size_t i = 0; i < N; ++i) {
        operator[](i) *= s;
    }
    return *this;
}

template<class T, size_t N>
template<class Scalar>
CONSTEXPR inline Vector<T, N>& Vector<T, N>::operator /= (Scalar s) {
    for (size_t i = 0; i < N; ++i) {
        operator[](i) /= s;
    }
    return *this;
}


template<class T, size_t N>
CONSTEXPR inline Vector<T, N>& operator - (const Vector<T, N>& v) {
    Vector<T, N> tmp;
    for (size_t i = 0; i < N; ++i) {
        tmp[i] = -v[i];
    }
    return tmp;
}

template<class T, size_t N>
constexpr inline Vector<T, N> operator - (Vector<T, N> lhs, const Vector<T, N>& rhs) {
    return lhs -= rhs;
}

template<class T, size_t N>
constexpr inline Vector<T, N> operator + (Vector<T, N> lhs, const Vector<T, N>& rhs) {
    return lhs += rhs;
}


template<class T, size_t N>
CONSTEXPR inline T operator * (const Vector<T, N>& lhs, const Vector<T, N>& rhs) {
    T res = 0;
    for (size_t i = 0; i < N; ++i) {
        res += lhs[i] * rhs[i];
    }
    return res;
}

template<class T, size_t N, class S>
CONSTEXPR inline Vector<T, N> operator * (Vector<T, N> lhs, S d) {
    return lhs *= d;
}

template<class T, size_t N, class S>
CONSTEXPR inline Vector<T, N> operator * (S d, Vector<T, N> lhs) {
    return lhs * d;
}

template<class T, size_t N, class S>
CONSTEXPR inline Vector<T, N> operator / (Vector<T, N> lhs, S d) {
    return lhs *= 1./d;
}

template<class T, size_t N>
constexpr inline T sqlen(const Vector<T, N>& vec) {
    return vec * vec;
}

template<class T, size_t N>
constexpr inline T len(const Vector<T, N>& vec) {
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
Vector<T, N>> general_cross(const Args<T, N>&... args) {
    static_assert(sizeof...(Args) == N - 1, "N-D space cross-product must have (N-1) arguments");
    static_assert(sizeof...(Args) == 2, "general_cross is not implemented for N != 3");
    return cross(args...);
}

template<typename Type, size_t N_x, size_t N_y>
struct Matrix : public Storage2D<Type, N_x, N_y> {
    using Storage2D<Type, N_x, N_y>::get;
public:

    constexpr size_t sizeX() const { return N_x; }
    constexpr size_t sizeY() const { return N_y; }
    constexpr size_t size() const { return N_x * N_y; }

    CONSTEXPR Type& operator () (size_t i, size_t j) { return get(i, j); }
    CONSTEXPR Type operator () (size_t i, size_t j) const { return get(i, j); }

    CONSTEXPR Matrix& operator += (const Matrix& v);
    CONSTEXPR Matrix& operator -= (const Matrix& v);

    template<class Scalar>
    CONSTEXPR Matrix& operator *= (Scalar s);

    template<class Scalar>
    CONSTEXPR Matrix& operator /= (Scalar s);


    using Storage2D<Type, N_x, N_y>::Storage2D;
    constexpr Matrix() : Storage2D<Type, N_x, N_y>() {}

    const Type* data() const { return Storage2D<Type, N_x, N_y>::plainData(); }
};

template<class T, size_t N_x, size_t N_y>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator += (const Matrix<T, N_x, N_y>& m) {
    for (size_t i = 0; i < size(); ++i) {
        get(i) += m.get(i);
    }
    return *this;
}

template<class T, size_t N_x, size_t N_y>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator -= (const Matrix<T, N_x, N_y>& m) {
    for (size_t i = 0; i < size(); ++i) {
        get(i) -= m.get(i);
    }
    return *this;
}

template<class T, size_t N_x, size_t N_y>
template<class Scalar>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator *= (Scalar s) {
    for (size_t i = 0; i < size(); ++i) {
        get(i) *= s;
    }
    return *this;
}

template<class T, size_t N_x, size_t N_y>
template<class Scalar>
CONSTEXPR inline Matrix<T, N_x, N_y>& Matrix<T, N_x, N_y>::operator /= (Scalar s) {
    for (size_t i = 0; i < size(); ++i) {
        get(i) /= s;
    }
    return *this;
}


template<class T, size_t N_i, size_t N_j>
constexpr inline Matrix<T, N_i, N_j> operator + (Matrix<T, N_i, N_j> m1, const Matrix<T, N_i, N_j>& m2) {
    return m1 += m2;
}

template<class T, size_t N_i, size_t N_j>
constexpr inline Matrix<T, N_i, N_j> operator - (Matrix<T, N_i, N_j> m1, const Matrix<T, N_i, N_j>& m2) {
    return m1 -= m2;
}

template<class T, size_t N_i, size_t N_j>
CONSTEXPR inline Matrix<T, N_i, N_j> operator - (const Matrix<T, N_i, N_j>& m) {
    Matrix<T, N_i, N_j> tmp;
    for (size_t i = 0; i < N_i; ++i) {
        for (size_t j = 0; j < N_j; ++j) {
            tmp(i, j) = -m(i, j);
        }
    }
    return tmp;
}

template<class T, size_t N_i, size_t N_j, class Scalar>
constexpr inline Matrix<T, N_i, N_j> operator * (Matrix<T, N_i, N_j> m, const Scalar& s) {
    return m *= s;
}

template<class T, size_t N_i, size_t N_j, class Scalar>
constexpr inline Matrix<T, N_i, N_j> operator * (const Scalar& s, const Matrix<T, N_i, N_j>& m) {
    return m * s;
}

template<class T, size_t N_i, size_t N_j, class Scalar>
constexpr inline Matrix<T, N_i, N_j> operator / (Matrix<T, N_i, N_j> m, const Scalar& s) {
    return m /= s;
}

// naive n^3
template<class T, size_t N_i, size_t N_j, size_t N_k>
CONSTEXPR inline Matrix<T, N_i, N_k> operator * (const Matrix<T, N_i, N_j>& m1, const Matrix<T, N_j, N_k>& m2) {
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
CONSTEXPR inline Vector<T, N_i> operator * (const Matrix<T, N_i, N_j>& m, const Vector<T, N_j>& v) {
    Vector<T, N_i> result;
    for (size_t i = 0; i < N_i; ++i) {
        for (size_t j = 0; j < N_j; ++j) {
            result[i] += m(i, j) * v[j];
        }
    }
    return result;
}

template<class T, size_t N_i, size_t N_j>
CONSTEXPR inline Vector<T, N_j> operator * (const Vector<T, N_i>& v, const Matrix<T, N_i, N_j>& m) {
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

#undef CONSTEXPR
#endif // !TMX_MATRIX_VECTOR
