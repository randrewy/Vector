#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

#if defined _MSC_VER
#define CONSTEXPR
#else
#define CONSTEXPR constexpr
#endif



template<typename... T>
  struct all_same : std::false_type { };

template<>
  struct all_same<> : std::true_type { };

template<typename T>
  struct all_same<T> : std::true_type { };

template<typename T, typename... Ts>
  struct all_same<T, T, Ts...> : all_same<T, Ts...> { };

template<typename Type, typename... T>
constexpr bool AllSame(){
	return all_same<std::is_same<T, Type>...>::value;
}


namespace tvec {

enum class Init {raw};
constexpr Init raw = Init::raw;



template<typename Type, size_t N, typename storage_type = std::conditional_t<(N<17), Type[N], std::vector<Type>>>
struct Storage {};

template<typename Type, size_t N>
struct Storage<Type, N, std::enable_if_t<(N<17), Type[N]>> {
protected:
    Type data[N];
public:
    template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N || sizeof...(Args) == 0, void>>
    constexpr Storage(const Args&&... args) : data{static_cast<Type>(args)...} {} // TODO: NOTE HERE!!
    Storage(const Init&) {}
};

template<typename Type, size_t N>
struct Storage<Type, N, std::vector<Type>> {
protected:
    std::vector<Type> data;
public:
    constexpr Storage() : data(N) {}
    template<typename... Args>
    constexpr Storage(const Args&&... args) : data{static_cast<Type>(args)...} {} // TODO: NOTE HERE!!
};



template<class T, size_t N>
struct Vector : public Storage<T, N> {
	using Storage<T, N>::data;

	constexpr size_t size() const noexcept { return N;}

	T& operator[](size_t i) noexcept { return data[i];}
	T operator[](size_t i) const noexcept { return data[i];}

	CONSTEXPR Vector& operator += (const Vector& v) noexcept;
	CONSTEXPR Vector& operator -= (const Vector& v) noexcept;

	template<class Scalar>
	CONSTEXPR Vector& operator *= (Scalar s) noexcept;

	template<class Scalar>
	CONSTEXPR Vector& operator /= (Scalar s) noexcept;

	using Storage<T, N>::Storage;
	constexpr Vector() noexcept {}

	Vector(const Vector& v) noexcept;
};

template<class T, size_t N>
Vector<T, N>::Vector(const Vector<T, N>& v) noexcept {
	for (size_t i = 0; i < N; ++i) {
		data[i] = v[i];
	}
}


template<class T, size_t N>
CONSTEXPR Vector<T, N>& Vector<T, N>::operator += (const Vector<T, N>& v) noexcept {
	for (size_t i = 0; i < N; ++i) {
		data[i] += v.data[i];
	}
	return *this;
}

template<class T, size_t N>
CONSTEXPR Vector<T, N>& Vector<T, N>::operator -= (const Vector<T, N>& v) noexcept {
	for (size_t i = 0; i < N; ++i) {
		data[i] -= v.data[i];
	}
	return *this;
}

template<class T, size_t N>
template<class Scalar>
CONSTEXPR Vector<T, N>& Vector<T, N>::operator *= (Scalar s) noexcept {
	for (size_t i = 0; i < N; ++i) {
		data[i] *= s;
	}
	return *this;
}

template<class T, size_t N>
template<class Scalar>
CONSTEXPR Vector<T, N>& Vector<T, N>::operator /= (Scalar s) noexcept {
	for (size_t i = 0; i < N; ++i) {
		data[i] /= s;
	}
	return *this;
}


template<class T, size_t N>
constexpr inline Vector<T, N>& operator - (const Vector<T, N>& v) noexcept {
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
constexpr inline T operator * (const Vector<T, N>& lhs, const Vector<T, N>& rhs) noexcept {
	T res = 0;
	for (size_t i = 0; i < N; ++i) {
		res += lhs[i] * rhs[i];
	}
	return res;
}

template<class T, size_t N, class S>
constexpr inline Vector<T, N> operator * (Vector<T, N> lhs, S d) noexcept {
	return lhs *= d;
}

template<class T, size_t N, class S>
constexpr inline Vector<T, N> operator * (S d, Vector<T, N> lhs) noexcept{
	return lhs * d;
}


template<class T, size_t N>
constexpr inline T sqlen (const Vector<T, N>& vec) noexcept{
	return vec * vec;
}

template<class T, size_t N>
constexpr inline T len (const Vector<T, N>& vec) noexcept {
	return sqrt(sqlen(vec));
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


template<class T, size_t N, template<typename, size_t> class... Args>
constexpr inline std::enable_if_t<AllSame<Vector<T, N>, Args<T, N>...>(),
Vector<T, N>> cross(const Args<T, N>&...) noexcept {
    static_assert(sizeof...(Args) == N - 1, "N-D space cross-product must have (N-1) arguments");
	return {};
}


//template<typename Type, size_t N_x, size_t N_y>
//struct Matrix { 
//	Type data[N_x * N_y];

//	constexpr Type& opearator () (int i, int j) { return data[i*N_y + j];}
	

//};

}

