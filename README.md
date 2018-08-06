## Simple Vector ##

This is a simple C++14 header-only library for matrix computation. It includes only the smallest set
of vital operations.

Both Matrix and Vector have fixed size, defined by the template parameters. Vectors and matrices up 
to 16 elements are stored without heap allocation and almost all operations for them are constexpr.

### Installation ###
No build is required. Just pull repository wherever you want and include `Vector.h` into your project.

### Tests ###
You'll need [googletest](https://github.com/google/googletest/blob/master/googletest/).
Compile and run tests.cpp. There are some basic Vector-Vector, Vector-Matrix and Matrix-Matrix tests.

Just for now there is no build system supported and it has to be done manualy, 
so don't forget to link gtest library and add proper include paths.

### Named Members ###
Vectors of size 2 and 3 has named members `.x`, `.y` and `.z` for Vector3. This implizes a little overhead if
optimization level is not high, but all branching is eliminated completely on high level. If you
don't need this feature or you are afraid of this tiny perfomance overhead in debug define `TMX_DONT_USE_NAMED_MEMBERS`
before you include library header.

### A taste of Vector.h ###

```c++
using Vector3 = tmx::Vector<double, 3>;
using Matrix3d = tmx::Matrix<double, 3, 3>;
constexpr Vector3 v1;
constexpr Vector3 v2 {1., 2., -1.};

constexpr Matrix3d m1 {
    1.,2.,3.,
    1.,2.,3.,
    1.,2.,3.,
};
```
A bit more precise matrix creation. Note double `{{`:
```c++
constexpr Matrix3d m2 {{
    {1., 0., 0.},
    {0., 1., 0.},
    {0., 0., 1.}
}};
```
Let's rotate a vector and print result:
```c++
constexpr Matrix3d rot {{
    { 0., -1., 0.},
    { 1.,  0., 0.},
    { 0.,  0., 1.}
}};

constexpr Vector3 direction = {1., 0., 0.};

constexpr auto result = rot * direction;
std::cout << result;
```
Should print `(0, 1, 0)`