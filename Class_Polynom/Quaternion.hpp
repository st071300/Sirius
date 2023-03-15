#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <complex>

using namespace std;

template<typename T>

class Quaternion
{
public:

    T scalar = 0, x = 0, y = 0, z = 0;
    Quaternion() = default;

    Quaternion<T> operator+(const Quaternion<T> &q1)
    {
        return Quaternion<T>{scalar + q1.scalar, x + q1.x, y + q1.y, z + q1.z};
    }


    Quaternion<T> operator-(const Quaternion<T> &q1)
    {
        return Quaternion<T>{scalar - q1.scalar, x - q1.x, y - q1.y, z - q1.z};
    }


    Quaternion<T> operator*(const Quaternion<T> &q1)
    {
        return Quaternion<T>{scalar * q1.scalar - x * q1.x - y * q1.y - z * q1.z,
                       scalar * q1.x + x * q1.scalar + y * q1.z - z * q1.y,
                       scalar * q1.y - x * q1.z + y * q1.scalar + z * q1.x,
                       scalar * q1.z + x * q1.y - y * q1.x + z * q1.scalar};
    }

    Quaternion<T> operator*(const T &sc)
    {
        return Quaternion<T>{scalar * sc + x * sc + y * sc + z * sc};
    }


};


template <class T> std::ostream& operator<<(std::ostream& os, const Quaternion<T>& quat)
{
    os << "{" << quat.scalar << " [" << quat.x<< " " << quat.y<< " " << quat.z<<"]}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::complex<double> & cplx)
{
    os << std::real(cplx) << " + " << std::imag(cplx) <<"i";
    return os;
}