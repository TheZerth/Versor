//
// Created by zerth on 12/12/24.
//

#include "Versor.h"

#include <cmath>

namespace VRSR {

    enum VectorComponents { X = 0, Y = 1, Z = 2};

    //These are the basis vectors used for the space. All multivectors are combinations of e1 and e2.
    static const std::vector<float> e1 = {1.0f,0.0f,0.0f};
    static const std::vector<float> e2 = {0.0f,1.0f,0.0f};
    static const std::vector<float> e3 = {0.0f,0.0f,1.0f};

    //--------------------Reverse--------------------
    Versor Versor::reverse() const {
        return Versor{a, x, y, z, -b, -c, -d, -e};
    }

    //--------------------Dual--------------------
    Versor Versor::dual() const {
        return Versor{-e, -d, c, -b, z, -y, x, a};
    }

    //--------------------Conjugate--------------------
    Versor Versor::conjugate() const {
        return {a, -x, -y, -z, -b, -c, -d, e};
    }

    //--------------------Involute--------------------
    Versor Versor::involute() const {
        return {a, -x, -y, -z, b, c, d, -e};
    }

    //--------------------SqNorm--------------------
    double Versor::sqNorm() const {
        const Versor res = this->mul(this->conjugate());
        return std::sqrt(std::abs(res.a));
    }

    //--------------------Normalized--------------------
    Versor Versor::normalize() const {
        const double n = sqNorm();
        if (n == 0) return *this;
        return {a / n, x / n, y / n, z / n, b / n, c / n, d / n, e / n};
    }

    //--------------------Magnitude--------------------
    double Versor::magnitude() const {
        return std::sqrt((a * a) + (x * x) + (y * y) + (z * z) + (b * b) + (c * c) + (d * d) + (e * e));
    }

    //--------------------Inverse--------------------
    Versor Versor::inverse() const {
        const double n = this->magnitude();
        if (n == 0) return *this;
        return {a / n, x / n, y / n, z / n, -b / n, -c / n, -d / n, -e / n};
    }

    //--------------------Addition--------------------
    template <typename T>
    Versor Versor::add(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a + t.a, x + t.x, y + t.y, z + t.z, b + t.b, c + t.c, d + t.d, e + t.e};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a + t, x, y, z, b, c, d, e};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 3:
                    return {a, x + t[X], y + t[Y], z + t[Z], b, c, d, e};
                case 2:
                    return {a, x + t[X], y + t[Y], z, b, c, d, e};
                case 1:
                    return {a, x + t[X], y, z, b, c, d, e};
                default:
                    throw std::invalid_argument("Invalid size for addition, maximum size is 3.");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for addition.");
        }
    }
    template Versor Versor::add<Versor>(const Versor &t) const;
    template Versor Versor::add<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::add<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::add<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::add<double>(const double &t) const;
    template Versor Versor::add<float>(const float &t) const;
    template Versor Versor::add<int>(const int &t) const;

    //--------------------Subtraction--------------------
    template <typename T>
    Versor Versor::sub(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a - t.a, x - t.x, y - t.y, z - t.z, b - t.b, c - t.c, d - t.d, e - t.e};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a - t, x, y, z, b, c, d, e};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 3:
                    return {a, x - t[X], y - t[Y], z - t[Z], b, c, d, e};
                case 2:
                    return {a, x - t[X], y - t[Y], z, b, c, d, e};
                case 1:
                    return {a, x - t[X], y, z, b, c, d, e};
                default:
                    throw std::invalid_argument("Invalid size for subtraction, maximum size is 3.");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for subtraction.");
        }
    }
    template Versor Versor::sub<Versor>(const Versor &t) const;
    template Versor Versor::sub<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::sub<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::sub<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::sub<double>(const double &t) const;
    template Versor Versor::sub<float>(const float &t) const;
    template Versor Versor::sub<int>(const int &t) const;

    //--------------------Left Contraction--------------------
    template <typename T>
    Versor Versor::lco(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {(a * t.a) + (x * t.x) + (y * t.y) - (b * t.b) - (c * t.c) - (d * t.d) - (e * t.e),
                (a * t.x) - (y * t.b) - (z * t.c) - (d * t.e),
                (a * t.y) + (x * t.b) - (z * t.d) - (c * t.e),
                (a * t.z) + (x * t.c) + (y * t.d) - (b * t.e),
                (a * t.b) + (z * t.e),
                (a * t.c) - (y * t.e),
                (a * t.d) + (x * t.e),
                (a * t.e)};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 3:
                    return { (x * t[X]) + (y * t[Y]) + (z * t[Z]),
                                (a * t[X]),
                                (a * t[Y]),
                                (a * t[Z])};
                case 2:
                    return { (x * t[X]) + (y * t[Y]),
                                (a * t[X]),
                                (a * t[Y])};
                case 1:
                    return { (x * t[X]),
                                (a * t[X])};
                default:
                    throw std::invalid_argument("Invalid size for left contraction, maximum size is 3.");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for left contraction");
        }
    }
    template Versor Versor::lco<Versor>(const Versor &t) const;
    template Versor Versor::lco<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::lco<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::lco<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::lco<double>(const double &t) const;
    template Versor Versor::lco<float>(const float &t) const;
    template Versor Versor::lco<int>(const int &t) const;

    //--------------------Right Contraction--------------------
    template <typename T>
    Versor Versor::rco(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {(a * t.a) + (x * t.x) + (y * t.y) - (b * t.b) - (c * t.c) - (d * t.d) - (e * t.e),
                (x * t.a) + (b * t.y) + (c * t.z) - (e * t.d),
                (y * t.a) - (b * t.x) + (d * t.z) + (e * t.c),
                (z * t.a) - (c * t.x) - (d * t.y) - (e * t.b),
                (b * t.a) + (e * t.z),
                (c * t.a) - (e * t.y),
                (d * t.a) + (e * t.x),
                (e * t.a)};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, x * t, y * t, z * t, b * t, c * t, d * t, e * t};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 3:
                    return { (x * t[X]) + (y * t[Y]) + (z * t[Z]),
                                (b * t[Y]) + (c * t[Z]),
                                -(b * t[X]) + (d * t[Z]),
                                -(c * t[X]) - (d * t[Y]),
                                (e * t[Z]),
                                -(e * t[Y]),
                                (e * t[X])};
                case 2:
                    return { (x * t[X]) + (y * t[Y]),
                                (b * t[Y]),
                                -(b * t[X]),
                                -(c * t[X]) - (d * t[Y]),
                                0,
                                -(e * t[Y]),
                                (e * t[X])};
                case 1:
                    return { (x * t[X]),
                                0,
                                -(b * t[X]),
                                -(c * t[X]),
                                0,
                                0,
                                (e * t[X])};
                default:
                    throw std::invalid_argument("Invalid size for right contraction, maximum size is 3.");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for right contraction.");
        }
    }
    template Versor Versor::rco<Versor>(const Versor &t) const;
    template Versor Versor::rco<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::rco<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::rco<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::rco<double>(const double &t) const;
    template Versor Versor::rco<float>(const float &t) const;
    template Versor Versor::rco<int>(const int &t) const;

    //--------------------Inner Product--------------------
    template <typename T>
    Versor Versor::inr(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {(a * t.a) + (x * t.x) + (y * t.y) + (z * t.z) - (b * t.b) - (c * t.c) - (d * t.d) - (e * t.e),
                    (x * t.a) + (a * t.x) + (b * t.y) + (c * t.z) - (y * t.b) - (z * t.c) - (e * t.d) - (d * t.e),
                    (y * t.a) - (b * t.x) + (a * t.y) + (d * t.z) + (x * t.b) + (e * t.c) - (z * t.d) + (c * t.e),
                    (z * t.a) - (c * t.x) - (d * t.y) + (a * t.z) - (e * t.b) + (x * t.c) + (y * t.d) - (b * t.e),
                    (b * t.a) + (e * t.z) + (a * t.b) + (z * t.e),
                    (c * t.a) - (e * t.y) + (a * t.c) - (y * t.e),
                    (d * t.a) + (e * t.x) + (a * t.d) + (x * t.e),
                    (e * t.a) + (a * t.e)};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, x * t, y * t, z * t, b * t, c * t, d * t, e * t};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 3:
                    return { (x * t[X]) + (y * t[Y]) + (z * t[Z]),
                                (a * t[X]) + (b * t[Y]) + (c * t[Z]),
                                -(b * t[X]) + (a * t[Y]) + (d * t[Z]),
                                -(c * t[X]) - (d * t[Y]) + (a * t[Z]),
                                (e * t[Z]),
                                -(e * t[Y]),
                                (e * t[X])};
                case 2:
                    return { (x * t[X]) + (y * t[Y]),
                                (a * t[X]) + (b * t[Y]),
                                -(b * t[X]) + (a * t[Y]),
                                -(c * t[X]) - (d * t[Y]),
                                0,
                                -(e * t[Y]),
                                (e * t[X])};
                case 1:
                    return { (x * t[X]),
                                (a * t[X]),
                                -(b * t[X]),
                                -(c * t[X]),
                                0,
                                0,
                                (e * t[X])};
                default:
                    throw std::invalid_argument("Invalid size for inner product. Maximum size is 3");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for inner product");
        }
    }
    template Versor Versor::inr<Versor>(const Versor &t) const;
    template Versor Versor::inr<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::inr<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::inr<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::inr<double>(const double &t) const;
    template Versor Versor::inr<float>(const float &t) const;
    template Versor Versor::inr<int>(const int &t) const;

    //--------------------Outer Product--------------------
    template <typename T>
    Versor Versor::ext(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {(a * t.a),
                    (x * t.a) + (a * t.x),
                    (y * t.a) + (a * t.y),
                    (z * t.a) + (a * t.z),
                    (b * t.a) - (y * t.x) + (x * t.y) + (a * t.b),
                    (c * t.a) - (z * t.x) + (x * t.z) + (a * t.c),
                    (d * t.a) - (z * t.y) + (y * t.z) + (a * t.d),
                    (e * t.a) + (d * t.x) - (c * t.y) + (b + t.z) + (z * t.b) - (y * c) + (x * t.d) + (a * t.e)};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, x * t, y * t, z * t, b * t, c * t, d * t, e * t};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 3:
                    return { 0,
                                (a * t[X]),
                                (a * t[Y]),
                                (a * t[Z]),
                                -(y * t[X]) + (x * t[Y]),
                                -(z * t[X]) + (x * t[Z]),
                                -(z * t[Y]) + (y * t[Z]),
                                (d * t[X]) - (c * t[Y]) + (b * t[Z])};
                case 2:
                    return { 0,
                                (a * t[X]),
                                (a * t[Y]),
                                0,
                                -(y * t[X]) + (x * t[Y]),
                                -(z * t[X]),
                                -(z * t[Y]),
                                (d * t[X]) - (c * t[Y])};
                case 1:
                    return { 0,
                                (a * t[X]),
                                0,
                                0,
                                -(y * t[X]),
                                -(z * t[X]),
                                0,
                                (d * t[X])};
                default:
                    throw std::invalid_argument("Invalid size for outer product, maximum size is 3.");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for outer product.");
        }
    }
    template Versor Versor::ext<Versor>(const Versor &t) const;
    template Versor Versor::ext<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::ext<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::ext<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::ext<double>(const double &t) const;
    template Versor Versor::ext<float>(const float &t) const;
    template Versor Versor::ext<int>(const int &t) const;

    //--------------------Geometric Product--------------------
    template <typename T>
    Versor Versor::mul(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return { (a * t.a) + (x * t.x) + (y * t.y) + (z * t.z) - (b * t.b) - (c * t.c) - (d * t.d) - (e * t.e),
                        (x * t.a) + (a * t.x) + (b * t.y) + (c * t.z) - (y * t.b) - (z * t.c) - (e * t.d) - (d * t.e),
                        (y * t.a) - (b * t.x) + (a * t.y) + (d * t.z) + (x * t.b) + (e * t.c) - (z * t.d) + (c * t.e),
                        (z * t.a) - (c * t.x) - (d * t.y) + (a * t.z) - (e * t.b) + (x * t.c) + (y * t.d) - (b * t.e),
                        (b * t.a) - (y * t.x) + (x * t.y) + (e * t.z) + (a * t.b) + (d * t.c) - (c * t.d) + (z * t.e),
                        (c * t.a) - (z * t.x) - (e * t.y) + (x * t.z) - (d * t.b) + (a * t.c) + (b * t.d) - (y * t.e),
                        (d * t.a) + (e * t.x) - (z * t.y) + (y * t.z) + (c * t.b) - (b * t.c) + (a * t.d) + (x * t.e),
                        (e * t.a) + (d * t.x) - (c * t.y) + (b * t.z) + (z * t.b) - (y * t.c) + (x * t.d) + (a * t.e)};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, x * t, y * t, z * t, b * t, c * t, d * t, e * t};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 3:
                    return { (x * t[X]) + (y * t[Y]) + (z * t[Z]),
                                (a * t[X]) + (b * t[Y]) + (c * t[Z]),
                                -(b * t[X]) + (a * t[Y]) + (d * t[Z]),
                                -(c * t[X]) - (d * t[Y]) + (a * t[Z]),
                                -(y * t[X]) + (x * t[Y]) + (e * t[Z]),
                                -(z * t[X]) - (e * t[Y]) + (x * t[Z]),
                                -(c * t[X]) - (z * t[Y]) + (y * t[Z]),
                                (d * t[X]) - (c * t[Y]) + (b * t[Z])};
                case 2:
                    return { (x * t[X]) + (y * t[Y]),
                                (a * t[X]) + (b * t[Y]),
                                -(b * t[X]) + (a * t[Y]),
                                -(c * t[X]) - (d * t[Y]),
                                -(y * t[X]) + (x * t[Y]),
                                -(z * t[X]) - (e * t[Y]),
                                -(c * t[X]) - (z * t[Y]),
                                (d * t[X]) - (c * t[Y])};
                case 1:
                    return { (x * t[X]),
                                (a * t[X]),
                                -(b * t[X]),
                                -(c * t[X]),
                                -(y * t[X]),
                                -(z * t[X]),
                                -(c * t[X]),
                                (d * t[X])};
                default:
                    throw std::invalid_argument("Invalid size for outer product, maximum size is 3.");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for geometric product.");
        }
    }
    template Versor Versor::mul<Versor>(const Versor &t) const;
    template Versor Versor::mul<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::mul<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::mul<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::mul<double>(const double &t) const;
    template Versor Versor::mul<float>(const float &t) const;
    template Versor Versor::mul<int>(const int &t) const;
    Versor Versor::negate() const
    {
        return {a, x, y, z, -b, -c, -d, -e};
    }

    //--------------------Division----------------------
    template <typename T>
    Versor Versor::div(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return this->mul(t.inverse());
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return this->sdiv(t);
        }
        else {
            throw std::invalid_argument("Invalid type for division");
        }
    }
    // Explicit instantiation of the div method template for Versor
    template Versor Versor::div<Versor>(const Versor &t) const;
    template Versor Versor::div<double>(const double &t) const;
    template Versor Versor::div<float>(const float &t) const;
    template Versor Versor::div<int>(const int &t) const;

    //--------------------Division----------------------
    template <typename T>
    Versor Versor::sdiv(const T &t) const
    {
        if constexpr (std::is_arithmetic_v<T>) {
            return {a / t, x / t, y / t, z / t, b / t, c / t, d / t, e / t};
        }
        else {
            throw std::invalid_argument("Invalid type for division.");
        }
    }
    // Explicit instantiation of the div method template for Versor
    template Versor Versor::sdiv<double>(const double &t) const;
    template Versor Versor::sdiv<float>(const float &t) const;
    template Versor Versor::sdiv<int>(const int &t) const;

    //--------------------IO--------------------
    std::string Versor::toString() const {
        std::string result;
        if (a != 0) result += std::to_string(a);
        if (x != 0) {
            if (!result.empty()) result += " + ";
            result += std::to_string(x) + "e1";
        }
        if (y != 0) {
            if (!result.empty()) result += " + ";
            result += std::to_string(y) + "e2";
        }
        if (b != 0) {
            if (!result.empty()) result += " + ";
            result += std::to_string(b) + "e1^e2";
        }
        return result.empty() ? "0" : result;
    }
} // VRSR