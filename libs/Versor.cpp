//
// Created by zerth on 12/12/24.
//

#include "Versor.h"

namespace VRSR {

    enum VectorComponents { X = 0, Y = 1, Z = 2, W = 3 };

    //These are the basis vectors used for the space. All multivectors are combinations of e1 and e2.
    static const std::vector<float> e1 = {1.0f,0.0f};
    static const std::vector<float> e2 = {0.0f,1.0f};

    //--------------------Addition--------------------
    // ADDITION of multivectors is simply the individual sums of the compononents.
    template <typename T>
    Versor Versor::add(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a + t.a, x + t.x, y + t.y, b + t.b};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a + t, x, y, b};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {a + t[0], x + t[1], y + t[2], b + t[3]};
                case 3:
                    return {a + t[0], x + t[1], y + t[2], b};
                case 2:
                    return {a + t[0], x + t[1], y, b};
                case 1:
                    return {a + t[0], x, y, b};
                default:
                    throw std::invalid_argument("Invalid size for addition");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for addition");
        }
    }
    // Explicit instantiation of the add method template for Versor
    template Versor Versor::add<Versor>(const Versor &t) const;
    template Versor Versor::add<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::add<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::add<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::add<double>(const double &t) const;
    template Versor Versor::add<float>(const float &t) const;
    template Versor Versor::add<int>(const int &t) const;

    //--------------------Subtraction--------------------
    // SUBTRACTION is defined as the sum of the inverse of the input versor.
    // A - B = A + (-B), here is the literal implementation of this.
    template <typename T>
    Versor Versor::sub(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a - t.a, x - t.x, y - t.y, b - t.b};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a - t, x, y, b};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {a - t[0], x - t[1], y - t[2], b - t[3]};
                case 3:
                    return {a - t[0], x - t[1], y - t[2], b};
                case 2:
                    return {a - t[0], x - t[1], y, b};
                case 1:
                    return {a - t[0], x, y, b};
                default:
                    throw std::invalid_argument("Invalid size for subtraction");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for subtraction");
        }
    }
    // Explicit instantiation of the sub method template for Versor
    template Versor Versor::sub<Versor>(const Versor &t) const;
    template Versor Versor::sub<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::sub<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::sub<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::sub<double>(const double &t) const;
    template Versor Versor::sub<float>(const float &t) const;
    template Versor Versor::sub<int>(const int &t) const;

    //--------------------Left Contraction--------------------
    // Contraction can be thought of the act of removing one object from another and leaving the result.
    // It has many elagent uses in physics and mathematics.
    // Given the two orthogonal basis vectores e1 and e2, contraction can be described by the following:
    // e1 << e1 = 1, e1 << e2 = 0
    // e1 << (e1 ^ e2) = e2, e1 << (e2 ^ e1) = -e2
    // (e1 ^ e2) << (e1 ^ e2) = e1 << (e2 << (e1 ^ e2)) = -1
    // a << rhs = arhs
    // lhs << a = 0 if lhs grade > 0
    template <typename T>
    Versor Versor::lco(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a * t.a + x * t.x + y * t.y - b * t.b, 0.0f, 0.0f, 0.0f};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, 0.0f, 0.0f, 0.0f};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {a * t[0] + x * t[1] + y * t[2] - b * t[3], 0.0f, 0.0f, 0.0f};
                case 3:
                    return {a * t[0] + x * t[1] + y * t[2], 0.0f, 0.0f, 0.0f};
                case 2:
                    return {a * t[0] + x * t[1], 0.0f, 0.0f, 0.0f};
                case 1:
                    return {a * t[0], 0.0f, 0.0f, 0.0f};
                default:
                    throw std::invalid_argument("Invalid size for left contraction");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for left contraction");
        }
    }
    // Explicit instantiation of the lco method template for Versor
    /*
    template Versor Versor::lco<Versor>(const Versor &t) const;
    template Versor Versor::lco<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::lco<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::lco<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::lco<double>(const double &t) const;
    template Versor Versor::lco<float>(const float &t) const;
    template Versor Versor::lco<int>(const int &t) const;
    */

    //--------------------Right Contraction--------------------
    template <typename T>
    Versor Versor::rco(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a * t.a + x * t.x + y * t.y - b * t.b, 0.0f, 0.0f, 0.0f};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, 0.0f, 0.0f, 0.0f};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {a * t[0] + x * t[1] + y * t[2] - b * t[3], 0.0f, 0.0f, 0.0f};
                case 3:
                    return {a * t[0] + x * t[1] + y * t[2], 0.0f, 0.0f, 0.0f};
                case 2:
                    return {a * t[0] + x * t[1], 0.0f, 0.0f, 0.0f};
                case 1:
                    return {a * t[0], 0.0f, 0.0f, 0.0f};
                default:
                    throw std::invalid_argument("Invalid size for right contraction");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for right contraction");
        }
    }
    // Explicit instantiation of the rco method template for Versor
    /*
    template Versor Versor::rco<Versor>(const Versor &t) const;
    template Versor Versor::rco<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::rco<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::rco<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::rco<double>(const double &t) const;
    template Versor Versor::rco<float>(const float &t) const;
    template Versor Versor::rco<int>(const int &t) const;
    */

    //--------------------Exterior Product--------------------
    // Wedge Product, used to create bivectors and trivectors out of vectors. Also known as the exterior product.
    // This quantity represents the oriented area or volume constructed by the components.
    // A is a versor containing a scalar and vector components; a+xe1+ye2
    // A ^ B = (a_x * e1)(a_y * e2) ^ (b_x * e1)(b_y * e2), expanding with FOIL
    // (a_x * e1) ^ (b_x * e1) + (a_x * e1) ^ (b_y * e2) + (a_y * e2) ^ (b_x * e1) + (a_y * e2) ^ (b_y * e2) =
    // since e1_a = e1_b and a ^ a = 0, we can simplify this to:
    // 0 + (a_x * e1) ^ (b_y * e2) + (a_y * e2) ^ (b_x * e1) + 0 =
    // Using a ^ b = -b ^ a, we can simplify this to:
    // (a_x * e1) ^ (b_y * e2) - (a_y * e1) ^ (b_x * e2) =
    // Through factoring a_x, a_y, b_x and b_y, we can simplify this to:
    // (a_x * b_y - a_y * b_x) * e1 ^ e2
    // e1 ^ e2 is the unit bivector that represents the oriented volume of our basis vectors.
    // The rest is a scalar quantity.
    // If the wedge product is zero, then the two vectors are parallel.
    template <typename T>
    Versor Versor::ext(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            float tempa = (x * t.y) - (y * t.x);
            float tempx = (b * t.y) - (a * t.x);
            float tempy = (a * t.y) - (b * t.x);
            float tempb = (x * t.y) - (y * t.x);
            return { tempa, tempx, tempy, tempb };
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {0.0f, a * t, x * t, y * t};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {0.0f, x * t[3], y * t[3], b * t[3]};
                case 3:
                    return {0.0f, x * t[2], y * t[2], b};
                case 2:
                    return {0.0f, x * t[1], y, b};
                case 1:
                    return {0.0f, x * t[0], y, b};
                default:
                    throw std::invalid_argument("Invalid size for exterior product");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for exterior product");
        }
    }
    // Explicit instantiation of the ext method template for Versor
    /*
    template Versor Versor::ext<Versor>(const Versor &t) const;
    template Versor Versor::ext<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::ext<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::ext<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::ext<double>(const double &t) const;
    template Versor Versor::ext<float>(const float &t) const;
    template Versor Versor::ext<int>(const int &t) const;
    */

    //--------------------Interior Product--------------------
    // Dot Product, used to create scalars out of vectors. Also known as the inner product or scalar product.
    // Represents the projection of one vector onto another multiplied by the product of their magnitudes.
    // A . B = (a_x * e1 + a_y * e2) . (b_x * e1 + b_y * e2), expanding with FOIL
    // (a_x * e1) . (b_x * e1) + (a_x * e1) . (b_y * e2) + (a_y * e2) . (b_x * e1) + (a_y * e2) . (b_y * e2) =
    // Since e1 and e2 are orthogonal, e1 . e2 = 0 we can simplify this to:
    // (a_x * e1) . (b_x * e1) + 0 + 0 + (a_y * e2) . (b_y * e2) =
    // Through factoring a_x, a_y, b_x and b_y, we can simplify this to:
    // (a_x * b_x)(e1 . e1) + (a_y * b_y)(e2 . e2) =
    // Since e1 . e1 = 1 and e2 . e2 = 1, we can simplify this to:
    // (a_x * b_x) + (a_y * b_y), just a scalar quantity.
    // The inner product of itself is the square of the magnitude of the vector. A . A = |A|^2
    template <typename T>
    Versor Versor::inr(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a * t.a + x * t.x + y * t.y, 0.0f, 0.0f, 0.0f};
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, 0.0f, 0.0f, 0.0f};
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {a * t[0] + x * t[1] + y * t[2], 0.0f, 0.0f, 0.0f};
                case 3:
                    return {a * t[0] + x * t[1] + y * t[2], 0.0f, 0.0f, 0.0f};
                case 2:
                    return {a * t[0] + x * t[1], 0.0f, 0.0f, 0.0f};
                case 1:
                    return {a * t[0], 0.0f, 0.0f, 0.0f};
                default:
                    throw std::invalid_argument("Invalid size for interior product");
            }
        }
        else {
            throw std::invalid_argument("Invalid type for interior product");
        }
    }
    // Explicit instantiation of the inr method template for Versor
    /*
    template Versor Versor::inr<Versor>(const Versor &t) const;
    template Versor Versor::inr<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::inr<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::inr<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::inr<double>(const double &t) const;
    template Versor Versor::inr<float>(const float &t) const;
    template Versor Versor::inr<int>(const int &t) const;
    */

    //--------------------Multiplication--------------------
    // Versor multiplication is not standard multiplication, it is what is known as a geometric product.
    // If we start by assuming one can multiply two multivectors together, we can define the geometric product through the following:
    // ab = (ab + ab) / 2 = (ab + ba + ab - ba) / 2
    // (ab+ba)/2 + (ab-ba)/2
    // (ab+ba)/2 is the dot product, and (ab-ba)/2 is the wedge product.
    // The dot product is the projection of one vector onto another multiplied by the product of their magnitudes.
    // The wedge product is the oriented area or volume constructed by the components.
    // The geometric product is the sum of the dot and wedge products.
    // The geometric product is associative, distributive and the wedge component anti commutative.
    // e1e2 = e1 ^ e2, the geometric product of two orthogonal basis vectors is the wedge product.
    // e1e1 = 1
    // e1(e1e2) = (e1e1)e2 = e2s
    template <typename T>
    Versor Versor::mul(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            float tempa = (a * t.a) + (x * t.x) + (y * t.y) + (-b * t.b);
            float tempx = (x * t.a) + (b * t.y) + (a * t.x) + (-y * t.b);
            float tempy = (a * t.y) + (x * t.b) + (y * t.a) + (-b * t.x);
            float tempb = (b * t.a) + (a * t.b) + (x * t.y) + (-y * t.x);
            return { tempa, tempx, tempy, tempb };
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {a * t[0], x * t[1], y * t[2], b * t[3]};
                case 3:
                    return {a * t[0], x * t[1], y * t[2], b};
                case 2:
                    return {a * t[0], x * t[1], y, b};
                case 1:
                    return {a * t[0], x, y, b};
                default:
                    throw std::invalid_argument("Invalid size for multiplication");
            }
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a * t, x * t, y * t, b * t};
        }
        else {
            throw std::invalid_argument("Invalid type for multiplication");
        }
    }
    // Explicit instantiation of the mul method template for Versor
    /*
    template Versor Versor::mul<Versor>(const Versor &t) const;
    template Versor Versor::mul<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::mul<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::mul<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::mul<double>(const double &t) const;
    template Versor Versor::mul<float>(const float &t) const;
    template Versor Versor::mul<int>(const int &t) const;
    */

    //--------------------Division----------------------
    // Division with multivectors while simple at first is complicated.
    // Typically only division of scalars is possible. Each component of the multivector must be divided by the scalar.
    // This is not a true division, but a scaling of the multivector.
    // Traditionally it is thought you can not divide a vector by another vector, however this is not true.
    template <typename T>
    Versor Versor::div(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            if (t.a != 0.0f) {
                return {a / t.a, x / t.a, y / t.a, b / t.a};    //If the input versor is a scalar
            } else if (t.b == 0.0f) {
                return {0.0f, x / t.x, y / t.y, 0.0f};          //If the input versor is a vector
            } else {
                return {0.0f, 0.0f, 0.0f, b / t.b};             //If the input versor is a bivector
            }
        }
        else if constexpr (is_vector_v<T> && std::is_arithmetic_v<typename T::value_type>) {
            switch (t.size()) {
                case 4:
                    return {a / t[0], x / t[1], y / t[2], b / t[3]};
                case 3:
                    return {a / t[0], x / t[1], y / t[2], b};
                case 2:
                    return {a / t[0], x / t[1], y, b};
                case 1:
                    return {a / t[0], x, y, b};
                default:
                    throw std::invalid_argument("Invalid size for division");
            }
        }
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a / t, x / t, y / t, b / t};
        }
        else {
            throw std::invalid_argument("Invalid type for division");
        }
    }
    // Explicit instantiation of the div method template for Versor
    /*
    template Versor Versor::div<Versor>(const Versor &t) const;
    template Versor Versor::div<std::vector<double>>(const std::vector<double> &t) const;
    template Versor Versor::div<std::vector<float>>(const std::vector<float> &t) const;
    template Versor Versor::div<std::vector<int>>(const std::vector<int> &t) const;
    template Versor Versor::div<double>(const double &t) const;
    template Versor Versor::div<float>(const float &t) const;
    template Versor Versor::div<int>(const int &t) const;
    */

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