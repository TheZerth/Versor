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
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a + t, x, y, b};
        }
        else {
            throw std::invalid_argument("Invalid type for addition");
        }
    }

    //--------------------Subtraction--------------------
    // SUBTRACTION is defined as the sum of the inverse of the input versor.
    // A - B = A + (-B), here is the literal implementation of this.
    template <typename T>
    Versor Versor::sub(const T &t) const
    {
        if constexpr (std::is_same_v<T, Versor>) {
            return {a - t.a, x - t.x, y - t.y, b - t.b};
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
        else if constexpr (std::is_arithmetic_v<T>) {
            return {a - t, x, y, b};
        }
        else {
            throw std::invalid_argument("Invalid type for subtraction");
        }
    }
    /*
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
    Versor Versor::mul(const Versor &v) const {
        float tempa = (a * v.a) + (x * v.x) + (y * v.y) + (-b * v.b);
        float tempx = (x * v.a) + (b * v.y) + (a * v.x) + (-y * v.b);
        float tempy = (a * v.y) + (x * v.b) + (y * v.a) + (-b * v.x);
        float tempb = (b * v.a) + (a * v.b) + (x * v.y) + (-y * v.x);
        return { tempa, tempx, tempy, tempb };
    }
    Versor Versor::mul(const std::vector<double> &d) const {
        float tempa = (x * d[X]) + (y * d[Y]);
        float tempx = (b * d[Y]) + (a * d[X]);
        float tempy = (a * d[Y]) + (-b * d[X]);
        float tempb = (x * d[Y]) + (-y * d[X]);
        return { tempa, tempx, tempy, tempb };
    }

    //--------------------Division----------------------
    // Division with multivectors while simple at first is complicated.
    // Typically only division of scalars is possible. Each component of the multivector must be divided by the scalar.
    // This is not a true division, but a scaling of the multivector.
    // Traditionally it is thought you can not divide a vector by another vector, however this is not true.

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

    //--------------------Left Contraction--------------------
    // Contraction can be thought of the act of removing one object from another and leaving the result.
    // It has many elagent uses in physics and mathematics.
    // Given the two orthogonal basis vectores e1 and e2, contraction can be described by the following:
    // e1 << e1 = 1, e1 << e2 = 0
    // e1 << (e1 ^ e2) = e2, e1 << (e2 ^ e1) = -e2
    // (e1 ^ e2) << (e1 ^ e2) = e1 << (e2 << (e1 ^ e2)) = -1
    // a << rhs = arhs
    // lhs << a = 0 if lhs grade > 0
    */
} // VRSR