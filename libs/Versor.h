//
// Created by zerth on 12/12/24.
//

#ifndef VERSOR_H
#define VERSOR_H
#include <iomanip>
#include <ios>
#include <string>
#include <vector>

namespace VRSR {

class Versor {
//--------------------Variables--------------------
public:
    std::vector<double> e1, e2;    //These are the basis vectors.
    double a; //This represents the coefficient to the unit scalar, 1.
    double x; //This represents the coefficient to the basis vector, e1.
    double y; //This represents the coefficient to the basis vector, e2.
    double b; //This represents the oriented area or volume of the bivector. Will be zero for vectors and scalars.
public:
    Versor(const double a, const double x, const double y, const double b) : a(a), x(x), y(y), b(b) {}
    ~Versor() = default;
//--------------------OPERATORS--------------------
public:
    //--------------------Addition--------------------
    // ADDITION of multivectors is simply the individual sums of the compononents.
    Versor operator+(const Versor &v) const {
        return add(v);
    }
    Versor operator+(const std::vector<double> &d) const {
        return add(d);
    }
    Versor operator+(const std::vector<float> &f) const {
        return add(f);
    }
    Versor operator+(const std::vector<int> &i) const {
        return add(i);
    }
    Versor operator+(const double &d) const {
        return add(d);
    }
    Versor operator+(const float &f) const {
        return add(f);
    }
    Versor operator+(const int &i) const {
        return add(i);
    }

    //--------------------Subtraction--------------------
    // SUBTRACTION is defined as the sum of the inverse of the input versor.
    // A - B = A + (-B), here is the litteral implementation of this.
    // I will not use this due to the fact that it is not as efficient as it could be.
    /*Versor operator-(const Versor& v) const {
        Versor temp{-v.a, -v.x, -v.y, -v.b};
        return *this + temp;
    }*/
    Versor operator-(const Versor &v) const {
        return sub(v);
    }
    Versor operator-(const std::vector<double> &v) const {
        return sub(v);
    }
    Versor operator-(const std::vector<float> &v) const {
        return sub(v);
    }
    Versor operator-(const std::vector<int> &v) const {
        return sub(v);
    }
    Versor operator-(const double &d) const {
        return sub(d);
    }
    Versor operator-(const float &f) const {
        return sub(f);
    }
    Versor operator-(const int &i) const {
        return sub(i);
    }

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
    //--------------------Multiplication--------------------
    // Versor Multiplication (Multivector Multiplication)
    Versor operator*(const Versor &v) const {
        return mul(v);
    }
    // Vector Multiplication
    Versor operator*(const std::vector<double> &d) const {
        return mul(d);
    }
    Versor operator*(const std::vector<float> &f) const {
        return mul(f);
    }
    Versor operator*(const std::vector<int> &i) const {
        return mul(i);
    }
    // Scalar Multiplication
    Versor operator*(const double &d) const {
        return mul(d);
    }
    Versor operator*(const float &f) const {
        return mul(f);
    }
    Versor operator*(const int &i) const {
        return mul(i);
    }

    // Division with multivectors while simple at first is complicated.
    // Typically only division of scalars is possible. Each component of the multivector must be divided by the scalar.
    // This is not a true division, but a scaling of the multivector.
    // Traditionally it is thought you can not divide a vector by another vector, however this is not true.
    //--------------------Division--------------------
    // Vector Division
    Versor operator/(const std::vector<double> &d) const {
        return div(d);
    }
    Versor operator/(const std::vector<float> &f) const {
        return div(f);
    }
    Versor operator/(const std::vector<int> &i) const {
        return div(i);
    }
    // Versor Division
    Versor operator/(const Versor &v) const {
        return div(v);
    }
    // Scalar Division
    Versor operator/(const double &d) const {
        return div(d);
    }
    Versor operator/(const float &f) const {
        return div(f);
    }
    Versor operator/(const int &i) const {
        return div(i);
    }

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
    Versor operator^(const Versor &v) const {
        return ext(v);
    }
    Versor operator^(const std::vector<double> &v) const {
        return ext(v);
    }
    Versor operator^(const std::vector<float> &v) const {
        return ext(v);
    }
    Versor operator^(const std::vector<int> &v) const {
        return ext(v);
    }
    Versor operator^(const double &d) const {
        return ext(d);
    }
    Versor operator^(const float &f) const {
        return ext(f);
    }
    Versor operator^(const int &i) const {
        return ext(i);
    }

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
    Versor operator|(const Versor &v) const {
        return inr(v);
    }
    Versor operator|(const std::vector<double> &v) const {
        return inr(v);
    }
    Versor operator|(const std::vector<float> &v) const {
        return inr(v);
    }
    Versor operator|(const std::vector<int> &v) const {
        return inr(v);
    }
    Versor operator|(const double &d) const {
        return inr(d);
    }
    Versor operator|(const float &f) const {
        return inr(f);
    }
    Versor operator|(const int &i) const {
        return inr(i);
    }

    // Contraction can be thought of the act of removing one object from another and leaving the result.
    // It has many elagent uses in physics and mathematics.
    // Given the two orthogonal basis vectores e1 and e2, contraction can be described by the following:
    // e1 << e1 = 1, e1 << e2 = 0
    // e1 << (e1 ^ e2) = e2, e1 << (e2 ^ e1) = -e2
    // (e1 ^ e2) << (e1 ^ e2) = e1 << (e2 << (e1 ^ e2)) = -1
    // a << rhs = arhs
    // lhs << a = 0 if lhs grade > 0
    //--------------------Left Contraction--------------------
    Versor operator<<(const Versor &v) const {
        return (v.lco(*this));
    }
    Versor operator<<(const std::vector<double> &v) const {
        return (v.lco(*this));
    }
    Versor operator<<(const std::vector<float> &v) const {
        return (v.lco(*this));
    }
    Versor operator<<(const std::vector<int> &v) const {
        return (v.lco(*this));
    }
    Versor operator<<(const double &d) const {
        return (d.lco(*this));
    }
    Versor operator<<(const float &f) const {
        return (f.lco(*this));
    }
    Versor operator<<(const int &i) const {
        return (i.lco(*this));
    }

    //--------------------Right Contraction--------------------
    Versor operator>>(const Versor &v) const {
        return (v.rco(*this));
    }
    Versor operator>>(const std::vector<double> &v) const {
        return (v.rco(*this));
    }
    Versor operator>>(const std::vector<float> &v) const {
        return (v.rco(*this));
    }
    Versor operator>>(const std::vector<int> &v) const {
        return (v.rco(*this));
    }
    Versor operator>>(const double &d) const {
        return (d.rco(*this));
    }
    Versor operator>>(const float &f) const {
        return (f.rco(*this));
    }
    Versor operator>>(const int &i) const {
        return (i.rco(*this));
    }

    //--------------------Negate--------------------
    Versor operator!() const {
        return negate();
    }

    //--------------------Reverse--------------------
    Versor operator~() const {
        return reverse();
    }

    //--------------------IO-STREAM-FUNCTIONS--------------------
    // Console Outputd
    friend std::ostream &operator<<(std::ostream &os, const Versor &v) {
        os << std::fixed << std::setprecision(3) << std::setfill('0')
           << "[" << v.toString() << "]";
        return os;
    }
    // Console Input, create a Versor using A X Y Z input.
    friend std::istream &operator>>(std::istream &is, Versor &v) {
        is >> v.a >> v.x >> v.y >> v.b;
        return is;
    }
    // Return if two Versors are identical.
    bool operator==(const Versor & versor) const {
        return a == versor.a && x == versor.x && y == versor.y && b == versor.b;
    }
//--------------------FUNCTIONS--------------------
public:
    Versor normalize() const;
    Versor sqNorm() const;
    Versor negate() const;
    Versor reverse() const;
    std::string toString() const;
    std::string toLatex() const;
private:
    Versor add(const Versor &v) const;
    Versor add(const std::vector<double> &v) const;
    Versor add(const std::vector<float> &v) const;
    Versor add(const std::vector<int> &v) const;
    Versor add(const double &d) const;
    Versor add(const float &f) const;
    Versor add(const int &i) const;
    Versor sub(const Versor &v) const;
    Versor sub(const std::vector<double> &v) const;
    Versor sub(const std::vector<float> &v) const;
    Versor sub(const std::vector<int> &v) const;
    Versor sub(const double &d) const;
    Versor sub(const float &f) const;
    Versor sub(const int &i) const;
    Versor mul(const Versor &v) const;
    Versor mul(const std::vector<double> &v) const;
    Versor mul(const std::vector<float> &v) const;
    Versor mul(const std::vector<int> &v) const;
    Versor mul(const double &d) const;
    Versor mul(const float &f) const;
    Versor mul(const int &i) const;
    Versor div(const Versor &v) const;
    Versor div(const std::vector<double> &v) const;
    Versor div(const std::vector<float> &v) const;
    Versor div(const std::vector<int> &v) const;
    Versor div(const double &d) const;
    Versor div(const float &f) const;
    Versor div(const int &i) const;
    Versor inr(const Versor &v) const;
    Versor inr(const std::vector<double> &v) const;
    Versor inr(const std::vector<float> &v) const;
    Versor inr(const std::vector<int> &v) const;
    Versor inr(const double &d) const;
    Versor inr(const float &f) const;
    Versor inr(const int &i) const;
    Versor ext(const Versor &v) const;
    Versor ext(const std::vector<double> &v) const;
    Versor ext(const std::vector<float> &v) const;
    Versor ext(const std::vector<int> &v) const;
    Versor ext(const double &d) const;
    Versor ext(const float &f) const;
    Versor ext(const int &i) const;
};

} // VRSR

#endif //VERSOR_H
