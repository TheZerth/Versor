//
// Created by zerth on 12/12/24.
//

#ifndef VERSOR_H
#define VERSOR_H
#include <iomanip>
#include <ios>
#include <string>
#include <thread>
#include <vector>
#include <type_traits>

namespace VRSR {

// Helper template to check if a type is a std::vector
template <typename T>
struct is_vector : std::false_type {};
template <typename T, typename A>
struct is_vector<std::vector<T, A>> : std::true_type {};
template <typename T>
inline constexpr bool is_vector_v = is_vector<T>::value;

class Versor {
//-------------------------------------------------
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
//-------------------------------------------------
//--------------------OPERATORS--------------------
    //--------------------Addition----------------------
    template <typename T>
    Versor operator+(const T &t) const { return add(t); }
    //--------------------Subtraction--------------------
    template <typename T>
    Versor operator-(const T &t) const { return sub(t); }
    //--------------------Left Contraction--------------------
    template <typename T>
    Versor operator<<(const T &t) const { return lco(t); }
    //--------------------Right Contraction--------------------
    template <typename T>
    Versor operator>>(const T &t) const { return rco(t); }
    //--------------------Wedge Product--------------------
    template <typename T>
    Versor operator^(const T &t) const { return ext(t); }
    //--------------------Dot Product--------------------
    template <typename T>
    Versor operator|(const T &t) const { return inr(t); }
    //--------------------Multiplication-----------------
    template <typename T>
    Versor operator*(const T &t) const { return mul(t); }
    //--------------------Division----------------------
    template <typename T>
    Versor operator/(const T &t) const { return div(t); }
    //--------------------Negate--------------------
    Versor operator!() const { return negate(); }
    //--------------------Reverse--------------------
    Versor operator~() const { return reverse(); }
    //--------------------Boolean Operations--------------------
    bool operator==(const Versor & versor) const    { return a == versor.a &&   x == versor.x &&    y == versor.y &&    b == versor.b; }
    bool operator!=(const Versor & versor) const    { return a != versor.a ||   x != versor.x ||    y != versor.y ||    b != versor.b; }
    bool operator>(const Versor & versor) const     { return a > versor.a &&    x > versor.x &&     y > versor.y &&     b > versor.b; }
    bool operator>=(const Versor & versor) const    { return a >= versor.a &&   x >= versor.x &&    y >= versor.y &&    b >= versor.b; }
    bool operator<(const Versor & versor) const     { return a < versor.a &&    x < versor.x &&     y < versor.y &&     b < versor.b; }
    bool operator<=(const Versor & versor) const    { return a <= versor.a &&   x <= versor.x &&    y <= versor.y &&    b <= versor.b; }
    //--------------------IO-STREAM-FUNCTIONS--------------------
    // Console Outputd
    friend std::ostream &operator<<(std::ostream &os, const Versor &v)  { os << std::fixed << std::setprecision(3) << std::setfill('0') << "[" << v.toString() << "]"; return os; }
    // Console Input, create a Versor using A X Y Z input.
    friend std::istream &operator>>(std::istream &is, Versor &v)        { is >> v.a >> v.x >> v.y >> v.b; return is; };
//---------------------------------------------------------------
//---------------------------FUNCTIONS---------------------------
    Versor normalize() const;
    Versor sqNorm() const;
    Versor negate() const;
    Versor reverse() const;
    std::string toString() const;
    std::string toLatex() const;
private:
    template <typename T>
    Versor add(const T &t) const;
    template <typename T>
    Versor sub(const T &t) const;
    template <typename T>
    Versor lco(const T &t) const;
    template <typename T>
    Versor rco(const T &t) const;
    template <typename T>
    Versor inr(const T &t) const;
    template <typename T>
    Versor ext(const T &t) const;
    template <typename T>
    Versor mul(const T &t) const;
    template <typename T>
    Versor div(const T &t) const;
};

} // VRSR

#endif //VERSOR_H
