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
    std::vector<double> e1, e2, e3;    //These are the basis vectors.
    double a; //This represents the coefficient to the unit scalar, 1.
    double x; //This represents the coefficient to the basis vector, e1.
    double y; //This represents the coefficient to the basis vector, e2.
    double z; //This represents the coefficient to the basis vector, e3.
    double b; //This represents the oriented area or volume of the bivector e1^e2.
    double c; //This represents the oriented area or volume of the bivector e1^e3.
    double d; //This represents the oriented area or volume of the bivector e2^e3.
    double e; //This represents the oriented area or volume of the trivector e1^e2^e3.
public:
    Versor(const double a = 0, const double x = 0, const double y = 0, const double z = 0,
           const double b = 0, const double c = 0, const double d = 0, const double e = 0)
        : a(a), x(x), y(y), z(z), b(b), c(c), d(d), e(e) {}
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
    bool operator==(const Versor & versor) const    { return a == versor.a &&   x == versor.x &&    y == versor.y &&    z == versor.z && b == versor.b && c == versor.c && d == versor.d && e == versor.e; }
    bool operator!=(const Versor & versor) const    { return a != versor.a ||   x != versor.x ||    y != versor.y ||    z != versor.z || b != versor.b || c != versor.c || d != versor.d || e != versor.e; }
    bool operator>(const Versor & versor) const     { return a > versor.a &&    x > versor.x &&     y > versor.y &&     z > versor.z && b > versor.b && c > versor.c && d > versor.d && e > versor.e; }
    bool operator>=(const Versor & versor) const    { return a >= versor.a &&   x >= versor.x &&    y >= versor.y &&    z >= versor.z && b >= versor.b && c >= versor.c && d >= versor.d && e >= versor.e; }
    bool operator<(const Versor & versor) const     { return a < versor.a &&    x < versor.x &&     y < versor.y &&     z < versor.z && b < versor.b && c < versor.c && d < versor.d && e < versor.e; }
    bool operator<=(const Versor & versor) const    { return a <= versor.a &&   x <= versor.x &&    y <= versor.y &&    z <= versor.z && b <= versor.b && c <= versor.c && d <= versor.d && e <= versor.e; }
    //--------------------IO-STREAM-FUNCTIONS--------------------
    // Console Outputd
    friend std::ostream &operator<<(std::ostream &os, const Versor &v)  { os << std::fixed << std::setprecision(3) << std::setfill('0') << "[" << v.toString() << "]"; return os; }
    // Console Input, create a Versor using A X Y Z input.
    friend std::istream &operator>>(std::istream &is, Versor &v)        { is >> v.a >> v.x >> v.y >> v.z >> v.b >> v.c >> v.d >> v.e; return is; };
//---------------------------------------------------------------
//---------------------------FUNCTIONS---------------------------
    Versor normalize() const;
    double sqNorm() const;
    Versor negate() const;
    Versor reverse() const;
    Versor dual() const;
    Versor conjugate() const;
    Versor involute() const;
    double magnitude() const;
    Versor inverse() const;
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
    template <typename T>
    Versor sdiv(const T &t) const;
};

} // VRSR

#endif //VERSOR_H

