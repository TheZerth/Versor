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

    //--------------------Multiplication-----------------
    template <typename T>
    Versor operator*(const T &t) const { return mul(t); }

    //--------------------Division----------------------
    template <typename T>
    Versor operator/(const T &t) const { return div(t); }

    //--------------------Exterior Product--------------------
    template <typename T>
    Versor operator^(const T &t) const { return ext(t); }

    //--------------------Interior Product--------------------
    template <typename T>
    Versor operator|(const T &t) const { return inr(t); }

    //--------------------Left Contraction--------------------
    template <typename T>
    Versor operator<<(const T &t) const { return lco(t); }

    //--------------------Right Contraction--------------------
    template <typename T>
    Versor operator>>(const T &t) const { return rco(t); }

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
    Versor add(const Versor &v) const;
    Versor add(const std::vector<double> &d) const;
    Versor add(const std::vector<float> &f) const;
    Versor add(const std::vector<int> &i) const;
    Versor add(const double &d) const;
    Versor add(const float &f) const;
    Versor add(const int &i) const;
    Versor sub(const Versor &v) const;
    Versor sub(const std::vector<double> &d) const;
    Versor sub(const std::vector<float> &f) const;
    Versor sub(const std::vector<int> &i) const;
    Versor sub(const double &d) const;
    Versor sub(const float &f) const;
    Versor sub(const int &i) const;
    Versor mul(const Versor &v) const;
    Versor mul(const std::vector<double> &d) const;
    Versor mul(const std::vector<float> &f) const;
    Versor mul(const std::vector<int> &i) const;
    Versor mul(const double &d) const;
    Versor mul(const float &f) const;
    Versor mul(const int &i) const;
    Versor div(const Versor &v) const;
    Versor div(const std::vector<double> &d) const;
    Versor div(const std::vector<float> &f) const;
    Versor div(const std::vector<int> &i) const;
    Versor div(const double &d) const;
    Versor div(const float &f) const;
    Versor div(const int &i) const;
    Versor inr(const Versor &v) const;
    Versor inr(const std::vector<double> &d) const;
    Versor inr(const std::vector<float> &f) const;
    Versor inr(const std::vector<int> &i) const;
    Versor inr(const double &d) const;
    Versor inr(const float &f) const;
    Versor inr(const int &i) const;
    Versor ext(const Versor &v) const;
    Versor ext(const std::vector<double> &d) const;
    Versor ext(const std::vector<float> &f) const;
    Versor ext(const std::vector<int> &i) const;
    Versor ext(const double &d) const;
    Versor ext(const float &f) const;
    Versor ext(const int &i) const;
    Versor lco(const Versor &v) const;
    Versor lco(const std::vector<double> &d) const;
    Versor lco(const std::vector<float> &f) const;
    Versor lco(const std::vector<int> &i) const;
    Versor lco(const double &d) const;
    Versor lco(const float &f) const;
    Versor lco(const int &i) const;
    Versor rco(const Versor &v) const;
    Versor rco(const std::vector<double> &d) const;
    Versor rco(const std::vector<float> &f) const;
    Versor rco(const std::vector<int> &i) const;
    Versor rco(const double &d) const;
    Versor rco(const float &f) const;
    Versor rco(const int &i) const;
};

} // VRSR

#endif //VERSOR_H
