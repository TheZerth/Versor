//
// Created by zerth on 12/12/24.
//

#include "Versor.h"

namespace VRSR {

    //These are the basis vectors used for the space. All multivectors are combinations of e1 and e2.
    static const std::vector<float> e1 = {1.0f,0.0f};
    static const std::vector<float> e2 = {0.0f,1.0f};

    //Versor Addition
    Versor Versor::add(const Versor &v) const { return                  {a + v.a,   x + v.x,    y + v.y,    b + v.b}; }
    Versor Versor::add(const std::vector<double> &d) const { return     {a + d[0],  x + d[1],   y + d[2],   b + d[3]}; }
    Versor Versor::add(const std::vector<float> &f) const { return      {a + f[0],  x + f[1],   y + f[2],   b + f[3]}; }
    Versor Versor::add(const std::vector<int> &i) const { return        {a + i[0],  x + i[1],   y + i[2],   b + i[3]}; }
    Versor Versor::add(const double &d) const { return                  {a + d,     x,          y,          b}; }
    Versor Versor::add(const float &f) const { return                   {a + f,     x,          y,          b}; }
    Versor Versor::add(const int &i) const { return                     {a + i,     x,          y,          b}; }

    //Versor Subtraction
    Versor Versor::sub(const Versor &v) const { return                  {a - v.a,   x - v.x,    y - v.y,    b - v.b}; }
    Versor Versor::sub(const std::vector<double> &d) const { return     {a - d[0],  x - d[1],   y - d[2],   b - d[3]}; }
    Versor Versor::sub(const std::vector<float> &f) const { return      {a - f[0],  x - f[1],   y - f[2],   b - f[3]}; }
    Versor Versor::sub(const std::vector<int> &i) const { return        {a - i[0],  x - i[1],   y - i[2],   b - i[3]}; }
    Versor Versor::sub(const double &d) const { return                  {a - d,     x,          y,          b}; }
    Versor Versor::sub(const float &f) const { return                   {a - f,     x,          y,          b}; }
    Versor Versor::sub(const int &i) const { return                     {a - i,     x,          y,          b}; }

} // VRSR