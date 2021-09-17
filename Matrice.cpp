#include "Matrice.h"
#include <vector>

using namespace std;

Matrice::Matrice(vector<double> Vec1, vector<double> Vec2, vector<double> Vec3, vector<double> Vec4, vector<double> Vec5, vector<double> Vec6, vector<double> Vec7, vector<double> Vec8, vector<double> Vec9, vector<double> Vec10):m1(Vec1),m2(Vec2),m3(Vec3),m4(Vec4),m5(Vec5),m6(Vec6),m7(Vec7),m8(Vec8),m9(Vec9),m10(Vec10)
{

}

vector<double> Matrice::getVec1() const
{
    return m1;
}
vector<double> Matrice::getVec2() const
{
    return m2;
}
vector<double> Matrice::getVec3() const
{
    return m3;
}
vector<double> Matrice::getVec4() const
{
    return m4;
}
vector<double> Matrice::getVec5() const
{
    return m5;
}
vector<double> Matrice::getVec6() const
{
    return m6;
}
vector<double> Matrice::getVec7() const
{
    return m7;
}
vector<double> Matrice::getVec8() const
{
    return m8;
}
vector<double> Matrice::getVec9() const
{
    return m9;
}
vector<double> Matrice::getVec10() const
{
    return m10;
}
