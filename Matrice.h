#ifndef MATRICE_H
#define MATRICE_H
#include <vector>

using namespace std;

class Matrice
{
    public:
    Matrice(vector<double> Vec1, vector<double> Vec2, vector<double> Vec3, vector<double> Vec4, vector<double> Vec5, vector<double> Vec6, vector<double> Vec7, vector<double> Vec8, vector<double> Vec9, vector<double> Vec10);
    vector<double> getVec1() const;
    vector<double> getVec2() const;
    vector<double> getVec3() const;
    vector<double> getVec4() const;
    vector<double> getVec5() const;
    vector<double> getVec6() const;
    vector<double> getVec7() const;
    vector<double> getVec8() const;
    vector<double> getVec9() const;
    vector<double> getVec10() const;

    private:
    vector<double> m1;
    vector<double> m2;
    vector<double> m3;
    vector<double> m4;
    vector<double> m5;
    vector<double> m6;
    vector<double> m7;
    vector<double> m8;
    vector<double> m9;
    vector<double> m10;
};

#endif // MATRICE_H
