#ifndef TRIANGLE_CLASS_H
#define TRIANGLE_CLASS_H

#include<iostream>
#include"v3.h"
using namespace std;

class Triangle{

  public:
  
    Triangle(v3 n, v3 v0, v3 v1, v3 v2);

    void setZMin (float z) ;

    void setZMax (float z) ;

    // Triangle& operator-=(const v3 &pt) ;

    // bool operator <(const Triangle &t) ;

    // friend ostream& operator<<(ostream& os, const Triangle& t);
  
  public:

    v3 v[3];
    v3 normal;
    float zMin;
    float zMax;
};  

ostream& operator<<(ostream& os, const Triangle& t);

#endif