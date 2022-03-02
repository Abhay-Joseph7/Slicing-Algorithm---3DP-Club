#ifndef TRIANGLE_MESH_CLASS_H
#define TRIANGLE_MESH_CLASS_H

#include <vector>
#include<iostream>
#include"TriangleClass.h"

using namespace std;

class TriangleMesh {

  public:

    TriangleMesh();

    size_t size() const ;

    void push_back(Triangle &t) ;

    v3 meshAABBSize() const ;

    const vector<Triangle>& getvTriangle() const;

    v3 getBottomLeftVertex() const ;

    v3 getUpperRightVertex() const ;

  private:

    int meshSize;
    v3 bottomLeftVertex;
    v3 upperRightVertex;
    vector<Triangle> vTriangle;
};

#endif