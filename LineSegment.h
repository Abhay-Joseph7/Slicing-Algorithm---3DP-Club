#ifndef LINE_SEGMENT_H
#define LINE_SEGMENT_H

#include<iostream>
#include"v3.h"

using namespace std;

class LineSegment {
    
    public:

      LineSegment (v3 p0=v3(), v3 p1=v3(), int i=0);

      // bool operator==(const LineSegment &ls) const ;

      // friend ostream& operator<<(ostream& os, const LineSegment& ls) ;
      
      v3 v[2];

    private:

      double a;

      double b;

      bool vertical;

      int index;

};

#endif