#include<iostream>
#include<fstream>
#include<math.h>
#include<string.h>
#include<assert.h>
#include<unordered_map>
#include<vector>
#include<algorithm>

#include"v3.h"
#include"TriangleClass.h"
#include"TriangleMesh.h"
#include"LineSegment.h"
//#include"Mesh_Triangle_List.h"

using namespace std;

#define MAX_READ_CHUNK 80
#define FACET_READ_SIZE 12

template<typename T> inline void hash_combine(size_t &seed, const T &v) {
  hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct HashV3 {
  size_t operator() (const v3 &v) const {
    size_t h = hash<float>()(v.x);
    hash_combine (h, v.y);
    hash_combine (h, v.z);
    return h;
  }
};

typedef unordered_map<v3, vector<v3>, HashV3> PointMesh;

/*Structures*/
typedef struct node {
  Triangle t;
  struct node *next;
  struct node *prev;
} Mesh_Triangle_Node_t;

typedef struct _list {
   Mesh_Triangle_Node_t *head;
   Mesh_Triangle_Node_t *tail;
} Mesh_Triangle_List_t;

/*-----------------------------------------------------------------------*/
Mesh_Triangle_List_t* Mesh_Triangle_List_create (void) {
  Mesh_Triangle_List_t *L = (Mesh_Triangle_List_t *)malloc(sizeof(Mesh_Triangle_List_t));
  L->head = NULL;
  L->tail = NULL;
  return L;
}

/*-----------------------------------------------------------------------*/
void Mesh_Triangle_List_insert (Triangle t, Mesh_Triangle_List_t *L) {
  Mesh_Triangle_Node_t *node = (Mesh_Triangle_Node_t *)malloc(sizeof(Mesh_Triangle_Node_t));
  node->t = t;
  node->next = L->head;
  node->prev = NULL;
  if (L->head == NULL) {
    /*New head*/
    L->head = L->tail = node;
  }
  else {
    L->head->prev = node;
    L->head = node;
  }
}

/*-----------------------------------------------------------------------*/
void Mesh_Triangle_List_union (Mesh_Triangle_List_t *L1, Mesh_Triangle_List_t *L2) {
  if ( (L1->head != NULL) && (L2->head != NULL) ) {
    L1->tail->next = L2->head;
    L2->head->prev = L1->tail;
    L1->tail = L2->tail;;
  }
  else if (L2->head != NULL) {
    L1->head = L2->head;
    L1->tail = L2->tail;
  }
}

/*-----------------------------------------------------------------------*/
Mesh_Triangle_Node_t* Mesh_Triangle_List_remove (Mesh_Triangle_List_t *L, Mesh_Triangle_Node_t *node) {
  if ((node->prev == NULL) && (node->next == NULL)) {
    free (node);
    L->head = NULL;
    L->tail = NULL;
    return NULL;
  }
  else if (node->prev == NULL) {
    node->next->prev = NULL;
    L->head = node->next;
    free (node);
    return L->head;
  }
  else if (node->next == NULL) {
    node->prev->next = NULL;
    L->tail = node->prev;
    free (node);
    return NULL;
  }
  else {
    Mesh_Triangle_Node_t *next = node->next;
    node->next->prev = node->prev;
    node->prev->next = next;
    free (node);
    return next;
  }
}

v3 R3_Mesh_Side_slice (v3 vi, v3 vj, float Z) {
   double dx = vj.x - vi.x;
   double dy = vj.y - vi.y;
   double dz = vj.z - vi.z;
   assert(dz != 0);
   double frac = (Z - vi.z)/dz;
   float xint = (float)(frac*dx + (double)vi.x);
   float yint = (float)(frac*dy + (double)vi.y);
   return (v3){ xint, yint, Z };
}


typedef struct _contour {
   bool external;
   bool clockwise;
   vector<v3> P;
} contour;

//Rounds x to an integer multiple of eps
float xround (float x, double eps, int mod, int rem) {
  double y = round((double)x/(mod * eps));
  double z = (y * mod + rem) * eps;
  return (float)z;
}

// rounding a vertex (x,y,z) to an even multiple of eps
// v3 v3_round(float x,float y,float z,double eps){
//   v3 p;
//   p.x = xround(x, eps, 2, 0);
//   p.y = xround(y, eps, 2, 0);
//   p.z = xround(z, eps, 2, 0);
//   return p;
// }
v3 v3_round(v3 point,double eps){
  v3 p;
  p.x = xround(point.x, eps, 2, 0);
  p.y = xround(point.y, eps, 2, 0);
  p.z = xround(point.z, eps, 2, 0);
  return p;
}

//returns a triangle given 3 vertices and normal after rounding off easch vertice to even multiple of eps
// Triangle make_triangle(float v12[FACET_READ_SIZE],double eps){
//     return Triangle (v3(v12[0],v12[1],v12[2]), v3_round(v12[3],v12[4],v12[5], eps), v3_round(v12[6],v12[7],v12[8], eps), v3_round(v12[9],v12[10],v12[11], eps));
// }
Triangle make_triangle(v3 normal, v3 vertex1, v3 vertex2, v3 vertex3 ,double eps){
    return Triangle (normal, v3_round(vertex1, eps), v3_round(vertex2, eps), v3_round(vertex3, eps));
}

// checks if any vertices of a triangle are coinciding
bool degenerate (Triangle t) {
  if (t.v[0].distTo(t.v[1]) < 0.000001) { return true; }
  if (t.v[1].distTo(t.v[2]) < 0.000001) { return true; }
  if (t.v[2].distTo(t.v[0]) < 0.000001) { return true; }
  return false;
}

LineSegment R3_Mesh_Triangle_slice (Mesh_Triangle_Node_t *t, float Z) {
   assert((t->t.zMin < Z) && (t->t.zMax > Z));
   int np = 0; /* Number of segment endpoints found */
   LineSegment seg;
   for (int i = 0; i < 3; i++) {
      /* Get side {i} of triangle: */
      int j = (i == 2 ? 0 : i+1);
      v3 vi = (t->t.v[i]);
      v3 vj = (t->t.v[j]);
      /* Check for intersection of plane with {vi--vj}. */
      /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
      float vzMin = (vi.z < vj.z ? vi.z : vj.z);
      float vzMax = (vi.z > vj.z ? vi.z : vj.z);
      if ((vzMin <= Z) && (vzMax > Z)) {
         v3 p = R3_Mesh_Side_slice (vi, vj, Z);
         assert(np < 2);
         seg.v[np] = p;
         np++;
      }
   }
   //cout<<np;
   assert(np == 2);
   return seg; 
}

/*-----------------------------------------------------------------------*/
// /*Gets an arbitrary segment from {H}, removes it from {H} and returns it as a trivial chain. */
vector<v3> IncrementalStartLoop(vector<PointMesh> &H) {
   vector<v3> P;
   auto it = (H[0]).begin();
   v3 u = (*it).first;
   vector<v3> vw = (*it).second;
   v3 v = vw.at(0);
   P.push_back(u);
   P.push_back(v);
   auto temp = std::remove((H[0][u]).begin(), (H[0][u]).end(), v);
   (H[0][u]).erase(temp, (H[0][u]).end());
   if (H[0][u].size() == 0) { (H[0]).erase(u); }
   (H[0][v]).erase(std::remove((H[0][v]).begin(), (H[0][v]).end(), u), (H[0][v]).end());
   if (H[0][v].size() == 0) { (H[0]).erase(v); }
   return P;
}

/*-----------------------------------------------------------------------*/
/*Extends the chain {P} wih segments from {H}, removing them, while possible. */
void IncrementalExtendLoop(vector<v3> &P, vector<PointMesh> &H) { 
  int index = 0;
  int n = P.size();
  v3 first = P.front();  
  v3 current = P.back(); 
  v3 last;
         
  /* Collect other vertices: */
  while (true) {
    auto it = (H[0]).find(current);
    if (it == (H[0]).end()) { /*Vertex {current} is a dead end:*/ break; }
    v3 key1 = (*it).first; assert(key1 == current);  /*Paranoia check.*/
            
    /*Get {next}, the first unused neighbor of {current}:*/
    vector<v3> vw = (*it).second; /*Unused neighbors of {current}.*/
    assert (vw.size() != 0); 
    v3 next = vw.at(0); /*First unused neighbor of {current}.*/

    /*Append the segment {(current,next)} to {P} and delete from {H}:*/
    P.push_back(next);

    /*Remove the segment {(current,next)} from {H}:*/
    (H[0][current]).erase(std::remove((H[0][current]).begin(), (H[0][current]).end(), next), (H[0][current]).end());
    if (H[0][current].size() == 0) { (H[0]).erase(current); } 
    (H[0][next]).erase(std::remove((H[0][next]).begin(), (H[0][next]).end(), current), (H[0][next]).end());
    if (H[0][next].size() == 0) { (H[0]).erase(next); } 

    if (next == first) { /*We closed a loop:*/ break; }

    /*Move on:*/
    current = next;
  }
}

/*Reverses the chain {P}.*/
void IncrementalReverseLoop(vector<v3> &P) { 
  std::reverse(P.begin(),P.end());
}

void ContourConstruction (vector<LineSegment> segs, vector<contour> polygons[], int plane) {

   bool verbose = false;

   //clock_t contour_begin = clock();

   /*Creating the hash table.*/
   vector<PointMesh> H(1);

   /*Rounding vertices and filling the hash table.*/
   double eps = 1/128.0;
   for (std::vector<LineSegment>::iterator i = segs.begin(); i != segs.end(); i++) {
      LineSegment q = *i;
      q.v[0].x = round(q.v[0].x / eps) * eps;
      q.v[0].y = round(q.v[0].y / eps) * eps;
      q.v[0].z = plane;
      q.v[1].x = round(q.v[1].x / eps) * eps;
      q.v[1].y = round(q.v[1].y / eps) * eps;
      q.v[1].z = plane;
      if (q.v[0].distTo(q.v[1]) > 0.0001) {
         (H[0][q.v[0]]).push_back(q.v[1]);
         (H[0][q.v[1]]).push_back(q.v[0]);
      }
   }

   /* Count vertices by degree: */
   if (verbose) {
     int degmax = 10;
     int ctdeg[degmax+1];
     for (int deg = 0; deg <= degmax; deg++) { ctdeg[deg] = 0; }
     for (auto i = (H[0]).begin(); i != (H[0]).end(); i++) {
        vector<v3> L = (*i).second;
        int deg = L.size();
        if (deg > degmax) { deg = degmax; }
        ctdeg[deg]++;
     }
     assert(ctdeg[0] == 0);
     bool closedSlice = true;
     for (int deg = 1; deg <= degmax; deg++) { 
       if (((deg % 2) != 0) && (ctdeg[deg] > 0)) { closedSlice = false; }
       if ((verbose || (deg != 2)) && (ctdeg[deg] != 0))
         { cout << "there are " << ctdeg[deg] << " vertices of degree " << deg << " on plane " << plane << endl; }
     }
     if (!closedSlice) { cout << "** contours of plane " << plane << " are not closed" << endl; }
   }

   /*Contour construction.*/
   bool maximal = true;
   while (!(H[0]).empty()) {
     if (maximal) {
       vector<v3> P = IncrementalStartLoop(H);
       IncrementalExtendLoop(P,H);
       if (P.front() != P.back()) { //Chain {P} is open
         IncrementalReverseLoop(P);
         IncrementalExtendLoop(P,H);
       }
       polygons[plane].push_back({false, false, P});
     }
     else {
       vector<v3> P = IncrementalStartLoop(H);
       IncrementalExtendLoop(P,H);
       polygons[plane].push_back({false, false, P});
     }
   }
   //clock_t contour_end = clock();
   //loopclosure_time += double(contour_end - contour_begin)/CLOCKS_PER_SEC;
}

Mesh_Triangle_List_t** IncrementalSlicing_buildLists (bool srt, double delta, const TriangleMesh *mesh, vector<float> P) {

  int k = P.size(); /* Number of planes. */

  Mesh_Triangle_List_t **L = (Mesh_Triangle_List_t **)malloc((k+1) * sizeof(Mesh_Triangle_List_t *));

  for (size_t p = 0; p <= k; p++) { L[p] = Mesh_Triangle_List_create(); }

  const vector<Triangle> &T = mesh->getvTriangle();

  int n = T.size(); /* Number of triangles. */

    /* Uniform slicing - compute list index: */
    for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
      Triangle t = *it;
      int p;
      if (t.zMin < P[0]) {
        p = 0;
      }
      else if (t.zMin > P[k-1]) {
        p = k;
      }
      else {
        p = floor((t.zMin - P[0])/delta) + 1;
      }
      Mesh_Triangle_List_insert (t, L[p]);
    }

  //MY EDIT: check if this is correct;
  return L;
}
  
void IncrementalSlicing (const TriangleMesh *mesh, vector<float> P, float delta, bool srt, vector<contour> polygons[], bool chaining, bool orienting) {

  /*Slicing*/

  int k = P.size();

  vector<LineSegment> segs[k];

  /* Classify triangles by the plane gaps that contain their {zMin}: */
  Mesh_Triangle_List_t **L = IncrementalSlicing_buildLists (srt, delta, mesh, P);
  /* Now perform a plane sweep from bottom to top: */

  Mesh_Triangle_List_t *A = Mesh_Triangle_List_create(); /* Active triangle list. */
  for (int p = 0; p < k; p++) {
    /* Add triangles that start between {P[p-1]} and {P[p]}: */
    Mesh_Triangle_List_union (A, L[p]);
    /* Scan the active triangles: */
    Mesh_Triangle_Node_t *aux = A->head;
    while (aux != NULL) {
      Mesh_Triangle_Node_t *next = aux->next;
      if (aux->t.zMax < P[p]) {
        /* Triangle is exhausted: */
        Mesh_Triangle_List_remove (A, aux);
      } else {
        /* Compute intersection: */
        if ((aux->t.zMin < P[p]) && (aux->t.zMax > P[p])) {
          LineSegment seg = R3_Mesh_Triangle_slice (aux, P[p]);
          segs[p].push_back(seg);
          //intersections++;
        }
      }
      aux = next;
    }
  }
  free(L);
  /*End-Slicing*/

    if (chaining) {
    /*Contour construction:*/
    for (size_t p = 0; p < k; p++) {
      if (!segs[p].empty()) {
          ContourConstruction (segs[p], polygons, p);
         #ifdef DEBUG
           char fname[256];
           sprintf(fname, "slice_%03d.txt", (int)p);
           FILE *fslice = fopen(fname, "w");
           fprintf(fslice, "----------- Segmentos ---------------\n");
           for (int ii = 0; ii < segs[p].size(); ii++) {
              LineSegment ss = segs[p].at(ii); 
              fprintf(fslice, "%f %f   %f %f\n", ss.v[0].x, ss.v[0].y, ss.v[1].x, ss.v[1].y);           
           }
           fprintf(fslice, "---------- End Segmentos --------------\n");
           record_polygons (polygons[p], fslice);
           fclose(fslice);
         #endif
        //  if (orienting) {   
        //    ray_casting (polygons[p]);
        //  }
         segs[p].clear(); 
      }
    }
    /*End construction.*/
  }
  else {
    // export_svg_no_chaining ("segments.svg", segs, k, mesh->meshAABBSize());
  }
}

/*Compute uniform and adaptive z-plane coordinates!*/
vector<float> compute_planes (const TriangleMesh *mesh, float max_thickness, double eps, float *delta) {

  bool rounding = true; /*To avoid that the Z-coordinates of all planes are distinct from the Z-coordinates of all vertices.*/

  /* Vector to keep the plane coordinates: */
  vector<float> Planes;

  /* Assuming the model as a 3D axis-aligned bounding-box: */
  double model_zmax = std::max(mesh->getUpperRightVertex().z, mesh->meshAABBSize().z);

  double model_zmin = mesh->getBottomLeftVertex().z;

    double spacing = (rounding ? xround (max_thickness, eps, 2, 0) : max_thickness); /*Plane spacing even multiple of {eps}*/

    double P0 = xround (model_zmin - spacing, eps, 2, 1); /*First plane odd multiple of {eps}.*/

    int no_planes = 1 + (int)((model_zmax + spacing - P0)/spacing); /* Number of planes: */

    cout << "eps = " << eps << endl;
    cout << "max thickness = " << max_thickness << endl;
    cout << "rounded plane spacing spacing = " << spacing << endl;
    cout << "model zmin = " << model_zmin << ", model zmax = " << model_zmax << ", first plane Z = " << P0 << ", number of planes = " << no_planes << endl;

    for (size_t i = 0; i < no_planes; i++) {
      /* Building the vector with the slice z coordinates: */
      float Pi = (float)(P0 + i * spacing);
      if ((Pi > model_zmin) && (Pi < model_zmax)) {
          Planes.push_back ((float)(P0 + i * spacing));
      }
    }
    *delta = (float)(spacing);

    //MY EDIT: check if this is correct;
    return Planes;
  }

  
float parse_float(std::ifstream& s) {
  char f_buf[sizeof(float)];
  s.read(f_buf, 4);
  float* fptr = (float*) f_buf;
  return *fptr;
}

v3 parse_point(std::ifstream& s) {
  float x = parse_float(s);
  float y = parse_float(s);
  float z = parse_float(s);
  return v3(x, y, z);
}

//function to read triangles from STL file and store it in the TriangleMesh class
void stlToMeshInMemory (const char *stlFile, TriangleMesh *mesh, double eps){

  // int no_of_triangles = 0;                        //to store no of triangles
  char no_of_triangles[4];                        //to store no of triangles
  char var[MAX_READ_CHUNK];                       //to read and store header from input file
  float v12[FACET_READ_SIZE];                    //to store the vector read from input file
  char attribute[2];                                 //to store the attribute
  int ndegenerated = 0;

  cout<<"[MY PERSONAL DEBUGGER]" << "index 000 : " << "Entry point Hit!!!" << endl;
  ifstream stlRead(stlFile, ios::binary);
  cout<<"[MY PERSONAL DEBUGGER]" << "index 000 : " << "Reader initialized..." << endl;
  
  stlRead.read(var,MAX_READ_CHUNK);
  cout<<"[MY PERSONAL DEBUGGER]" << "index 0 : " << var << endl;
  stlRead.read(no_of_triangles,4);
  unsigned int* r = (unsigned int*) no_of_triangles;
  unsigned int num_triangles = *r;
  cout<<"[MY PERSONAL DEBUGGER]" << "index 1 : " << num_triangles << endl;

  for(int i = 0;i < num_triangles;i++){
    //stlRead.read((char*)v12,sizeof(float)*FACET_READ_SIZE);
    v3 normal = parse_point(stlRead);
    v3 point1 = parse_point(stlRead);
    v3 point2 = parse_point(stlRead);
    v3 point3 = parse_point(stlRead);
    // cout<<"[MY PERSONAL DEBUGGER]" << "read partially successfully : " << normal.x << " | " << normal.y << " | " << normal.z << endl;
    // cout<<"[MY PERSONAL DEBUGGER]" << "read partially successfully : " << point1.x << " | " << point1.y << " | " << point1.z << endl;
    // cout<<"[MY PERSONAL DEBUGGER]" << "read partially successfully : " << point2.x << " | " << point2.y << " | " << point2.z << endl;
    // cout<<"[MY PERSONAL DEBUGGER]" << "read partially successfully : " << point3.x << " | " << point3.y << " | " << point3.z << endl;
    stlRead.read(attribute,2);
    // cout<<"[MY PERSONAL DEBUGGER]" << "read successfull : " << i << endl;

    Triangle t = make_triangle(normal, point1, point2, point3 ,eps);

    if (!degenerate(t)) {
       mesh->push_back (t);
    }
    else {
      ndegenerated++;
    }
  }

  stlRead.close();
  cout << "number of degenerated triangles = " << ndegenerated << endl;
  cout << "[MY PERSONAL DEBUGGER]" << "SIZE " << mesh->size() << endl;
  cout << "[MY PERSONAL DEBUGGER]" << "TOTAL " << (mesh->size()+ndegenerated) << endl;
}

int main(int argc, char** argv)
{
  bool chaining = true;

  double eps = 0.004;

  char *model;

  cout<<"[MY PERSONAL DEBUGGER]" << "index 000 : " << "Reading arguments" << endl;
  cout << argc << " | ";
  for (int i=0; i<argc; i++)
    cout << argv[i] << " | ";
  cout << endl << endl;

    if (strcmp(argv[2], "-model") == 0) { // takes model label as second  argument to command line
     model = argv[3];                     // takes path to stl as third argument to command line
    }

  float max_thickness, delta;

  if (strcmp(argv[4], "-thickness") == 0) { // takes thickness label as second  argument to command line
    max_thickness = atof(argv[5]);          // takes thickness params as third argument to command line
  }  
  else {
    printf("Error: specify the slicing spacing in mm (thickness)!!!\n");
  }

  //   char *adaptive;

  // if (strcmp(argv[6], "-adaptive") == 0) {
  //   adaptive = argv[7];
  // }  

  char *write_option = argv[8];

  // char *rotate;

  // if (strcmp(argv[9], "-rotate") == 0) {
  //    rotate = argv[10];
  // }

  bool orienting;

  if (strcmp(argv[11], "-orienting_contours") == 0) {
    if (strcmp(argv[12],"true") == 0) {
      orienting = true;
    }
    else {
      orienting = false;
    }
  }

  TriangleMesh mesh;

  cout<<"[MY PERSONAL DEBUGGER]" << "index 000 : " << "Entering the universe..........." << endl;

  stlToMeshInMemory (model, &mesh, eps);

  vector<float> P = compute_planes (&mesh, max_thickness, eps, &delta);

  int nplanes = P.size();

  cout<<"[MY PERSONAL DEBUGGER]" << "Number of planes : " << nplanes << endl;

  vector<contour> polygons[nplanes];
 
  bool srt = false;

    if (strcmp(write_option, "-No") == 0) {
     chaining = false;
  }

  IncrementalSlicing (&mesh, P, delta, srt, polygons, chaining, orienting);

  return 0;
}
