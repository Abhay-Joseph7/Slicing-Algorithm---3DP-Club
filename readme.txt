Execution Commands(on MinGW compiler):

STEP I : in terminal navigate to SLICER_3DP_PROJECT-MASTER

STEP II : For compiling into .exe file:
    g++ main.cpp v3.cpp TriangleClass.cpp TriangleMesh.cpp LineSegment.cpp Mesh_Triangle_List.cpp -o final_executable.exe

STEP III : For running compiled versions:
    ./final_executable.exe testing -model "stl_models/01.liver.stl" -thickness 0.02 -inni -mimmi -Yes chini mo -orienting_contours true chimapanzee
    
    Tryout different models:
     * ./final_executable.exe testing -model "stl_models/02.femur.stl" -thickness 0.02 -inni -mimmi -Yes chini mo -orienting_contours true chimapanzee
     * ./final_executable.exe testing -model "stl_models/03.bunny.stl" -thickness 0.02 -inni -mimmi -Yes chini mo -orienting_contours true chimapanzee
     * ./final_executable.exe testing -model "stl_models/cube_continuous.stl" -thickness 0.02 -inni -mimmi -Yes chini mo -orienting_contours true chimapanzee

TODO : 
* start from line 396: keep on uncommenting functions until line 206, main.cpp
* What does "chaining" argument does? Is there any unidentified arguments? 
* learn to use "DEBUG" mode, "verbose"
* what is .first and .second in "vector.h" library ? 
* Assert fails when you run "stl_models/01.liver.stl" model. Find its root cause...