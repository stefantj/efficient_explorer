//
//  utilities.h
//  FMT
//
//  Created by Megamind on 7/23/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#ifndef FMT_utilities_h
#define FMT_utilities_h


#include <functional>
#include <chrono>
#include <math.h>
#include <string>

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

// Definitions for map states
#define MAP_OCCU        127    // This cell is occupied
#define MAP_GOCC_U      126    // Guess that this cell is occupied (actual state is unknown)
#define MAP_GOCC_F      125    // Guess that this cell is occupied (actual state is free) Used for collision avoidance.
#define MAP_OTHER       66     // Treat others as obstacles
#define MAP_OCC_THRESH  65     // Used to determine whether cell is occupied
#define MAP_UNKN        64     // This cell is unknown
#define MAP_GOAL_U      63     // Marks goal centroid
#define MAP_FREE        59     // This and anything less is free
#define MAP_GOAL_F      58
#define MAP_TRACE       6
#define MAP_SELF        5
#define MAP_CENT        2      // Marks cluster centroid
#define MAP_FRONT       1      // Marks frontier (always empty if frontier)


#define ANSI_COLOR_BLACK    "\x1b[30m"
#define ANSI_COLOR_RED      "\x1b[31m"
#define ANSI_COLOR_GREEN    "\x1b[32m"
#define ANSI_COLOR_YELLOW   "\x1b[33m"
#define ANSI_COLOR_BLUE     "\x1b[34m"
#define ANSI_COLOR_MAGENTA  "\x1b[35m"
#define ANSI_COLOR_CYAN     "\x1b[36m"
#define ANSI_COLOR_GRAY     "\x1b[37m"
#define ANSI_COLOR_RESET    "\x1b[0m"
#define ANSI_BG_CYAN        "\x1b[46m"
#define ANSI_BG_RED         "\x1b[41m"
#define ANSI_BG_GREEN       "\x1b[42m"
#define ANSI_BG_YELLOW      "\x1b[43m"
#define ANSI_BG_BLUE        "\x1b[44m"

// Used to encode a point in 3D space
struct Point {
    float x;
    float y;
    float z;
};

inline float dist(Point x1, Point x2){
    return sqrtf((x1.x-x2.x)*(x1.x-x2.x) + (x1.y-x2.y)*(x1.y-x2.y) + (x1.z-x2.z)*(x1.z-x2.z));
}

inline float norm(Point x){
    return sqrtf((x.x*x.x + x.y*x.y + x.z*x.z));
}

float nice_angle(float theta);

// Removes element from a sorted collection. Assumes that first element is smallest
int remove_element(int* collection, int collection_size, int element);

// Returns index of element. Returns -1 if not found.
int find_element(int* collection, int collection_size, int element);

// Joins two sorted sets and puts the union into set1
int join_sets(int* set1, int set1_size, int* set2, int set2_size);

// Used to represent obstacles
struct Map {
    char* data = nullptr;
    std::size_t X_dim = 1;
    std::size_t Y_dim = 1;
    std::size_t Z_dim = 1;
    
    void set(std::size_t idx, std::size_t idy, std::size_t idz, char val){
        if(idx < X_dim && idy < Y_dim && idz < Z_dim){
            data[((idx*Y_dim + idy)*Z_dim + idz)] = val;
        }
    }

    void set(Point loc, char val){
        if(int(loc.x) < X_dim && int(loc.y) < Y_dim && int(loc.z) < Z_dim){
            data[((int(loc.x)*Y_dim + int(loc.y))*Z_dim + int(loc.z))] = val;
        }
    }

    bool is_free(Point loc){
        if(int(loc.x) < X_dim && int(loc.y) < Y_dim && int(loc.z) < Z_dim){
            return(data[((int(loc.x)*Y_dim + int(loc.y))*Z_dim + int(loc.z))] < MAP_OCC_THRESH);
        }
        return false;
    }
    
    std::size_t ind2num(std::size_t idx, std::size_t idy, std::size_t idz){
        if(idx < X_dim && idy < Y_dim && idz < Z_dim){
            return ((idx*Y_dim + idy)*Z_dim + idz);
        }else{
            return 0;
        }
    }
    
    void num2ind(Point* pt, std::size_t num){
        pt->z = num % Z_dim;
        pt->y = float(int((num-pt->z)/Z_dim)%int(Y_dim));
        pt->x = float(int((num-pt->z -Z_dim*pt->y)/Y_dim));
    }
    
    char operator()(std::size_t idx, std::size_t idy, std::size_t idz){
        if(idx < X_dim && idy < Y_dim && idz < Z_dim){
            return (data[((idx*Y_dim + idy)*Z_dim + idz)]);
        }else{
            return 255;
        }
    };
};

// Helper for collision checker
bool check_map(const float* pt,  Map* map, bool* confirmed = NULL);

// Returns true if a collision found
bool collision(const Point pt1, const Point pt2,  Map* map, bool* confirmed = NULL);


struct Cluster {
    int size;
    int number;
    int* members;
    Point median;
};

// Place the contents of both clusters into the smalleter numbered one.
void merge_clusters(Cluster* C1, Cluster* C2);


void save_julia_var(FILE* f, std::string var_name, float* variable, int num_vars);
void save_julia_var(FILE* f, std::string var_name, Point* variable, int num_vars);
void save_julia_var(FILE* f, std::string var_name, Map* variable);



#endif
