//
//  fmt.h
//  FMT
//  Fed up with vectors, going with arrays. Forgive my C-ness.
//  Created by Megamind on 7/18/16.
//  MIT licence

#ifndef FMT_fmt_h
#define FMT_fmt_h

//Used to toggle debug text
//#define FMT_DEBUG
//#define FMT_TIMING
#define FMT_JULIA_DEBUG

// Connection types
#define RAD_CON 1
#define KNN_CON 2

#define CACHE_UNK    0  // Untested connection
#define CACHE_FREE   1  // Tested connection, guaranteed free path
#define CACHE_FREE_U 2  // Tested connection, free path until state update
#define CACHE_OCC    3  // Tested connection, guaranteed occupied path
#define CACHE_OCC_U  4  // Tested connection, occupied path until state update

#include <math.h>
#include <iostream>
#include <functional>
#include <queue>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <time.h>
#include "smoother.h"

#include "utilities.h"

class PathSmoother;

// Used for describing neighborhoods of a point
typedef struct Neighborhoods {
    uint8_t* cache;   // Whether or not the collision has been checked
    int* indices; // Indices of neighbor points
    float* costs; // Costs of traveling to neighbor points
    int size;     // Number of neighbors
} Neighborhood;

// Used to set behavior of FMT*
typedef struct FMT_Parameters {
    int X_limit;            // Size of box in x dimension
    int Y_limit;            // Size of box in y dimension
    int Z_limit;            // Size of box in z dimension (set to 0 for planning in 2D)
    int num_pts;            // Number of points to sample
    int connection_type;    // 1 = r-connected, 2 = k-nearest neighbors
    float connection_param; // either radius or number of neighbors
} FMT_Params;


class FMT {
public:
    /*** Constructor methods  ***/
    // Default Constructor
    FMT(){
        // Call explicit constructor with default parameters
        printf("Error: default constructor doesn't work yet\n");
    }
    
    // Destructor
    ~FMT();
    
    // Explicit constructor
    FMT(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param);
    
    // Construct new FMT instance from saved parameter file
    FMT(std::string filename){
        is_initialized = false;
        printf("Creating FMT instance from file is not supported yet\n");
    }
    
    /*** Planning methods ***/
    //Plan with a fresh instance
    void fmtstar(Point start, Point goal, Map* map, Point* path);
    
    //Plan a path between clusters
    void fmtstar(Cluster* start, Cluster* goal, Map* map, Point* path);
    
    //Plan a path from point to cluster
    void fmtstar(Point start, Cluster* goal, Map* map, Point* path);
    
    void print_datafile(Point* path, int path_length, Map* map);

    
    void clear_cache();
    
    void print_stats();
    
    PathSmoother Smoothed;
    
private:
    
    bool collision_inds(int ind1, int ind2, Map* map);
    
    void plan_path(Map* map, Point* path);
    
    // Sets the parameters of the planner
    void set_parameters(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param);
    
    // Samples new points and stores in points container
    void sample_points();
    
    // Compute neighborhoods for sampled points
    bool compute_neighborhoods();
    
    // Compute single neighborhood for given point at index
    bool compute_neighborhood(Point pt, int index, Map* map);

    // Compute single neighborhood for cluster of points point at index
    bool compute_neighborhood(Cluster* pt, int index, Map* map);
    
    // Undoes work of compute_neighborhood(Point pt, int index)
    // resets the state
    bool reset_neighborhood(int index);
    
    // Return true if points[index] is within goal region
    bool is_goal_pt(int index);
    
    
    void filter(Neighborhood* filtered_neighborhood, int index, int* filter_vec,int filter_size, int exclude);

    // Printing utilities
    void print_parameters();
    void print_neighborhood(Neighborhood n);
    void print_point(Point pt);
    void print_array(int* vec, int size);
    
    // Variables:
    FMT_Params parameters;                      // Parameters for FMT
    int max_neighborhood_size;
    int start_pt_ind;                           // Helper for startpt index
    int goal_pt_ind;                            // Helper for endpt index
    Point* points;                              // Sampled points
    Neighborhood* neighborhoods;    // Neighborhoods of sampled points
    bool is_initialized;                        // Must be true for computation to happen
    
    // Containers used for planning. Put them here so we only allocate memory once
    int* Windex;
    int* Active;
    int* Parents;
    float* Costs;
    int* H_new;
    Neighborhood x_near;
    Neighborhood y_near;
    
    // Tracking
    int num_skipped_checks;
    int num_total_checks;
};


#endif
