//
//  fmt.h
//  Created by Stefan Jorgensen
//  MIT licence

#ifndef FMT_fmt_h
#define FMT_fmt_h

//Used to toggle debug text
//#define FMT_DEBUG
#define FMT_TIMING
//#define FMT_JULIA_DEBUG

// Connection types
#define RAD_CON 1  // Radial
#define KNN_CON 2  // KNN (not implemented)

#define FMT_ORD 3  // Sampling order. 1 => positions only, 2 => pos + vel, 3 => pos,vel,acc etc.
#define POLY_ORD 3 // Order for polynomial fitting.

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

#include "polynomial.h"
#include "utilities.h"

class PathSmoother;

// Used to represent kino points for FMT
typedef struct FMT_Points {
    Point state[FMT_ORD];
} FMT_Point;


// Used for describing neighborhoods of a point
typedef struct Neighborhoods {
    uint8_t* cache;   // Whether or not the collision has been checked
    int* indices; // Indices of neighbor points
    PolyState* paths; // Polynomial paths which connect neighbors, in state space
    float* costs; // Costs of traveling to neighbor points
    int size;     // Number of neighbors
} Neighborhood;

// Used to set behavior of FMT*
typedef struct FMT_Parameters {
    int X_limit;            // Size of box in x dimension
    int Y_limit;            // Size of box in y dimension
    int Z_limit;            // Size of box in z dimension (set to 0 for planning in 2D)
    float X_vel_limit;        // Max vel in x dimension
    float Y_vel_limit;        // Max vel in y dimension
    float Z_vel_limit;        // Max vel in z dimension
    float X_acc_limit;        // Max acc in x dimension
    float Y_acc_limit;        // Max acc in y dimension
    float Z_acc_limit;        // Max acc in z dimension
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
    FMT(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param, Map* map_structure);
    
    /*** Planning methods ***/
    //Plan with a fresh instance
    int fmtstar(Point start, Point goal, Map* map, PolyState* path);
    
    // Clears cached collision checks
    void clear_cache();

    // Sets the initial state used for planning
    void set_initial_state(Point* vel, Point* acc);
    
    /*** Debug interfaces ***/
    void print_datafile(Point* path, int path_length, Map* map);

    void print_stats();
    
private:
    // Checks for collision using cache to speed up.
    bool collision_inds(int ind1, int ind2, Map* map);
    
    // Run the main loop of FMT, uses internal tree variables
    int plan_path(Map* map, PolyState* path);
    
    // Sets the parameters of the planner
    void set_parameters(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param);
    
    // Samples new points and stores in points container
    void sample_points();
    
    // Compute single neighborhood for given point at index
    bool compute_neighborhood(Point pt, int index, Map* map);
    
    // Undoes work of compute_neighborhood(Point pt, int index)
    // resets the state
    bool reset_neighborhood(int index);

    // Compute neighborhoods for sampled points
    bool compute_neighborhoods(Map* map_structure);
    
    // Return true if points[index] is within goal region
    bool is_goal_pt(int index);

    
    // Returns elements in a neighborhood which are in another list.
    void filter(Neighborhood* filtered_neighborhood, int index, int* filter_vec,int filter_size, int exclude);

    // Returns cached polynomial connecting start and goal
    PolyState* get_neighborhood_poly(int start_ind, int end_ind);
    
    // Printing utilities
    void print_parameters();
    void print_neighborhood(Neighborhood n);
    void print_point(Point pt);
    void print_array(int* vec, int size);
    void plot_neighborhood(int i);
    
    // Variables:
    bool is_initialized;                        // Must be true for computation to happen
    FMT_Params parameters;                      // Parameters for FMT

    int start_pt_ind;                           // Helper for startpt index
    int goal_pt_ind;                            // Helper for endpt index
    FMT_Point* points;                          // Sampled points
    int max_neighborhood_size;
    Neighborhood* neighborhoods;                // Neighborhoods of sampled points
    
    // Containers used for planning. Put them here so we only allocate memory once
    int* Windex;
    int* Active;
    int* Parents;
    float* Costs;
    int* H_new;
    Neighborhood x_near;
    Neighborhood y_near;
    PolynomialSmoother* smoother;
    
    Point init_vel;
    Point init_acc;
    
    // Tracking
    int num_skipped_checks;
    int num_total_checks;
};


#endif
