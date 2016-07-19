//
//  fmt.h
//  FMT
//
//  Created by Megamind on 7/18/16.
//  MIT licence

#ifndef FMT_fmt_h
#define FMT_fmt_h

//Used to toggle debug text
#define FMT_DEBUG
#define FMT_TIMING

// Connection types
#define RAD_CON 1
#define KNN_CON 2


#include <functional>
#include <queue>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <time.h>

// Used to encode a point in 3D space
typedef struct Points {
    float x;
    float y;
    float z;
} Point;

// Used for describing neighborhoods of a point
typedef struct Neighborhoods {
    std::vector<int> indices; //Indices of neighbor points
    std::vector<float> costs; //Costs of traveling to neighbor points
    int size;                 //Number of neighbors
} Neighborhood;

// Used to set behavior of FMT*
typedef struct FMT_Parameters {
    int X_limit;            // Size of box in x dimension
    int Y_limit;            // Size of box in y dimension
    int Z_limit;            // Size of box in z dimension (set to 0 for planning in 2D)
    int num_pts;            // Number of points to sample
    int connection_type;    // 1 = r-connected, 2 = k-nearest neighbors
    float connection_param; // either radius or number of neighbors
    float goal_radius;
} FMT_Params;


class FMT {
public:
    /*** Constructor methods  ***/
    // Default Constructor
    FMT(){
        // Call explicit constructor with default parameters
        printf("Error: default constructor doesn't work yet\n");
    }

    // Explicit constructor
    FMT(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param, float goal_radius);
    
    // Construct new FMT instance from saved parameter file
    FMT(std::string filename){
        is_initialized = false;
        std::printf("Creating FMT instance from file is not supported yet\n");
    }
    
    /*** Planning methods ***/
    //Plan with a fresh instance
    std::vector<int> fmtstar(Point start, Point goal);
    
private:

    // Sets the parameters of the planner
    void set_parameters(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param, float goal_radius);
    
    // Samples new points and stores in points container
    void sample_points();

    // Compute neighborhoods for sampled points
    bool compute_neighborhoods();
    
    // Compute single neighborhood for given point at index
    bool compute_neighborhood(Point pt, int index);
    
    // Undoes work of compute_neighborhood(Point pt, int index)
    // resets the state
    bool reset_neighborhood(int index);
    
    // Return true if points[index] is within goal region
    bool is_goal_pt(int index);
    
    // Return true if collision between points
    bool collision(Point pt1, Point pt2);
    
    // utility functions:
    // Returns the distance between two points
    inline float dist(Point x1, Point x2){
        return sqrtf((x1.x-x2.x)*(x1.x-x2.x) + (x1.y-x2.y)*(x1.y-x2.y) + (x1.z-x2.z)*(x1.z-x2.z));
    }
    
    void remove_element(std::vector<int> collection, int element);
    Neighborhood filter(int index, std::vector<int> filter_vec, int exclude);
    void remove_neighbor(int index, int neighbor_index);
    void print_neighborhood(Neighborhood n);
    void print_point(Point pt);
    void print_vec(std::vector<int> vec);
    
    // Variables:
    FMT_Params parameters;                      // Parameters for FMT
    int start_pt_ind;                           // Helper for startpt index
    int goal_pt_ind;                            // Helper for endpt index
    std::vector<Point> points;                  // Sampled points
    std::vector<Neighborhood> neighborhoods;    // Neighborhoods of sampled points
    bool is_initialized;                        // Must be true for computation to happen
};


#endif
