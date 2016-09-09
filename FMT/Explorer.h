//
//  Explorer.h
//  FMT
//
//  Handles creation of frontiers, reduction to goals, and assignment.
//
//  Created by Megamind on 7/22/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#ifndef __FMT__Explorer__
#define __FMT__Explorer__

#include <iostream>
#include <vector>
#include "fmt.h"
#include "simulator.h"
#include "utilities.h"
#include "Hungarian.h"

#define MAX_CLUSTERS    256
#define MAX_CLUSTERSIZE 1048
#define MAX_TEAMSIZE    32
#define F_SPARSE_PTS    400
#define F_DENSE_PTS     1024

//#define EXPLORE_DEBUG
//#define EXPLORE_STATUS
//#define EXPLORE_TIMING

class Explorer{
public:
    Explorer(int agent_ID, int map_x, int map_y);
//    void get_updated_map();
    ~Explorer();
    
    // Measurement functions
    void update_state(Point* positions, int num_agents);
    void update_map(Map* measurement, Point lower, Point upper);
    void merge_map(Map* other_map);
    Map* get_map(){return &exploration_map;}
    void force_map_update(Map* map);
    void clear_map();
    
    // Planning functions
    bool cluster_frontiers();
    int  compute_costs();
    void compute_cost(Point start, Point goal, Point* path, int goal_num);
    void compute_cost(Point start, Cluster* goal, Point* path, int goal_num);
    void assign(Point* waypoint);
    
    // Plotting utilities
    void plot_clusters();
    void print_costs();
    void print_call_rates();
    
    
    // Service functions
    int get_clusters(Point* cluster_centroids, float* costs, int max_num_clusters);
    void get_neighbor_clusters(int agent_id, Point* cluster_centroids, float* costs, int n_clusters);
    
    Point current_plan[F_DENSE_PTS+2];    // Holds the current plan
private:
    
    void relabel_frontiers(int clust1_ind, int clust2_ind);
    
    void print_cluster(Cluster* clust);
    
    
    // Tracks the size of the team (as best known to the agent)
    int self_id;
    int team_size;
    Point swarm_locs[MAX_TEAMSIZE];
    
    // Variables for clustering frontiers
    int* frontiers;    // Index of frontiers
    int num_frontiers; // Number of frontiers
    int* cluster_ids;  // List of ids for each frontier

    float merge_dist;   // Distance to cluster over
    Cluster clusters[MAX_CLUSTERS];  // List of active clusters
    int num_clusters;   // Number of active clusters
    int active_clusters[MAX_CLUSTERS];    // Index of non-squashed clusters
    int num_active_clusters; // Number of active clusters

    float goal_costs[MAX_TEAMSIZE][MAX_CLUSTERS];
    Point goal_locations[MAX_CLUSTERS]; // Holds the goal locations of the swarm
    int num_goals;                      // Holds the number of goals for the swarm
    
    Point paths[MAX_CLUSTERS][F_DENSE_PTS+2];     // Holds the paths to goals

    
    //Map for exploration
    Map exploration_map;
    
    // FMT planners for computing distances/trajectories
    FMT* F_sparse;
    FMT* F_dense;
    int f_sparse_calls;
    int f_dense_calls;
    
};


#endif /* defined(__FMT__Explorer__) */
