//
//  simulator.h
//  FMT
//
//  Handles simple physics and message passing
//
//  Created by Stefan Jorgensen
//  MIT licence
//

#ifndef __FMT__simulator__
#define __FMT__simulator__

#include <iostream>
#include "utilities.h"
#include <stdlib.h>
#include <string>
class Explorer;

//#define SIM_DEBUG
//#define SIM_PLOT

#define MAX_NUM_AGENTS  32
#define MAX_NUM_GOALS  256

struct Agent{
    Point location;
    float bearing;
    Point target;
    Explorer* planner;
};


class Simulator{
public:
    Simulator();
    void run_simulator(int iters, int team_size);
    
private:
    void raytrace(int agent);
    void move_base(int agent, int step);
    
    void init_video_file(FILE* f);
    void print_video_file(FILE* f, std::string filename, float time);
    
    Map true_map;
    Map global_map;

    Point agent_location_buffer[MAX_NUM_AGENTS][MAX_NUM_AGENTS]; // Holds buffer for location data of agents

    float agent_costs_buffer[MAX_NUM_AGENTS][MAX_NUM_GOALS];
    Point agent_clusters_buffer[MAX_NUM_AGENTS][MAX_NUM_GOALS];
    int agent_num_clusters[MAX_NUM_AGENTS];
    Agent agents[MAX_NUM_AGENTS]; // List of agents
    
    float delta_time;
    int plot_delay;
    int replan_delay;
    int num_agents;
    
    bool comms[MAX_NUM_AGENTS][MAX_NUM_AGENTS];
};

#endif /* defined(__FMT__simulator__) */
