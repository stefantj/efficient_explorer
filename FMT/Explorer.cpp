//
//  Explorer.cpp
//  FMT
//
//  Created by Stefan Jorgensen
//  MIT licence
//


#include "Explorer.h"

// Constructor, creates a new agent
Explorer::Explorer(int agentID, int map_x, int map_y){
    self_id = agentID;
    team_size = 6;
    merge_dist = 15;
    
    // Initialize swarm_locs to shut valgrind up:
    for(int i =0; i < MAX_TEAMSIZE; i++){
        swarm_locs[i].x = 0.0;
        swarm_locs[i].y = 0.0;
        swarm_locs[i].z = 0.0;
    }

    target.x = swarm_locs[self_id].x;
    target.y = swarm_locs[self_id].y;
    target.z = swarm_locs[self_id].z;
    
    frontiers = new int[map_x*map_y];
    num_frontiers = 0;
    cluster_ids = new int[map_x*map_y];
    num_clusters = 0;
    num_active_clusters = 0;
    num_goals = 0;
    
    for(int i = 0; i < MAX_CLUSTERS; i++){
        active_clusters[i] = 0;
        clusters[i].members = new int[MAX_CLUSTERSIZE];
        clusters[i].size = 0;
        goal_locations[MAX_CLUSTERS].x = 0.0;
        goal_locations[MAX_CLUSTERS].y = 0.0;
        goal_locations[MAX_CLUSTERS].z = 0.0;
    }
    
    for(int i = 0; i < MAX_TEAMSIZE; i++){
        for(int j = 0; j <MAX_CLUSTERS; j++){
            goal_costs[i][j] = 0;
        }
    }

    
    // Initialize plans:
    for(int j=0; j < F_DENSE_PTS+2; j++){
        current_plan[j].cost = -1;
        current_plan[j].cells = new size_t[MAX_POLY_CELLS];
        
        tmp_path[j].cost=-1;
        tmp_path[j].cells = new size_t[MAX_POLY_CELLS];
    }

    exploration_map.X_dim = map_x;
    exploration_map.Y_dim = map_y;
    exploration_map.Z_dim = 1;
    exploration_map.data = new char[map_x*map_y*1];
    
    // Initialize to "unknown"
    for(int i = 0; i < map_x; i++)
        for(int j = 0; j < map_y; j++)
            exploration_map.set(i, j, 0, MAP_UNKN);
    
    // Fast planner used most of the time
    F_sparse = new FMT(map_x, map_y, 0, F_SPARSE_PTS, KNN_CON, 25, &exploration_map);
    f_sparse_calls = 0;
    // Slow(er) planner used when F_sparse fails
    F_dense  = new FMT(map_x, map_y, 0, F_DENSE_PTS, KNN_CON, 75, &exploration_map);
    f_dense_calls = 0;

}

Explorer::~Explorer(){
    F_dense->~FMT();
    F_sparse->~FMT();
    delete[] exploration_map.data;
    for(int i = 0; i < MAX_CLUSTERS; i++){
        delete [] clusters[i].members;
    }
    
    for(int i = 0; i < F_DENSE_PTS+2; i++){
        if( current_plan[i].cells != nullptr)
            delete [] current_plan[i].cells;
        if( tmp_path[i].cells != nullptr)
            delete [] tmp_path[i].cells;
    }
    delete [] frontiers;
    delete[] cluster_ids;
}


int Explorer::get_clusters(Point* cluster_centroids, float* costs, int max_num_clusters){
    int i = 0;
    for(i = 0; i < num_active_clusters && i < max_num_clusters; i++){
        costs[i] = goal_costs[self_id][i];
        cluster_centroids[i].x = clusters[active_clusters[i]].median.x;
        cluster_centroids[i].y = clusters[active_clusters[i]].median.y;
        cluster_centroids[i].z = clusters[active_clusters[i]].median.z;
    }
    return i;
}

void Explorer::get_neighbor_clusters(int agent_id, Point *cluster_centroids, float *costs, int n_clusters){
    // Test whether to add to list:
    if(agent_id == self_id){
        return;
    }
    
    for(int c_new = 0; c_new < n_clusters; c_new++){
        bool redundant = false;
        int new_old_index = 0;
        
        for(int c_old = 0; c_old < num_goals; c_old++){
            if(dist(cluster_centroids[c_new], goal_locations[c_old]) < merge_dist && !collision(cluster_centroids[c_new],goal_locations[c_old], &exploration_map)){
                redundant = true;
                new_old_index = c_old;
                break;
            }
        }
        if(redundant){
            //Add neighbor's cost to appropriate row
            goal_costs[agent_id][new_old_index] = costs[c_new];
        }else{
#ifdef EXPLORE_DEBUG
            printf("Adding a goal at (%f,%f,%f) with cost %f (agent %d)\n", cluster_centroids[c_new].x,cluster_centroids[c_new].y,cluster_centroids[c_new].z,costs[c_new],self_id);
#endif
            // Add new goal and fill in inf for unknown rewards
            // Put in inf for yourself if you didn't originally plan for it.
            // If you trust map merging then you can remove this and plan based on others' goals.
            goal_locations[num_goals].x = cluster_centroids[c_new].x;
            goal_locations[num_goals].y = cluster_centroids[c_new].y;
            goal_locations[num_goals].z = cluster_centroids[c_new].z;
//            compute_cost(swarm_locs[self_id], goal_locations[num_goals], &(paths[num_goals]), num_goals);
            goal_costs[self_id][num_goals]=1000000;
//            if(paths[num_goals][0].x == -1){
//                goal_costs[self_id][num_goals] = 1000000;
//            }else{
//                goal_costs[self_id][num_goals] = paths[num_goals][0].y;
//            }
            for(int i = 0; i < team_size; i++){
                if(i != self_id){
                    if(i != agent_id){
                        goal_costs[i][num_goals] = 1000000;
                    }else {
                        goal_costs[i][num_goals] = costs[c_new];
                    }
                }
            }
            num_goals++;
        }
    }
}

// Finishes assignments
void Explorer::assign(Point* waypoint){
    
    const int PLANNER_NOMINAL = 1;
    const int PLANNER_DEFAULT = 2;
    const int PLANNER_STUCK   = 3;
    int planner_status = PLANNER_NOMINAL;
    
#ifdef EXPLORE_DEBUG
    printf("Assigning %d agents to %d goals:\n", team_size, num_goals);
    printf("Agent %d  at (%f,%f,%f)\n", self_id, swarm_locs[self_id].x, swarm_locs[self_id].y, swarm_locs[self_id].z);
    printf("Goals are at \n");
    for(int i = 0; i < num_goals; i++)
        printf("(%f,%f,%f), %f \n", goal_locations[i].x,goal_locations[i].y,goal_locations[i].z, goal_costs[self_id][i]);
#endif
    if(num_goals > 0){
        // give yourself a discount on the goal you're headed toward:
        Point goal_pt;
        int goal_ind = -1;
        
        int target_assignments[team_size];
        for(int i =0; i < team_size; i++)
            target_assignments[i]=0;
        
        if(team_size > 1){
            std::vector<std::vector<double>> tmp_costs (team_size, std::vector<double>( num_goals, 0));
            for(int agent = 0; agent < team_size; agent++){
                for(int goal =0; goal< num_goals; goal++){
                    tmp_costs[agent][goal] = goal_costs[agent][goal];
                }
            }
    #ifdef EXPLORE_DEBUG
            printf("Agent %d cost matrix:\n",self_id);
            for(int i = 0; i < team_size; i++){
                printf("%d : ",i);
                for(int j=0; j <num_goals; j++){
                    printf(" %2.2f ", goal_costs[i][j]);
                }
                printf("\n");
            }
            
    #endif
            
            assign_min(tmp_costs, target_assignments, team_size);

                
            }
        
        if(team_size == 1 || target_assignments[self_id] >= num_goals){
            planner_status = PLANNER_DEFAULT;
            // Go to minimum:
            float mincost = goal_costs[self_id][0];
            int min_goal  = 0;
            for(int i = 1; i < num_goals; i++){
                if(goal_costs[self_id][i] < mincost){
                    mincost =goal_costs[self_id][i];
                    min_goal = i;
                }
            }
            target_assignments[self_id] = min_goal;
        }
        
        goal_ind = target_assignments[self_id];
        goal_pt.x = goal_locations[target_assignments[self_id]].x;
        goal_pt.y = goal_locations[target_assignments[self_id]].y;
        goal_pt.z = goal_locations[target_assignments[self_id]].z;
        if(exploration_map(int(goal_pt.x), int(goal_pt.y), int(goal_pt.z)) >=MAP_FREE){
            exploration_map.set(int(goal_pt.x), int(goal_pt.y), int(goal_pt.z), MAP_GOAL_U);
        }else{
            exploration_map.set(int(goal_pt.x), int(goal_pt.y), int(goal_pt.z), MAP_GOAL_F);
        }


#ifdef EXPLORE_DEBUG
        printf("Agent %d assigned to (%f,%f)\n", self_id, goal_pt.x,goal_pt.y);
#endif
        
        // Always replan.
        F_sparse->set_initial_state(&current_vel, &current_acc);
        F_dense->set_initial_state(&current_vel, &current_acc);

        curr_segs = F_sparse->fmtstar(swarm_locs[self_id], goal_pt, &exploration_map, current_plan);
        if(current_plan[0].cost == -1){
            curr_segs = F_dense->fmtstar(swarm_locs[self_id], goal_pt, &exploration_map, current_plan);
        }

        
#ifdef EXPLORE_DEBUG
        /*
        printf("Path has waypoints: (incorrect if reversed)");
        for(int seg = 0; seg <= curr_segs; seg++)
            printf(" ( %f,%f,%f) ", current_plan[seg].coefficients_x[0], current_plan[seg].coefficients_y[0], current_plan[seg].coefficients_z[0]);
        printf("\n");
        
        int k = 1;
        printf("Coefficients of current polynomial (order %d, rev = %d):\n", current_plan[k].order,current_plan[k].reverse);
        printf("x=[");
        for(int i = 0;  i< current_plan[k].order; i++)
            printf(" %f, ",current_plan[k].coefficients_x[i]);
        printf("]\ny=[");
        for(int i = 0;  i< current_plan[k].order; i++)
            printf(" %f, ",current_plan[k].coefficients_y[i]);
        printf("]\nz=[");
        for(int i = 0;  i< current_plan[k].order; i++)
            printf(" %f, ",current_plan[k].coefficients_z[i]);
        printf("]\n");
         */
#endif
        
        if(current_plan[0].cost == -1){
            current_vel.x = 0.0; current_vel.y = 0.0; current_vel.z = 0.0;
            current_acc.x = 0.0; current_acc.y = 0.0; current_acc.z = 0.0;
            
            planner_status=PLANNER_STUCK;
            waypoint->x = swarm_locs[self_id].x;
            waypoint->y = swarm_locs[self_id].y;
            waypoint->z = swarm_locs[self_id].z;
            
        }else{ // this is looking at the next point.
            if(curr_segs > 1){
                waypoint->x =current_plan[2].coefficients_x[0];
                waypoint->y =current_plan[2].coefficients_y[0];
                waypoint->z =current_plan[2].coefficients_z[0];
            }else{
                get_poly_der(current_plan[1], current_plan[1].duration, waypoint, 0);
            }
        }
        
        get_poly_der(current_plan[curr_segs], current_plan[curr_segs].duration, waypoint, 0);

    
#ifdef EXPLORE_STATUS
    if(planner_status == PLANNER_NOMINAL){
        printf(ANSI_COLOR_GREEN "AGENT %d IS GOOD. GOING TO GOAL\n" ANSI_COLOR_RESET, self_id);
    }else if(planner_status == PLANNER_DEFAULT){
        printf(ANSI_COLOR_YELLOW "AGENT %d NOT ASSIGNED. BEING GREEDY\n" ANSI_COLOR_RESET, self_id);
    }else{
        printf(ANSI_COLOR_RED "AGENT %d IS STUCK. STAYING STILL\n" ANSI_COLOR_RESET, self_id);
    }
#endif
        
    }
}


// Updates the map markers
void Explorer::update_state(Point* positions, int num_agents){
    
    // THIS IS DUMB
    // Clear the unknown states on the map, reset all goals etc.
    for(int i = 0; i < exploration_map.X_dim; i++){
        for(int j= 0; j < exploration_map.Y_dim; j++){
            char state = exploration_map(i,j,0);
            if(state==MAP_GOCC_U || state == MAP_GOAL_U){
                exploration_map.set(i,j,0,MAP_UNKN);
            }
            if(state==MAP_GOCC_F || state == MAP_GOAL_F){
                exploration_map.set(i,j,0,MAP_FREE);
            }
        }
    }
    

    team_size = num_agents;
    int clearing_radius = 3;
    // Mark other positions, self position, free/occupied/unknown cells
    for(int i = 0; i < team_size; i++){
        exploration_map.set(int(swarm_locs[i].x), int(swarm_locs[i].y), int(swarm_locs[i].z), MAP_TRACE);
        swarm_locs[i].x = positions[i].x;
        swarm_locs[i].y = positions[i].y;
        swarm_locs[i].z = positions[i].z;
        if(i == self_id){
            exploration_map.set(int(swarm_locs[i].x), int(swarm_locs[i].y), int(swarm_locs[i].z), MAP_SELF);
        }else{
            // Make new marks
            exploration_map.set(int(swarm_locs[i].x), int(swarm_locs[i].y), int(swarm_locs[i].z), MAP_OTHER);
            for(int di = -clearing_radius; di <= clearing_radius; di++){
                for(int dj = -clearing_radius; dj <= clearing_radius; dj++){
                    if(exploration_map(int(swarm_locs[i].x)+di, int(swarm_locs[i].y)+dj, int(swarm_locs[i].z)) == MAP_UNKN)
                        exploration_map.set(int(swarm_locs[i].x)+di, int(swarm_locs[i].y)+dj, int(swarm_locs[i].z),MAP_GOCC_U);
                    if(exploration_map(int(swarm_locs[i].x)+di, int(swarm_locs[i].y)+dj, int(swarm_locs[i].z)) <= MAP_FREE)
                        exploration_map.set(int(swarm_locs[i].x)+di, int(swarm_locs[i].y)+dj, int(swarm_locs[i].z),MAP_GOCC_F);
                }
            }

        }
    }
}

void Explorer::merge_map(Map* other_map){
    for(int i = 0; i < other_map->X_dim; i++){
        for(int j = 0; j < other_map->Y_dim; j++){
            char state = exploration_map(i,j,0);
            if(state == MAP_UNKN || state == MAP_GOCC_U || state == MAP_GOAL_U){ //all unknown states
                if((*other_map)(i,j,0) == MAP_UNKN){
                    
                }else{
                    if((*other_map)(i,j,0) <= MAP_FREE){
                        exploration_map.set(i,j,0,MAP_FREE);
                    }else if((*other_map)(i,j,0) == MAP_OCCU){
                        exploration_map.set(i,j,0, MAP_OCCU);
                    }
                }
            }
        }
    }
}

void Explorer::update_map(Map* measurement, Point lower, Point upper){
    int num_cells = 0;
    int num_skipped = 0;
    for(int i = lower.x; i <= upper.x; i++){
        for(int j = lower.y; j <= upper.y; j++){
            if(exploration_map(i,j,0) == MAP_UNKN || exploration_map(i,j,0)==MAP_GOCC_U){
                if((*measurement)(i,j,0) > MAP_UNKN){
                    exploration_map.set(i,j,0, MAP_OCCU);
                }else if((*measurement)(i,j,0) < MAP_UNKN){
                    if(exploration_map(i,j,0)==MAP_GOCC_U){
                        exploration_map.set(i,j,0, MAP_GOCC_F);
                    }else{
                        exploration_map.set(i,j,0,MAP_FREE);
                    }
                }
                num_cells++;
            }else{
                num_skipped++;
            }
        }
    }
    
#ifdef EXPLORE_DEBUG
//    printf("Agent %d: Skipped %d cells, updated %d cells!\n",self_id, num_skipped, num_cells);
#endif
}


//Overwrites map with given value
void Explorer::force_map_update(Map* map){
    exploration_map.X_dim = map->X_dim;
    exploration_map.Y_dim = map->Y_dim;
    exploration_map.Z_dim = map->Z_dim;
    
    for(int i = 0; i < exploration_map.X_dim*exploration_map.Y_dim*exploration_map.Z_dim; i++)
        exploration_map.data[i] = map->data[i];
    
    exploration_map.set(int(swarm_locs[self_id].x), int(swarm_locs[self_id].y), int(swarm_locs[self_id].z), MAP_SELF);
}

void Explorer::clear_map(){
    
    for(int i = 0; i < exploration_map.X_dim*exploration_map.Y_dim*exploration_map.Z_dim; i++)
        exploration_map.data[i] = MAP_UNKN;
    
    exploration_map.set(int(swarm_locs[self_id].x), int(swarm_locs[self_id].y), int(swarm_locs[self_id].z), MAP_SELF);
}



bool Explorer::cluster_frontiers(){
    if(team_size > MAX_TEAMSIZE){
        team_size = MAX_TEAMSIZE;
        printf("Error: Team_size exceeded maximum!\n");
    }
    F_sparse->set_initial_state(&current_vel, &current_acc);
    F_dense->set_initial_state(&current_vel, &current_acc);

    
// What should happen: Only reset nodes that might have been changed in the map
// Measurement should have bounds so we don't traverse the whole map
    
    // reset all cluster ids:
    for(int i = 0; i < exploration_map.X_dim*exploration_map.Y_dim; i++)
        cluster_ids[i] = -1;

    // Can do this more efficiently if so inclined
    // reset frontiers on map:
    for(int i = 0; i < num_frontiers; i++)
        if(exploration_map.data[frontiers[i]] == MAP_FRONT)
            exploration_map.data[frontiers[i]] = MAP_FREE;
    num_frontiers = 0;
    
    for(int i = 0; i < num_clusters; i++)
        clusters[i].size = 0;
    num_clusters = 0;
    
    
    
    // Expand map
    const int inflation_radius = 1;
    for(int i = 0; i <  exploration_map.X_dim; i++){
        for(int j = 0; j < exploration_map.Y_dim; j++){
            if(exploration_map(i,j,0)<  MAP_OCC_THRESH){
                 // set all valid neighbors to MAP_GOCCU
            // Check neighbors and set self to occu if neighbor is actually occupied.
                 for(int di = -inflation_radius; di <= inflation_radius; di++){
                     if(i + di >= 0 && i+di < exploration_map.X_dim){
                         for(int dj = -inflation_radius; dj<=inflation_radius; dj++){
                             if(j + dj >= 0 && j+dj < exploration_map.Y_dim){
                                 if(!(di==0 && dj==0) && exploration_map(i+di, j+dj, 0) > MAP_OCC_THRESH){
                                     if(exploration_map(i+di, j+dj, 0) != MAP_GOCC_F && exploration_map(i+di, j+dj, 0) != MAP_GOCC_U ){
                                         if(exploration_map(i,j,0) > MAP_FREE){
                                             exploration_map.set(i, j, 0, MAP_GOCC_U);
                                         }else{
                                             exploration_map.set(i, j, 0, MAP_GOCC_F);
                                         }
                                     }
                                 }
                             }
                         }
                     }
                 }
            }
        }
    }
    
    // First take a pass through the map and extend any lines of occupied cells
    int min_line = 3;
    bool is_line;
    for(int j = 0; j < exploration_map.Y_dim; j++){
        is_line = false;
        int num_consecutive = 0;
        for(int i = 0; i < exploration_map.X_dim; i++){
            if(! is_line){
                char state = exploration_map(i,j,0);
                if(state == MAP_OCCU){
                    num_consecutive++;
                }else{
                    num_consecutive=0;
                }
                // Check if we just became a line
                if(num_consecutive>min_line){
                    num_consecutive = 0;
                    is_line= true;
                    //Backpropogate until empty space is found:
                    for(int back_i = i; back_i >= 0; back_i--){
                        state = exploration_map(back_i,j,0);
                        if(state > MAP_FREE){
                            if(state == MAP_UNKN)
                                exploration_map.set(back_i,j,0,MAP_GOCC_U);
                        }else{
                            break;
                        }
                    }
                }
            }else{
                if(exploration_map(i,j,0) > MAP_FREE){
                    if(exploration_map(i,j,0)==MAP_UNKN)
                        exploration_map.set(i,j,0,MAP_GOCC_U);
                }else{
                    is_line = false;
                }
            }
        }
    }
    for(int i = 0; i < exploration_map.X_dim; i++){
        is_line = false;
        int num_consecutive = 0;
        for(int j = 0; j < exploration_map.Y_dim; j++){
            if(! is_line){
                char state = exploration_map(i,j,0);
                if(state == MAP_OCCU){
                    num_consecutive++;
                }else{
                    num_consecutive=0;
                }
                // Check if we just became a line
                if(num_consecutive>min_line){
                    num_consecutive = 0;
                    is_line= true;
                    //Backpropogate until empty space is found:
                    for(int back_j = j; back_j >= 0; back_j--){
                        if(exploration_map(i,back_j,0) > MAP_FREE){
                            if(exploration_map(i,back_j,0) == MAP_UNKN)
                                exploration_map.set(i,back_j,0,MAP_GOCC_U);
                        }else{
                            break;
                        }
                    }
                }
            }else{
                if(exploration_map(i,j,0) > MAP_FREE){
                    if(exploration_map(i,j,0) == MAP_UNKN)
                        exploration_map.set(i,j,0,MAP_GOCC_U);
                }else{
                    is_line = false;
                }
            }
        }
    }
    
    
#ifdef EXPLORE_DEBUG
    printf("Clustering points.\n");
#endif
    int cluster_radius = 3;
    // Step 0: find frontiers/connected components
    for(int i = 0; i < exploration_map.X_dim; i++){
        for(int j = 0; j < exploration_map.Y_dim; j++){
            // Check for frontier-ness:

            if(exploration_map(i,j,0) <= MAP_FREE ){
                bool is_frontier = false;
                // Check neighbors:
                for(int di = -1; di <= 1; di++){
                    if(is_frontier)
                        break;
                    for(int dj = -1; dj <= 1; dj++){
                        if(exploration_map(i+di, j+dj, 0) == MAP_UNKN ){
                            is_frontier = true;
//                            break;
                        }else if(exploration_map(i+di,j+dj,0)>=MAP_OCC_THRESH){
                            is_frontier = false;
                            break;
                        }
                    }
                }

                if(is_frontier){
                    exploration_map.set(i, j, 0, MAP_FRONT);
                    if(num_frontiers < exploration_map.X_dim*exploration_map.Y_dim){
                        int f = (int)exploration_map.ind2num(i, j, 0);
                        frontiers[num_frontiers] = f;
                        num_frontiers+=1;
                        // TODO: Half of these look forward to where we know there aren't frontiers. can trim
                        // Backtrack to see if you're near an existing frontier:
                        for(int di = -cluster_radius; di< cluster_radius; di++){
                            if(cluster_ids[f] != -1) // Found a cluster
                                break;
                            for(int dj = -cluster_radius; dj < cluster_radius; dj++){
                                if(cluster_ids[f] != -1) // Found a cluster
                                    break;
                                if(!(di == 0 && dj == 0) && (i+di >= 0 && i+di < exploration_map.X_dim && j+dj >= 0 && j+dj < exploration_map.Y_dim)){
                                    if(exploration_map(i+di, j+dj, 0) == MAP_FRONT){
                                        int f_n = (int)exploration_map.ind2num(i+di, j+dj, 0);
                                        if(f < f_n)
                                            printf("Assumption that frontiers are added in order is violated\n");

                                        if(cluster_ids[f_n] == -1){
#ifdef EXPLORE_DEBUG
                                            printf("Error - frontier assigned to cluster -1");
#endif
                                        }else{
                                            cluster_ids[f] = cluster_ids[f_n];
                                            clusters[cluster_ids[f_n]].members[clusters[cluster_ids[f_n]].size] = f;
                                            clusters[cluster_ids[f_n]].size++;
                                        }
#ifdef EXPLORE_DEBUG
#endif
                                    }
                                }
                            }
                        }
                        // If no neighbor found, assign a new cluster
                        if(cluster_ids[f]==-1){
                            cluster_ids[f] = num_clusters;
                            clusters[cluster_ids[f]].size = 1;
                            clusters[cluster_ids[f]].number = num_clusters;
                            clusters[cluster_ids[f]].members[0] = f;
                            num_clusters++;
                        }
                    }else{
#ifdef EXPLORE_DEBUG
                        printf("ERROR: Adding frontiers that shouldn't exist!\n");
#endif
                    }
                }
            }
        }
    }
    
    
    if(num_frontiers == 0){
#ifdef EXPLORE_DEBUG
        printf("Exploration complete!\n");
#endif
        return true;
    }
#ifdef EXPLORE_DEBUG
        printf("Found %d frontiers!\n", num_frontiers);
#endif
    
    //Compute median. We pick this as the location of a frontier to avoid moving the median into occupied/unknown space.
   for(int i = 0; i < num_clusters; i++){
       if(clusters[i].size!=0){
           int med_index = int(clusters[i].size/2);
           Point median;
           exploration_map.num2ind(&median, clusters[i].members[med_index]);
           clusters[i].median.x = median.x;
           clusters[i].median.y = median.y;
           clusters[i].median.z = median.z;
       }
    }
    
    // Time to here: ~0.1 ms
    
#ifdef EXPLORE_DEBUG
    printf("Formed %d clusters at ", num_clusters);
    for(int i = 0; i < num_clusters; i++)
        printf(" (%f %f %f) ", clusters[i].median.x, clusters[i].median.y, clusters[i].median.z);
    printf("\n");
#endif
    
    // If we have more clusters than agents, try to condense them
    int merge_partners[num_clusters];
    for(int i = 0; i < num_clusters; i++)
        merge_partners[i] = -1;
    if(num_clusters > team_size){
        
        // Step 2: Collision check each pair of components/merge if OK
        for(int c1 = num_clusters-1; c1 >= 0; c1--){
            if(clusters[c1].size != 0){
                for(int c2 = num_clusters-1; c2 >=0; c2--){
                    bool merge = false;

                    
                    if(c2 != c1 && clusters[c2].size != 0){
                        
                        // Check if should merge based on centroids:
                        if(collision(clusters[c1].median,clusters[c2].median , &exploration_map)){
                            if(dist(clusters[c1].median, clusters[c2].median) < merge_dist){
                                merge = true;
                            }
                        }
                        
                        // Check if should merge based on ends
                        if(!merge){
                            Point p1; Point p2;
                            int c1_step = 1;
                            int c2_step = 1;
                            if(clusters[c1].size > 1)
                                c1_step = int(clusters[c1].size-1);
                            if(clusters[c2].size > 30)
                                c2_step = int(clusters[c2].size-1);
                            
                            for(int c1_m = 0; c1_m < clusters[c1].size; c1_m+=c1_step){
                                if(merge)
                                    break;
                                for(int c2_m = 0; c2_m < clusters[c2].size; c2_m+=c2_step){
                                    exploration_map.num2ind(&p1, clusters[c1].members[c1_m]);
                                    exploration_map.num2ind(&p2, clusters[c2].members[c2_m]);
                                    if(dist(p1, p2)<= merge_dist && collision(p1, p2, &exploration_map)){
                                        merge = true;
                                        break;
                                    }
                                }
                            }
                        }
                        
                        if(merge)
                            merge_partners[c1] = c2;
                    }
                }
            }
        }
        
        for(int c1 = num_clusters-1; c1 >=0; c1--){
            if(clusters[c1].size != 0 && merge_partners[c1]!=-1){
                int c2 = merge_partners[c1];
                while(merge_partners[c2] != -1 && merge_partners[c2]<c2){
                    c2 = merge_partners[c2];
                }

                if(c1!=c2){
                    // Relabel the frontier container
                    relabel_frontiers(c1, c2);
                    merge_clusters(&clusters[c1], &clusters[c2]);
                }
            }
        }

    }
    
    for(int c = 0; c < MAX_CLUSTERS; c++){
        if(clusters[c].size != 0){
            Point cent;
            Point loc;
            exploration_map.num2ind(&cent, clusters[c].members[0]);
            cent.x /= float(clusters[c].size);
            cent.y /= float(clusters[c].size);
            cent.z /= float(clusters[c].size);
            for(int i = 1; i< clusters[c].size; i++){
                exploration_map.num2ind(&loc, clusters[c].members[i]);
                cent.x += loc.x/(float(clusters[c].size));
                cent.y += loc.y/(float(clusters[c].size));
                cent.z += loc.z/(float(clusters[c].size));
            }
            if(exploration_map(int(cent.x),int(cent.y),int(cent.z)) == MAP_FREE){
                clusters[c].median.x = round(cent.x);
                clusters[c].median.y = round(cent.y);
                clusters[c].median.z = round(cent.z);
            }
        }
    }


    num_active_clusters = 0;
    for(int i = 0; i < MAX_CLUSTERS; i++){
        if(clusters[i].size > 0){
            active_clusters[num_active_clusters] = i;
            num_active_clusters++;
            // If centroid lands on top of agent, move to the frontier.
            if(dist(swarm_locs[self_id], clusters[i].median) < 2.0){
                int mid_ind = clusters[i].size/2;
                Point mid_pt;
                exploration_map.num2ind(&mid_pt, clusters[i].members[mid_ind]);
                clusters[i].median.x =mid_pt.x;
                clusters[i].median.y =mid_pt.y;
                clusters[i].median.z =mid_pt.z;
                
            }
        }
    }
    
#ifdef EXPLORE_DEBUG
    printf("Post-processed cluster locations: ", num_active_clusters);
    for(int i = 0; i < num_active_clusters; i++)
        printf(" (%f %f %f) ", clusters[active_clusters[i]].median.x, clusters[active_clusters[i]].median.y, clusters[active_clusters[i]].median.z);
    printf("\n");

#endif
    
    // Exploration not done yet
    return false;
}

/*
void Explorer::compute_cost(Point start, Cluster* goal, PolyState* path, int goal_num){
    F_sparse->fmtstar(start, goal, &exploration_map, path);
    
    
    f_sparse_calls++;
    if(path[0].x == -1){
        F_dense->fmtstar(start, goal, &exploration_map, path);
        f_dense_calls++;
    }
    if(path[0].x == -1){
        goal_costs[self_id][goal_num] = 1000000;
    }else{
        goal_costs[self_id][goal_num] = path[0].y;
    }

}
 */

void Explorer::compute_cost(Point start, Point goal, Point goal_vel, PolyState* path, int goal_num){

    // Check for line of sight - we don't actually use these paths for anything.
    
    path_segs[goal_num] = F_sparse->fmtstar(start, goal, goal_vel, &exploration_map, path);
    

    f_sparse_calls++;
    if(path[0].cost == -1){
        path_segs[goal_num] = F_dense->fmtstar(start, goal, goal_vel,  &exploration_map, path);
        f_dense_calls++;
    }
    if(path[0].cost == -1){
#ifdef EXPLORE_DEBUG
        printf("Couldn't find a path to goal (%f,%f)\n", goal.x, goal.y);
#endif
        goal_costs[self_id][goal_num] = 1000000;
    }else{
#ifdef EXPLORE_DEBUG
        printf("Cost to goal (%f,%f) is %f\n", goal.x, goal.y, path[0].cost);
#endif
        goal_costs[self_id][goal_num] = path[0].cost;
    }
}




int Explorer::compute_costs(){
    num_goals = 0;
    
#ifdef EXPLORE_DEBUG
    printf("Computing costs to %d clusters\n", num_active_clusters);
#endif
    // Always clear cache. TODO: Fix this
    F_sparse->clear_cache();
    F_dense->clear_cache();
    

    for(int i = 0; i < num_active_clusters; i++){
        goal_locations[i].x = clusters[active_clusters[i]].median.x;
        goal_locations[i].y = clusters[active_clusters[i]].median.y;
        goal_locations[i].z = clusters[active_clusters[i]].median.z;
        num_goals+=1;
        float min_dist = 0;
        Point min_dir;
        Point loc;
        for(int c = 0; c < clusters[active_clusters[i]].size; c++){
            exploration_map.num2ind(&loc, clusters[active_clusters[i]].members[i]);
            if(dist(loc, goal_locations[i]) < min_dist){
                min_dir.x = loc.x - goal_locations[i].x;
                min_dir.y = loc.y - goal_locations[i].y;
                min_dir.z = loc.z - goal_locations[i].z;
            }
        }

        if(!exploration_map.is_free(goal_locations[i])){
            min_dir.x = -min_dir.x;
            min_dir.y = -min_dir.y;
            min_dir.z = - min_dir.z;
        }
        
        compute_cost(swarm_locs[self_id], clusters[active_clusters[i]].median, min_dir, tmp_path, i);
    }
    
    // Discount the goal you already chose
    /*
    int goal_ind = -1;
    if(current_plan[0].cost != -1){
        
        // Don't discount staying still.
            float mindist = MAXFLOAT;
            for(int i =0; i < num_goals; i++){
                float d= dist(goal_locations[i],target);
                if(d > merge_dist && d < mindist){
                    mindist = d;
                    goal_ind = i;
                }
            }
            
        if(goal_ind != -1){
            if(dist(target, swarm_locs[self_id]) > 1.0){
                goal_costs[self_id][goal_ind] -= 2000;
            }else{
                // Discourage staying still
                goal_costs[self_id][goal_ind] += 2000;
            }
        }
    }
     */
    

    return num_active_clusters;
}

                                                          


// Relabels members of the cluster with larger index to have smaller index.
void Explorer::relabel_frontiers(int clust1_ind, int clust2_ind){
    int other_clust = clust2_ind;
    int target_clust = clust1_ind;
    if(clusters[clust1_ind].number > clusters[clust2_ind].number){
        other_clust = clust1_ind;
        target_clust = clust2_ind;
    }
    for(int i = 0; i < clusters[other_clust].size; i++)
        cluster_ids[clusters[other_clust].members[i]] = clusters[target_clust].number;
}


void Explorer::print_costs(){
//    for(int i = 0; i < num_active_clusters; i++){
//        if(paths[i][0].cost != -1){
//            printf("cost to (%f,%f) is %f with %d segments\n", clusters[active_clusters[i]].median.x, clusters[active_clusters[i]].median.y, paths[i][0].cost, path_segs[i]);}
//    }
}

#ifdef EXPLORE_PLOT
// Plots clusters using Ascii
void Explorer::plot_clusters(){
    printf("\n");
    for(int j = 0; j < exploration_map.Y_dim; j++){
        for(int i = 0; i < exploration_map.X_dim; i++){
            if(exploration_map(i,j,0) == MAP_FREE){
                printf("  " );
            }else if(exploration_map(i,j,0) == MAP_GOCC_U || exploration_map(i,j,0) == MAP_GOCC_F ){
                printf(ANSI_COLOR_GRAY "##" ANSI_COLOR_RESET);
            }else if(exploration_map(i,j,0) == MAP_FRONT){
                if(cluster_ids[exploration_map.ind2num(i, j, 0)]<100){
                    if(cluster_ids[exploration_map.ind2num(i, j, 0)]%6==5){
                        printf(ANSI_COLOR_RED "**" ANSI_COLOR_RESET);
                    }else if(cluster_ids[exploration_map.ind2num(i, j, 0)]%6==1){
                        printf(ANSI_COLOR_BLUE "**" ANSI_COLOR_RESET);
                    }else if(cluster_ids[exploration_map.ind2num(i, j, 0)]%6==2){
                        printf(ANSI_COLOR_YELLOW "**" ANSI_COLOR_RESET);
                    }else if(cluster_ids[exploration_map.ind2num(i, j, 0)]%6==3){
                        printf(ANSI_COLOR_GREEN "**" ANSI_COLOR_RESET);
                    }else if(cluster_ids[exploration_map.ind2num(i, j, 0)]%6==4){
                        printf(ANSI_COLOR_MAGENTA "**" ANSI_COLOR_RESET);
                    }else if(cluster_ids[exploration_map.ind2num(i, j, 0)]%6==0){
                        printf(ANSI_COLOR_CYAN "**" ANSI_COLOR_RESET);
                    }else{
                        printf(" %d", cluster_ids[exploration_map.ind2num(i, j, 0)]);
                    }
                }else{
                    printf(" %d", cluster_ids[exploration_map.ind2num(i, j, 0)]);
                }
            }else if(exploration_map(i,j,0)== MAP_OCCU){
                printf("‡‡");
            }else if(exploration_map(i,j,0) == MAP_UNKN){
                printf(ANSI_BG_CYAN "  " ANSI_COLOR_RESET);
            }else if(exploration_map(i,j,0) == MAP_GOAL){
                printf(ANSI_BG_RED "  " ANSI_COLOR_RESET);
                exploration_map.set(i, j, 0, MAP_FREE);
            }else if(exploration_map(i,j,0) == MAP_PLAN){
                printf(ANSI_BG_GREEN "  " ANSI_COLOR_RESET);
                exploration_map.set(i, j, 0, MAP_FREE);
            }else if(exploration_map(i,j,0) == MAP_CENT){
                printf(ANSI_COLOR_BLUE "XX" ANSI_COLOR_RESET);
            }else if(exploration_map(i,j,0) == MAP_SELF){
                printf(ANSI_COLOR_GREEN "◊◊" ANSI_COLOR_RESET);
            }else if(exploration_map(i,j,0) == MAP_OTHER){
                printf(ANSI_COLOR_YELLOW "◊◊" ANSI_COLOR_RESET);
            }else if(exploration_map(i,j,0) == MAP_TRACE){
                printf(ANSI_BG_BLUE "  " ANSI_COLOR_RESET);
            }
            
        }
        printf("\n");
    }
}
#endif

void Explorer::get_vel(Point *vel){
    vel->x = current_vel.x;
    vel->y = current_vel.y;
    vel->z = current_vel.z;
}

void Explorer::get_acc(Point *acc){
    acc->x = current_acc.x;
    acc->y = current_acc.y;
    acc->z = current_acc.z;
}


void Explorer::print_call_rates(){
    printf("Sparse calls: %d  Dense calls: %d  Percent of first-successes: %f\n", f_sparse_calls, f_dense_calls, float(f_sparse_calls)/(float(f_sparse_calls+f_dense_calls)));
    printf("Sparse stats: ");
    F_sparse->print_stats();
    printf("Dense stats:");
    F_dense->print_stats();
}

void Explorer::print_cluster(Cluster* clust){
    printf("Cluster %d has %d members and centroid (%f,%f,%f)\n", clust->number, clust->size,clust->median.x,clust->median.y,clust->median.z);
}
