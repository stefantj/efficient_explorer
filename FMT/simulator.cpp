//
//  simulator.cpp
//  FMT
//
//  Created by Megamind on 7/22/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#include "simulator.h"
#include "Explorer.h"
#include <stdlib.h>
#include "afghan_village.h"

Simulator::Simulator(){
    delta_time = 0.1;
    plot_delay = 100;
}

void Simulator::init_video_file(FILE* f){
    fprintf(f,"println(\"Loading PyPlot\");\n");
    fprintf(f,"using PyPlot;\n");
    fprintf(f,"include(\"afghan_map.jl\")\n");
    fflush(f);
}

void Simulator::print_video_file(FILE* f, std::string filename, float time){
    fprintf(f, "include(\"%s.jl\")\n", filename.c_str());
    fprintf(f, "figure(1, figsize=(16,16));clf(); axis([0,400,0,300]);\n");
 //   fprintf(f, "imshow(afghan_map', cmap=\"gray\", interpolation=\"none\");\n");
    fprintf(f,"imshow(global_map', cmap=\"gray\", interpolation=\"none\");\n");
    fprintf(f,"axis(\"equal\")\n");
    
    fprintf(f,"scatter(targets[:,1], targets[:,2], marker=\"x\", color=:red)\n");
    /*
    for(int i =0; i < num_agents; i++){
        if(agents[i].planner->current_plan[0].x > 1){
            fprintf(f,"plot(trajectory_%d[2:end,1], trajectory_%d[2:end,2], color=:red, linestyle=\"-.\");\n",i,i);
        }
        for(int j = 0; j < num_agents; j++){
            if(i > 0 && comms[j][i] == true){
                fprintf(f,"plot([locations[%d,1],locations[%d,1]], [locations[%d,2],locations[%d,2]],color=:green,linestyle=\"-.\");\n",j+1,i+1,j+1,i+1);
            }
        }
    }
    */
    for(int i = 0; i < num_agents; i++){
        if(i%4 == 0){
            fprintf(f,"scatter(locations[%d,1], locations[%d,2], marker=\"o\", color=:green)\n",i+1,i+1);
        }else if(i%3 == 0){
            fprintf(f,"scatter(locations[%d,1], locations[%d,2], marker=\"o\", color=:yellow)\n",i+1,i+1);
        }else if(i%2 == 0){
            fprintf(f,"scatter(locations[%d,1], locations[%d,2], marker=\"o\", color=:cyan)\n",i+1,i+1);
        }else{
            fprintf(f,"scatter(locations[%d,1], locations[%d,2], marker=\"o\", color=:orange)\n",i+1,i+1);
        }
    }
    fprintf(f, "title(\"t = %.2f\");\n", time);
    fprintf(f,"println(\"Plotting frame: %s\");\n", filename.c_str());
    fprintf(f,"savefig(\"%s.png\");\n", filename.c_str());
//    fprintf(f,"savefig(\"%s.png\", dpi=720);\n", filename.c_str());
    fflush(f);
}

void Simulator::run_simulator(int iters, int team_size){
    num_agents = team_size;
    if(team_size > MAX_NUM_AGENTS){
        printf("Error: team size too large! Maximum team is %d\n", MAX_NUM_AGENTS);
        return;
    }

    // Make map:
    true_map.X_dim = 400;
    true_map.Y_dim = 300;
    true_map.Z_dim = 1;
    true_map.data  = new char[400*300*1];
    global_map.X_dim = 400;
    global_map.Y_dim = 300;
    global_map.Z_dim = 1;
    global_map.data = new char[400*300*1];
    for(int i = 0; i < true_map.X_dim; i++){
        for(int j = 0; j < true_map.Y_dim; j++){
            global_map.set(i,j,0,MAP_UNKN);
            if(map[0][i][j] == 1){
                true_map.set(i, j, 0, MAP_OCCU);
            }else if(map[0][i][j]== 0){
                true_map.set(i,j,0,MAP_FREE);
            }else{
                printf("Warning: Map has invalid value %f at (%d,%d)\n", map[0][i][j], i,j);
            }
        }
    }
    
    
    
    
    /*
        Point path[100];
        Point start, goal;
        start.x = 75.0; start.y = 80.0; start.z = 0.0;
        goal.x = 200.0; goal.y = 200.0; goal.z = 0.0;
        
        FMT F_dense = FMT(true_map.X_dim, true_map.Y_dim, 0, F_DENSE_PTS, RAD_CON, 0);
        F_dense.clear_cache();
        F_dense.fmtstar(start, goal, &true_map, path);
    
    start.x = 75.0; start.y = 80.0; start.z = 0.0;
    goal.x = 50.0; goal.y = 100.0; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    start.x = 75.0; start.y = 80.0; start.z = 0.0;
    goal.x = 50.0; goal.y = 200.0; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    start.x = 75.0; start.y = 80.0; start.z = 0.0;
    goal.x = 150,0; goal.y = 110.0; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    start.x = 75.0; start.y = 80.0; start.z = 0.0;
    goal.x = 40.0; goal.y = 25.0; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 54; goal.y = 108; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 82; goal.y = 113; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 44.5; goal.y = 88.5; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 26.2; goal.y = 82.5; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 27; goal.y = 64; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 45.6; goal.y = 41.5; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 103.5; goal.y = 32.0; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 144; goal.y = 33; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 153; goal.y = 33; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 172; goal.y = 20.7; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 172; goal.y = 67.8; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);

    goal.x = 185; goal.y = 75; goal.z = 0.0;
    F_dense.clear_cache();
    F_dense.fmtstar(start, goal, &true_map, path);


    return;
*/

    for(int i = 0; i < num_agents; i++){
        printf("Making agent %d\n",i);
        agents[i].location.x = 4.0 +(i/2)*7.0;
        agents[i].location.y = 4.0 +(i%2)*7.0;
        agents[i].location.z = 0;
        agents[i].bearing = M_PI;
        agents[i].target.x = 300.0;
        agents[i].target.y = 300.0;
        agents[i].target.z = 0.0;
        agents[i].planner = new Explorer(i, true_map.X_dim, true_map.Y_dim);

        // Shouldnt be necessary anymore. For some reason this is necessary for agents to maintain a map.
        agents[i].planner->force_map_update(&true_map);
        agents[i].planner->clear_map();
    }
    
    FILE* f = fopen("afghan_map.jl","w");
    save_julia_var(f, std::string("afghan_map"), &true_map);
    fclose(f);
    
    f = fopen("make_video.jl","w");
    init_video_file(f);
    fclose(f);
    
    // Now we run the simulator!
    //For now, all communications work perfectly
    
    TimeVar t_loop;
    TimeVar t_real = timeNow();
    float cluster_times[10000] = {0};
    float compute_times[10000] = {0};
    float assign_times[10000]  = {0};
    float num_targets[10000]     = {0};
    
    for(int iteration = 0; iteration < iters; iteration++){
#ifdef SIM_PLOT
        if(iteration%plot_delay==0)
            system("clear");
#endif
        t_loop = timeNow();
        float realtime =(float)duration(timeNow()-t_real)/float(1000000000);
        printf("Real time: %f s.  Sim time = %f s   (%f x)\n", realtime, float(iteration)*delta_time, (float(iteration)*delta_time)/realtime);
        
        Point true_location_buffer[num_agents];
        Point agent_targets_buffer[num_agents];
        
        // Communications handled here:
        for(int i = 0; i < num_agents; i++){
            for(int j = 0; j < num_agents; j++){
                
                if(i==j){
                    comms[i][j] = true;
                }else{
                    if(collision(agents[j].location, agents[i].location, &true_map)){
                        // P(communicate given wall) = 0.01
                        if(rand()%100 > 1){
                            comms[i][j] = false;
                        }else{
                            comms[i][j] = true;
                        }
                    }else{
                        // P(communicate given no wall) = 0.96
                        if(rand()%100 > 4){
                            comms[i][j] = true;
                        }else{
                            comms[i][j] = false;
                        }
                    }
                }
                if(comms[i][j]){
                    // Assume all communications work fine & immediate
                    agent_location_buffer[i][j].x = agents[j].location.x;
                    agent_location_buffer[i][j].y = agents[j].location.y;
                    agent_location_buffer[i][j].z = agents[j].location.z;
                }
            }
        }
        // Agents need to get their measurements first.
        for(int i = 0; i < num_agents; i++){
            agents[i].planner->update_state(agent_location_buffer[i], num_agents);
            
        }
        // Raytrace to get them measurements
        for(int i = 0; i < num_agents; i++){
            raytrace(i);
        }
        
        // Compute frontiers
        bool done= true;
        TimeVar t_0;
        for(int i = 0; i < num_agents; i++){
            t_0 = timeNow();
            done = agents[i].planner->cluster_frontiers();
            cluster_times[iteration] += float(duration(timeNow()-t_0)/num_agents);
            t_0 = timeNow();
            num_targets[iteration]+=(float)agents[i].planner->compute_costs()/float(num_agents);
            compute_times[iteration] += float(duration(timeNow()-t_0)/num_agents);
            agent_num_clusters[i] = agents[i].planner->get_clusters(agent_clusters_buffer[i], agent_costs_buffer[i], MAX_NUM_GOALS);
        }
        
        // Share costs
        for(int i = 0; i < num_agents; i++){
            for(int j= 0; j < num_agents; j++){
                if(comms[i][j]){
                    agents[i].planner->get_neighbor_clusters(j, agent_clusters_buffer[j], agent_costs_buffer[j], agent_num_clusters[i]);
                }
            }
        }

        // Have them assign, returning target waypoint
        for(int i = 0; i < num_agents; i++){
            t_0 = timeNow();
            agents[i].planner->assign(&agents[i].target);
            assign_times[iteration]+=float(duration(timeNow()-t_0)/num_agents);
        }
#ifdef SIM_DEBUG
        printf("Saving datafiles\n");
#endif
        
        FILE* datafile = fopen("timing.jl", "w");
        if(f == nullptr){
            printf("Error: File pointer is bad\n");
        }else{
            save_julia_var(datafile, std::string("avg_cluster_times"), cluster_times, iteration);
            save_julia_var(datafile, std::string("avg_compute_times"), compute_times, iteration);
            save_julia_var(datafile, std::string("avg_assign_times"), assign_times, iteration);
            save_julia_var(datafile, std::string("num_targets"), num_targets, iteration);
        }
        fclose(datafile);
        
#ifdef SIM_DEBUG
        printf("Moving agents\n");
        for(int i = 0; i < num_agents; i++){
            printf("Agent %d going to (%f,%f,%f)\n", i, agents[i].target.x,agents[i].target.y,agents[i].target.z);
        }
#endif
        
        
        // Here we save the data to the datafile:
        char filename[50];
        for(int i = 0; i < num_agents; i++){
            true_location_buffer[i].x = agents[i].location.x;
            true_location_buffer[i].y = agents[i].location.y;
            true_location_buffer[i].z = agents[i].location.z;
            agent_targets_buffer[i].x = agents[i].target.x;
            agent_targets_buffer[i].y = agents[i].target.y;
            agent_targets_buffer[i].z = agents[i].target.z;
        }
#ifdef SIM_DEBUG
        printf("Saving frame data (%d)\n",iteration);
#endif
        sprintf(filename, "frame%d.jl",iteration);
        datafile = fopen(filename, "w");
        if(datafile==nullptr){
            printf("Error: opened null file\n");
        }else{
            save_julia_var(datafile, std::string("locations"),true_location_buffer , num_agents);
            save_julia_var(datafile, std::string("targets"),agent_targets_buffer, num_agents);
            
            for(int i = 0; i < num_agents; i++){
                sprintf(filename, "trajectory_%d", i);
                save_julia_var(datafile, std::string(filename), agents[i].planner->current_plan, int(agents[i].planner->current_plan[0].x)+1);
            }
            if(iteration%plot_delay == 0){
                save_julia_var(datafile, std::string("global_map"), &global_map);
            }
        }
        fclose(datafile);

        sprintf(filename, "frame%d",iteration);
        datafile = fopen("make_video.jl","a");
        if(datafile==nullptr){
            
        }else{
            print_video_file(datafile, filename, (float)iteration*delta_time);
        }
        fclose(datafile);
        
        // Send agents toward their target destination:
        for(int i = 0; i < num_agents; i++)
            move_base(i);
        
        if(done){
            printf("Exploration complete!\n");
            for(int i = 0;i < num_agents; i++){
                agents[i].planner->print_call_rates();
            }
            return;
        }

        
#ifdef SIM_PLOT
        bool ontime = duration(timeNow()-t_loop) < num_agents*delta_time*(1000000000);
        if(iteration%plot_delay == 0){
            agents[0].planner->plot_clusters();
            }
            
            if(!ontime){
                printf(ANSI_COLOR_YELLOW "Status: Late" ANSI_COLOR_RESET);
            }else{
                printf(ANSI_COLOR_GREEN "Status: Good" ANSI_COLOR_RESET);
            }
            
            while(duration(timeNow()-t_loop) < delta_time*(1000000000)){}
#endif
    }
#ifdef SIM_PLOT
    agents[0].planner->plot_clusters();
#endif
//    printf("Telling julia to make the video\n");
//  system("julia make_video.jl");
    
    
    delete[] true_map.data;
    delete[] global_map.data;
    for(int i = 0; i < num_agents; i++){
        delete agents[i].planner;
    }
}


void Simulator::move_base(int agent){
    // Try using smoothed paths:
    // First, scale time:
    
#ifdef SIM_DEBUG
    Point old_loc;
    old_loc.x = agents[agent].location.x;
    old_loc.y = agents[agent].location.y;
    old_loc.z = agents[agent].location.z;
#endif
    
    float speed_limit = 2*delta_time;
    float grid_vel = dist(agents[agent].location, agents[agent].target);
    if(grid_vel == 0){
        return;
    }
    float dx = (agents[agent].target.x - agents[agent].location.x)/grid_vel;
    float dy = (agents[agent].target.y - agents[agent].location.y)/grid_vel;
    float dz = (agents[agent].target.z - agents[agent].location.z)/grid_vel;
    
    // Scale to be proper speed for sim:
    if(grid_vel*5.0/36.0 > speed_limit)
        grid_vel = (36.0/5.0)*speed_limit;
    
    agents[agent].bearing = atan2f(dy, dx);
    
    agents[agent].location.x += dx*grid_vel;
    agents[agent].location.y += dy*grid_vel;
    agents[agent].location.z += dz*grid_vel;
#ifdef SIM_DEBUG
    printf("Deltas are %f %f %f. Grid velocity is %f ", dx,dy,dz,grid_vel);
    printf("Moved agent %f distance\n", dist(agents[agent].location, old_loc));
#endif
    
}


// Traces a ray from the agents position and pushes the measurement to the agent's planner
// This is wasteful, but just a simulator so it doesn't matter much.
void Simulator::raytrace(int agent){
    // First we run the ray:
    float range = 60;
    agents[agent].bearing = nice_angle(agents[agent].bearing);
    float d_theta = 0.01;
    Point lower;
    Point upper;
    lower.x = true_map.X_dim;
    lower.y = true_map.Y_dim;
    lower.z = 0;
    upper.x = 0;
    upper.y = 0;
    upper.z = 0;
    
    Map measurement;
    measurement.X_dim = true_map.X_dim;
    measurement.Y_dim = true_map.Y_dim;
    measurement.Z_dim = true_map.Z_dim;
    measurement.data  = new char[measurement.X_dim*measurement.Y_dim];
    for(int i = 0; i < measurement.X_dim*measurement.Y_dim; i++)
        measurement.data[i] = MAP_UNKN;
    
    // Special case di=dj=1:
    if(nice_angle(agents[agent].bearing-M_PI_4) <= 3*M_PI_4 || nice_angle(agents[agent].bearing-M_PI_4) >= 5*M_PI_4){
        for(int n = 0; float(n)*M_SQRT2 < range; n++){
            Point pos;
            pos.x = agents[agent].location.x + float(n)*M_SQRT1_2;
            pos.y = agents[agent].location.y + float(n)*M_SQRT1_2;
            pos.z = 0;

            
            if(pos.x < true_map.X_dim && pos.y < true_map.Y_dim && pos.x >=0 && pos.y >= 0){
                
                if(pos.x < lower.x)
                    lower.x = pos.x;
                if(pos.y < lower.y)
                    lower.y = pos.y;
                if(pos.x > upper.x)
                    upper.x = pos.x;
                if(pos.y > upper.y)
                    upper.y = pos.y;
                
                char meas = true_map(int(pos.x), int(pos.y), int(pos.z));
                measurement.set(int(pos.x), int(pos.y),int(pos.z), true_map(int(pos.x), int(pos.y), int(pos.z)));
                global_map.set(int(pos.x), int(pos.y), int(pos.z), true_map(int(pos.x), int(pos.y), int(pos.z)));
                if(meas >= MAP_OCC_THRESH){
                    break;
                }
            }
        }
    }
    
    for(int di = -1; di <= 1; di+=2){
        for(float dj = -1; dj <= 1; dj+=d_theta){
            float theta = nice_angle(agents[agent].bearing - atan2(float(dj), float(di)));
            if(theta <= 3*M_PI_4 || theta >= 5*M_PI_4){
                for(int n = 0; (sqrtf((n*di)*(n*di) + (n*dj)*(n*dj)) < range); n++){
                    Point pos;
                    pos.x = agents[agent].location.x + n*di;
                    pos.y = agents[agent].location.y + n*dj;
                    pos.z = 0;
                    
                    if(pos.x < true_map.X_dim && pos.y < true_map.Y_dim && pos.x >=0 && pos.y >= 0){
                        
                        if(pos.x < lower.x)
                            lower.x = pos.x;
                        if(pos.y < lower.y)
                            lower.y = pos.y;
                        if(pos.x > upper.x)
                            upper.x = pos.x;
                        if(pos.y > upper.y)
                            upper.y = pos.y;
                        
                        char meas = true_map(int(pos.x), int(pos.y), int(pos.z));
                        measurement.set(int(pos.x), int(pos.y),int(pos.z), true_map(int(pos.x), int(pos.y), int(pos.z)));
                        global_map.set(int(pos.x), int(pos.y), int(pos.z), true_map(int(pos.x), int(pos.y), int(pos.z)));
                        if(meas >= MAP_OCC_THRESH){
                            break;
                        }
                        
                    }
                }
            }
        }
    }
    
    for(int dj = -1; dj <= 1; dj+=2){
        for(float di = -1; di <= 1; di+=d_theta){
            float theta = nice_angle(agents[agent].bearing - atan2(float(dj), float(di)));
            if(theta <= 3*M_PI_4 || theta >= 5*M_PI_4){
                for(int n = 0; (sqrtf((n*di)*(n*di) + (n*dj)*(n*dj)) < range); n++){
                    Point pos;
                    pos.x = agents[agent].location.x + n*di;
                    pos.y = agents[agent].location.y + n*dj;
                    pos.z = 0;
                    if(pos.x < true_map.X_dim && pos.y < true_map.Y_dim && pos.x >=0 && pos.y >= 0){
                        
                        if(pos.x < lower.x)
                            lower.x = pos.x;
                        if(pos.y < lower.y)
                            lower.y = pos.y;
                        if(pos.x > upper.x)
                            upper.x = pos.x;
                        if(pos.y > upper.y)
                            upper.y = pos.y;
                        
                        char meas = true_map(int(pos.x), int(pos.y), int(pos.z));
                        measurement.set(int(pos.x), int(pos.y),int(pos.z), true_map(int(pos.x), int(pos.y), int(pos.z)));
                        global_map.set(int(pos.x), int(pos.y), int(pos.z), true_map(int(pos.x), int(pos.y), int(pos.z)));

                        if(meas >= MAP_OCC_THRESH){
                            break;
                        }
                        
                    }
                }
            }
        }
    }
    
    // Now send the measurement
    agents[agent].planner->update_map(&measurement, lower, upper);
    for(int j = num_agents-1; j >=0; j--){
        if(comms[agent][j]){
            // Agent can hear j
            agents[agent].planner->merge_map(agents[j].planner->get_map());
        }
    }
    delete [] measurement.data;    
}




