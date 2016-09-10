//
//  main.cpp
//  FMT
//
//  Created by Megamind on 7/18/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#include "fmt.h"

// Explicit constructor method
FMT::FMT(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param){
    // set parameters
    set_parameters(xlimit, ylimit, zlimit, num_pts, connection_type, connection_param);
    max_neighborhood_size = int(num_pts/5);
    
    // Preallocate vectors (+2 because we add start and end points at the end)
    points = new Point[(parameters.num_pts+2)];
    sample_points();
    
    // Compute neighborhoods (return if error)
    if(!compute_neighborhoods())
        return;
#ifdef FMT_DEBUG
    print_parameters();
    printf("compute neighborhoods successful\n");
#endif
    
    
    // Initialize planning containers:
    Windex = new int[parameters.num_pts+1];
    Parents = new int[parameters.num_pts+2];
    Costs = new float[parameters.num_pts+2];
    Active = new int[parameters.num_pts+2];
    H_new = new int[parameters.num_pts+2];
    x_near.costs = new float[parameters.num_pts+2];
    x_near.indices = new int[parameters.num_pts+2];
    
    y_near.costs = new float[parameters.num_pts+2];
    y_near.indices = new int[parameters.num_pts+2];

    num_skipped_checks=0;
    num_total_checks=0;
    
    //TODO: Save precomputed neighborhoods to a parameter file
    is_initialized = true;
}

FMT::~FMT(){
    printf("Destructing FMT object\n");

    if(is_initialized){
        // Erase neighborhoods:
        for(int i = 0; i < parameters.num_pts+2; i++){
            delete [] neighborhoods[i].indices;
            delete [] neighborhoods[i].costs;
            delete [] neighborhoods[i].cache;
        }
        delete [] neighborhoods;
        
        // Erase planning containers
        delete [] Windex;
        delete [] Parents;
        delete [] Costs;
        delete [] Active;
        delete [] H_new;
        delete [] x_near.costs;
        delete [] x_near.indices;
        delete [] y_near.costs;
        delete [] y_near.indices;
        
    }
    delete [] points;
}

// TODO: make a destructor to handle memory issues.

void FMT::set_parameters(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param){
    parameters.X_limit = xlimit;
    parameters.Y_limit = ylimit;
    parameters.Z_limit = zlimit;
    parameters.num_pts = num_pts;
    parameters.connection_type = connection_type;
    parameters.connection_param = connection_param;
}


/*** Planning methods ***/

// Points used for heap in FMT*
typedef struct Heap_pts {
    int index;
    float cost;
}Heap_pt;

// Reversed operator for minheap
bool operator<(const Heap_pt& lhs, const Heap_pt& rhs){
    return (lhs.cost > rhs.cost);
}



// Finds (likely) shortest path from start point to point in goal cluster
void FMT::fmtstar(Point start, Cluster *goal, Map *map, Point* path){
    
    if(isnan(start.x) || isnan(start.y) || isnan(start.z) || isnan(goal->median.x) || isnan(goal->median.y) || isnan(goal->median.z)){
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        path[0].z = 0;
        return;
    }
    
    //push start, goal onto points list
    points[start_pt_ind].x = start.x;
    points[start_pt_ind].y = start.y;
    points[start_pt_ind].z = start.z;
    
    // Mark goal as median point:
    points[goal_pt_ind].x  = goal->median.x;
    points[goal_pt_ind].y  = goal->median.y;
    points[goal_pt_ind].z  = goal->median.z;
    
    // Check for trivial path
    // Find the closest point:
    Point pt;
    for(int ind = 0; ind < goal->size; ind++){
        map->num2ind(&pt, goal->members[ind]);
        if(dist(start, pt) < path[0].y){
            path[0].x = 2;
            path[0].y = dist(start, pt);
            path[0].z = 0;
            path[1].x = start.x;
            path[1].y = start.y;
            path[1].z = start.z;
            path[2].x = pt.x;
            path[2].y = pt.y;
            path[2].z = pt.z;
        }
    }
    
    if(!collision(path[1], path[2], map)){
        // Path is easy!
        return;
    }

    // Start by constructing neighborhoods for start/finish points
    if(! compute_neighborhood(start, start_pt_ind, map) ||! compute_neighborhood(goal, goal_pt_ind,map)){
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        path[0].z = 0;
        return;
    }
    return plan_path(map, path);
}


// Finds shortest path between members or centroids of start and goal clusters.
// Note - centroid behavior might be different than expected.
void FMT::fmtstar(Cluster *start, Cluster *goal, Map *map, Point* path){
    
    if(isnan(start->median.x) || isnan(start->median.y) || isnan(start->median.z) || isnan(goal->median.x) || isnan(goal->median.y) || isnan(goal->median.z)){
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        path[0].z = 0;
        return;
    }
    
    //push start, goal onto points list
    points[start_pt_ind].x = start->median.x;
    points[start_pt_ind].y = start->median.y;
    points[start_pt_ind].z = start->median.z;
    
    points[goal_pt_ind].x  = goal->median.x;
    points[goal_pt_ind].y  = goal->median.y;
    points[goal_pt_ind].z  = goal->median.z;

    // check closest point for line of sight:
    Point pt_g, pt_s;
    path[0].y = MAXFLOAT;
    for(int ind_g = 0; ind_g < goal->size; ind_g++){
        map->num2ind(&pt_g, goal->members[ind_g]);
        for(int ind_s=0; ind_s < start->size; ind_s++){
            map->num2ind(&pt_s, start->members[ind_s]);
            if(dist(pt_s, pt_g) < path[0].y){
                path[0].x = 2;
                path[0].y = dist(pt_s, pt_g);
                path[0].z = 0;
                path[1].x = pt_s.x;
                path[1].y = pt_s.y;
                path[1].z = pt_s.z;
                path[2].x = pt_g.x;
                path[2].y = pt_g.y;
                path[2].z = pt_g.z;
            }
        }
    }
    if(!collision(path[1], path[2], map)){
        // Path is easy!
        return;
    }
    
    // Start by constructing neighborhoods for start/finish points
    if(! compute_neighborhood(start, start_pt_ind,map) ||! compute_neighborhood(goal, goal_pt_ind,map)){
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        path[0].z = 0;
    }
    return plan_path(map, path);
}

void FMT::fmtstar(Point start, Point goal, Map* map, Point* path){
#ifdef FMT_TIMING
    TimeVar t1 = timeNow();
    double initialization_time = 0;
#endif
    
    //push start, goal onto points list
    points[start_pt_ind].x = start.x;
    points[start_pt_ind].y = start.y;
    points[start_pt_ind].z = start.z;
    
    points[goal_pt_ind].x  = goal.x;
    points[goal_pt_ind].y  = goal.y;
    points[goal_pt_ind].z  = goal.z;

 // First, check whether start and goal are trivially connected (check line of sight)
    if(isnan(start.x) || isnan(start.y) || isnan(start.z) || isnan(goal.x) || isnan(goal.y) || isnan(goal.z)){
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        path[0].z = 0;
        return;
    }
    
    
    
     if( !collision(start, goal, map)){
#ifdef FMT_DEBUG
         printf("Easy path!");
#endif
         path[0].x = 2;
         path[0].y = dist(start,goal);
         path[0].z = 0;
         path[1].x = start.x;
         path[1].y = start.y;
         path[1].z = start.z;
         path[2].x = goal.x;
         path[2].y = goal.y;
         path[2].z = goal.z;
         return;
     }
    

    // Start by constructing neighborhoods for start/finish points
    if(! compute_neighborhood(start, start_pt_ind, map) ||! compute_neighborhood(goal, goal_pt_ind, map)){
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        path[0].z = 0;
        return;
    }
    return plan_path(map, path);
}

void FMT::plan_path(Map* map, Point* path){
    
    // Initialize containers for search
    
    int Windex_size = parameters.num_pts+1;
    for(int i = 0; i < parameters.num_pts+1; i++)
        Windex[i] = i;
    
    int z = start_pt_ind;
    
    int Active_size = 1;
    Active[0] = z;
    
    for(int i = 0; i < parameters.num_pts+2; i++){
        Costs[i] = 0.0; //Cost to get somewhere should be very high unless we've found a path to it
        Parents[i] = -1;
    }
    Costs[start_pt_ind] = 0.0;

    // Containers are [index cost] tuples where cost is the shortest path found from start to points[index]
    std::priority_queue<Heap_pt> HHeap;
    
    int H_new_size = 0;
    
    // Form two empty neighborhoods for use in the loop:

    // Reset planning neighborhoods:
    x_near.size = 0;
    y_near.size = 0;
    
    
    while(! is_goal_pt(z)){

#ifdef FMT_DEBUG
//        printf("Removing %d from Windex\n", z);
#endif

        Windex_size = remove_element(Windex, Windex_size, z);

        // Reset containers:
        H_new_size = 0;
        
        //Construct neighborhood of z
        filter(&x_near, z, Windex, Windex_size, -1);
        
        int ind_x = 0;
        for(int i_x = 0; i_x < x_near.size; i_x++){
            ind_x = x_near.indices[i_x];
            filter(&y_near, ind_x, Active, Active_size, ind_x);
            
            if(y_near.size > 0){
                // Find minimum parent:
                float cmin = MAXFLOAT;
                int y_min = -1;
                for(int i_y = 0; i_y < y_near.size; i_y++){
                    if(ind_x == y_near.indices[i_y]) // Don't connect to self!
                        continue;
                    
                    float c = Costs[y_near.indices[i_y]]+y_near.costs[i_y];
                    if(c < cmin){
                        cmin  = c;
                        y_min = i_y;
                    }
                }
                
                if(y_min != -1 && !collision_inds(ind_x, y_near.indices[y_min], map)){

                    // Add parent
                    Parents[ind_x] = y_near.indices[y_min];
                    
                    // Update cost to go
                    Costs[ind_x] = cmin;
                    
                    // Mark as visited
                    H_new[H_new_size] = ind_x;
                    H_new_size += 1;

                    Windex_size = remove_element(Windex, Windex_size, ind_x);

                    // Add to heap
                    Heap_pt new_pt;
                    new_pt.index = ind_x;
                    new_pt.cost  = cmin;
                    HHeap.push(new_pt); // TODO: CHECK THIS!! Is it copied in, or passed by reference?
                }else{
                    //
                }
            }

        }

        // Update active list
        Active_size = join_sets(Active, Active_size, H_new, H_new_size);

        // Remove z
        Active_size = remove_element(Active,Active_size, z);

        if(HHeap.size() > 0){
            z = HHeap.top().index;
            HHeap.pop();
            
#ifdef FMT_DEBUG
//            printf("popping %d with parent %d and %d neighbors\n", z, Parents[z], neighborhoods[z].size);
#endif
        }else{
            // Try to connect to goal in case this is a sparse sampling issue
            if(!collision(points[z],points[goal_pt_ind], map)){
#ifdef FMT_DEBUG
                printf("Heap ran empty - trying to connect directly to goal\n");
#endif
                Parents[goal_pt_ind] = z;
                Costs[goal_pt_ind] = Costs[z] + dist(points[z], points[goal_pt_ind]);
                z = goal_pt_ind;
            }else{

#ifdef FMT_DEBUG
                std::cout << "Heap is empty - failure.\n";
#endif
                path[0].x = -1;
                path[0].y = MAXFLOAT;
                path[0].z = 0;
                return;
            }
            break;
        }
    }

    // Construct path by assembling parents from goal to start:
    // Going to end up backwards this way
    int path_length = 0;
    if(z!= goal_pt_ind){
        path[parameters.num_pts+1] = points[goal_pt_ind];
        path[parameters.num_pts] = points[z];
        path_length = 2;
    }else{
        path[parameters.num_pts+1] = points[z];
        path_length = 1;
    }

    // First point holds informtion about the path:
    // -1 => invalid path
    // x holds the number of elements in the path
    // y holds the cost of the path
    // z is empty
    path[0].x = -1;
    path[0].y = MAXFLOAT;
    path[0].z =0.0;

    // Insert points backwards from the back
    int curr_ind = z;
    while(Parents[curr_ind] >= 0){
        // Error checking:
        if(path_length >= parameters.num_pts+2) // Error in path - too many elements
        {path_length = -1; break;}
        if(curr_ind == Parents[curr_ind]) // Error in path - loop
        {path_length = -1; break;}
        
        //Actual code
        path[(parameters.num_pts+1) - path_length] = points[Parents[curr_ind]];
        curr_ind = Parents[curr_ind];
        path_length++;
    }
    
    if(path_length >= 0){
        // The path does not include start_pt yet. Add it then smooth
        path[1].x = points[start_pt_ind].x; path[1].y = points[start_pt_ind].y; path[1].z = points[start_pt_ind].z;
        path[0].x = 1;
        path[0].y = 0.0;

        // Index of raw path
        int raw_ind = parameters.num_pts+2-path_length;
        // Index of smoothed path
        int smooth_ind = 1;
        // Step through raw path and skip points if unnecessary
        while(raw_ind < parameters.num_pts+2){
            // Means we need the previous point
            if(collision(path[raw_ind], path[smooth_ind], map)){
                
                // Slide the point toward the following one as far as you can:
                float num_steps = 10;
                float dx = (path[raw_ind].x- path[raw_ind-1].x)/num_steps;
                float dy = (path[raw_ind].y- path[raw_ind-1].y)/num_steps;
                float dz = (path[raw_ind].z- path[raw_ind-1].z)/num_steps;
                Point new_pt;
                new_pt.x = path[raw_ind-1].x;
                new_pt.y = path[raw_ind-1].y;
                new_pt.z = path[raw_ind-1].z;
                
                for(int k = 0; k < num_steps; k++){
                    new_pt.x += dx;
                    new_pt.y += dy;
                    new_pt.z += dz;
                    if(collision(new_pt, path[smooth_ind], map)){
                        break;
                    }
                    // Can inch the path forward without collision
                    path[raw_ind-1].x+=dx;
                    path[raw_ind-1].y+=dy;
                    path[raw_ind-1].z+=dz;
                }
                
                smooth_ind++;
                path[smooth_ind].x = path[raw_ind-1].x;
                path[smooth_ind].y = path[raw_ind-1].y;
                path[smooth_ind].z = path[raw_ind-1].z;
                path[0].x++;
                path[0].y += dist(path[smooth_ind-1], path[smooth_ind]);
            }
            raw_ind++;
            if(raw_ind == parameters.num_pts+2){
                smooth_ind++;
                path[smooth_ind].x = path[raw_ind-1].x;
                path[smooth_ind].y = path[raw_ind-1].y;
                path[smooth_ind].z = path[raw_ind-1].z;
                path[0].x++;
                path[0].y += dist(path[smooth_ind-1], path[smooth_ind]);
            }
        }
    }else{
        // Path generation failed - return badpath.
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        path[0].z = 0;
        return;
    }

#ifdef FMT_DEBUG
    std::cout << "Finished at goal point (" << points[z].x << ", " << points[z].y << ", " << points[z].z << ") with cost " << Costs[z]<< "\n";
    std::cout<< "Percent of sites explored: " << float(parameters.num_pts - Windex_size)/parameters.num_pts << "\n";
    printf("Path has %d points. Status: %f,%f\n", path_length, path[0].x, path[0].y);
    for(int i = 0; i < path_length; i++){
        printf("(%f,%f)->", path[i+1].x,path[i+1].y);
    }
    printf("\n");
#endif
    
    
    // Good path! Smooth and check:
    bool path_has_collision = true;
    Point pp[100];
    pp[0].x = path[0].x;
    pp[0].y = path[0].y;
    int num_smooth_points = pp[0].x;
    if(4 < pp[0].x)
        num_smooth_points = 4;
    
    for(int k = 0; k < num_smooth_points; k++){ // Only bother smoothing first few points
        pp[k+1].x = path[k+1].x;
        pp[k+1].y = path[k+1].y;
        pp[k+1].z = path[k+1].z;
    }

    Point v_init;
    Point a_init;
    Point v_final;
    Point a_final;
    // V_init in direction of first segment:
    v_init.x = path[2].x - path[1].x;
    v_init.y = path[2].y - path[1].y;
    v_init.z = path[2].z - path[1].z;
    float len = norm(v_init);
    if(len != 0.0){
        v_init.x /= len;
        v_init.y /= len;
        v_init.z /= len;
    }
    a_init.x = 0; a_init.y = 0; a_init.z = 0;
    // V_final pointing toward last segment:
    v_final.x = 0; v_final.y = 0; v_final.z = 0;
    if(pp[0].x > num_smooth_points){
        v_final.x = path[num_smooth_points+1].x - path[num_smooth_points].x;
        v_final.y = path[num_smooth_points+1].y - path[num_smooth_points].y;
        v_final.z = path[num_smooth_points+1].z - path[num_smooth_points].z;
        len = norm(v_init);
        if(len != 0.0){
            v_final.x /= len;
            v_final.y /= len;
            v_final.z /= len;
        }
    }
    a_final.x = 0; a_final.y = 0; a_final.z = 0;
    
    
    
    Smoothed = PathSmoother();
    int fail_counter = 0;
    Point pos;
    while(path_has_collision){
        path_has_collision = false;
        // Smooth path
        Smoothed.smooth_path(pp, &v_init, &a_init, &v_final, &a_final);

        // check for collision
        // Dumb way: check at 100 points
        int num_sample = 100;
        int point_index = 1; // Index of the next point in the list
        
        float dt = Smoothed.T_split[(int)pp[0].x-1]/(float)num_sample;
        float t_eval = 0;
        for(int t = 0; t < num_sample; t++){
            t_eval += dt;
            Smoothed.get_point(t_eval, &pos);
            // mega fail.
            if(pos.x != pos.x || pos.y != pos.y){
                printf("t_eval = %f\n", t_eval);
                break;
            }

            // Check if close to next point
            if( dist(pos, pp[point_index+1]) < 0.1){
                point_index++;
            }
            
            // Check at point
            if(!map->is_free(pos)){
#ifdef FMT_DEBUG
                printf("Path has collision at (%f,%f) (seg %d)!\n", pos.x,pos.y, point_index);
#endif
                path_has_collision = true;
                break;
            }
        }
        
        // If collision, add point
        if(path_has_collision){
            // Shift points to make room
            pp[0].x += 1;
            for(int k = (int)pp[0].x; k > point_index; k--){
                pp[k+1].x = pp[k].x;
                pp[k+1].y = pp[k].y;
                pp[k+1].z = pp[k].z;
            }
            
            // Add a point at the projection of the collision onto the segment AB:
            //
            //
            // C            B
            //   .        *
            //     .    *
            //        o
            //      *
            //    *
            //  *
            //A
            //
            // The intersection is solved by doing a linear projection:
            //
            // (AC)•(AB)/(|AB|) * (AB)
            
            float dx_ab = pp[point_index+2].x - pp[point_index].x;
            float dy_ab = pp[point_index+2].y - pp[point_index].y;
            float dx_ac = pos.x - pp[point_index].x;
            float dy_ac = pos.y - pp[point_index].y;
            float coeff = (dx_ac*dx_ab + dy_ac*dy_ab)/(dx_ab*dx_ab + dy_ab*dy_ab);
            if(coeff <= 1.0 || coeff < 0.0){
                pp[point_index+1].x = coeff*dx_ab + pp[point_index].x;
                pp[point_index+1].y = coeff*dy_ab + pp[point_index].y;
            }else{
                pp[point_index+1].x = pp[point_index].x+dx_ab/2.0;
                pp[point_index+1].y = pp[point_index].y+dy_ab/2.0;
            }
            // Now we move the old point back as much as we can:
            if(point_index != 1){
                float num_steps = 30;
                dx_ab = (pp[point_index].x - pp[point_index-1].x)/num_steps;
                dy_ab = (pp[point_index].y - pp[point_index-1].y)/num_steps;
                for(int k = 0; k < num_steps; k++){
                    pp[point_index].x -= dx_ab;
                    pp[point_index].y -= dy_ab;
                    
                    if(collision(pp[point_index], pp[point_index+1], map)){
                        pp[point_index].x += dx_ab;
                        pp[point_index].y += dy_ab;
                        break;
                    }
                }
            }
            
            pp[point_index+1].z = (pp[point_index+2].z + pp[point_index].z)/2;
            fail_counter++;
        }
        if(pp[0].x > 12){
            break;
        }
    }
#ifdef FMT_JULIA_DEBUG
    printf("chomp(readline()); clf();\n");
    printf("println(\"Drawing new figure from point (%f,%f)\");\n", path[1].x, path[1].y);

    // Show original path
    printf("fmt_path = [");
    for(int k = 0; k < path[0].x; k++){
        printf("%f %f %f ", path[k+1].x, path[k+1].y, path[k+1].z);
        if(k == int(path[0].x)-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }
    
#endif
    // Replace path with sampled polynomial
    float t = 0;
    Point smoothed_vel[100];
    Point smoothed_acc[100];
    Point smoothed_jerk[100];
    float times[100];
    float dt = 1.0;
    
    for(path[0].x = 1; path[0].x < 100; path[0].x++){
        int k = (int)path[0].x;
        if(&path[k] == NULL)
            break;
        Smoothed.get_point(t,&(path[k]));
        Smoothed.get_der(t, 1, &(smoothed_vel[k-1]));
        Smoothed.get_der(t, 2, &(smoothed_acc[k-1]));
        Smoothed.get_der(t, 3, &(smoothed_jerk[k-1]));
        times[k-1] = t;
        t+=dt;
        if(t >= Smoothed.T_split[Smoothed.num_points-1])
            break;
        
    }
#ifdef FMT_JULIA_DEBUG
    // Show original path
    printf("smoothed = [");
    for(int k = 0; k < path[0].x; k++){
        printf("%f %f %f ", path[k+1].x, path[k+1].y, path[k+1].z);
        if(k == int(path[0].x)-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }
    
    printf("smoothed_vel = [");
    for(int k = 0; k < path[0].x; k++){
        printf("%f %f %f ", smoothed_vel[k].x, smoothed_vel[k].y, smoothed_vel[k].z);
        if(k == int(path[0].x)-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }

    printf("smoothed_acc = [");
    for(int k = 0; k < path[0].x; k++){
        printf("%f %f %f ", smoothed_acc[k].x, smoothed_acc[k].y, smoothed_acc[k].z);
        if(k == int(path[0].x)-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }
    printf("smoothed_jerk = [");
    for(int k = 0; k < path[0].x; k++){
        printf("%f %f %f ", smoothed_jerk[k].x, smoothed_jerk[k].y, smoothed_jerk[k].z);
        if(k == int(path[0].x)-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }
    printf("times = [");
    for(int k = 0; k < path[0].x; k++) {
        printf("%f",times[k]);
        if(k==int(path[0].x)-1){
            printf("];\n");
        }else{
            printf(";");
        }
    }
    
    // Show new points
    printf("fmt_path_pp = [");
    for(int k = 0; k < int(pp[0].x); k++){
        printf("%f %f %f ", pp[k+1].x, pp[k+1].y, pp[k+1].z);
        if(k == int(pp[0].x)-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }

    // Show smoothed path
    printf("coeff_x = [");
    for(int k = 0; k < Smoothed.degree; k++){
        printf("%f", Smoothed.x_coeffs[k]);
        if(k == Smoothed.degree-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }
    printf("coeff_y = [");
    for(int k = 0; k < Smoothed.degree; k++){
        printf("%f", Smoothed.y_coeffs[k]);
        if(k == Smoothed.degree-1){
            printf("];\n");
        }else{
            printf("; ");
        }
    }
//    printf("clf();\n");
//    printf("plot_poly(coeff_x,coeff_y,fmt_path_pp, %f, 3);\n", Smoothed.T_split[Smoothed.num_points-1]);
    printf("subplot(2,2,1);\n");
    printf("imshow(afghan_map', cmap=\"gray\", interpolation=\"none\");\n");
    printf("plot(fmt_path[:,1],fmt_path[:,2])\n");
    printf("scatter(%f,%f,marker=\"+\", color=:orange);\n", pos.x, pos.y);
    printf("scatter(fmt_path_pp[:,1],fmt_path_pp[:,2],marker=\"o\",color=:green);\n");
    printf("scatter(fmt_path[:,1], fmt_path[:,2],marker=\"x\",color=:red);\n");
    printf("plot(smoothed[:,1], smoothed[:,2],color=:red);\n");
    printf("axis([0,200,0,200]);\n");
    printf("subplot(2,2,2)\n");
    printf("plot(times, smoothed_vel[:,1], color=:red);\n");
    printf("plot(times, smoothed_vel[:,2], color=:blue);\n");
    printf("subplot(2,2,3)\n");
    printf("plot(times, smoothed_acc[:,1], color=:red)\n");
    printf("plot(times, smoothed_acc[:,2], color=:blue)\n");
    printf("subplot(2,2,4)\n");
    printf("plot(times, smoothed_jerk[:,1], color=:red)\n");
    printf("plot(times, smoothed_jerk[:,2], color=:blue)\n");

#endif 
    
    // Clean up planner containers
    reset_neighborhood(start_pt_ind);
    reset_neighborhood(goal_pt_ind);

    
    if(path_has_collision){
#ifdef FMT_DEBUG
        printf("Smoothing aborted!\n");
#endif
        path[0].x = -1;
        path[0].y = MAXFLOAT;
        return;
    }else{
        path[0].y = Smoothed.get_length();
    }
    
    return;
}

// Checks whether point is close to goal
bool FMT::is_goal_pt(int index){
    // Search neighborhood of goal to see:
    for(int i = 0; i < neighborhoods[goal_pt_ind].size; i++){
        if(neighborhoods[goal_pt_ind].indices[i] == index)
            return true;
    }
    return false;
}


// Returns the neighborhood of index filtered by filter
// Since the filter is sorted, can speed up
// Maintains same sorting as original neighborhood
// Assumes that neighborhoods[index].indices and filter_vec are sorted, smallest first.
void FMT::filter(Neighborhood* filtered_neighborhood, int index, int* filter_vec, int filter_size, int exclude){
    int i_n = 0; int i_f = 0; int ind_n = 0;
    filtered_neighborhood->size = 0;
    float c = 0;
    
    while(i_n < neighborhoods[index].size && i_f < filter_size){
        ind_n = neighborhoods[index].indices[i_n];
        if(ind_n == filter_vec[i_f]){
            if(ind_n != exclude){
                c = neighborhoods[index].costs[i_n];
                filtered_neighborhood->indices[filtered_neighborhood->size] = ind_n;
                filtered_neighborhood->costs[filtered_neighborhood->size]   = c;
                filtered_neighborhood->size+=1;
            }
            i_n+=1; i_f += 1;
        }else if( ind_n < filter_vec[i_f]){
            i_n += 1;
        }else{
            i_f += 1;
        }
    }
}


/*** Private methods ***/

void FMT::sample_points(){
#ifdef FMT_DEBUG
    printf("Sampling points\n");
#endif
    
    goal_pt_ind  = parameters.num_pts;
    start_pt_ind = parameters.num_pts+1;
    
    // Generate random numbers
    float x_scaling = float(parameters.X_limit)/RAND_MAX;
    float y_scaling = float(parameters.Y_limit)/RAND_MAX;
    for(int i = 0; i < parameters.num_pts; i++)
    {
        points[i].x = rand()*x_scaling;
        points[i].y = rand()*y_scaling;
    }
    
    // Only call rand() if we're planning in Z
    if(parameters.Z_limit == 0){
        for(int i = 0; i < parameters.num_pts; i++)
            points[i].z = 0.0;
    }else{
        float z_scaling = parameters.Z_limit/RAND_MAX;
        for(int i = 0; i < parameters.num_pts; i++)
            points[i].z = rand()*z_scaling;
    }
}

// Computes a single neighborhood
bool FMT::compute_neighborhood(Point pt, int index, Map* map){
    
    if(!is_initialized){
        printf("Call to compute_neighborhood before initialization is complete. \n");
        return false;
    }
    
    if(&(neighborhoods[index]) == NULL){
        printf("Call to FMT with uninitialized neighborhood for point %d\n", index);
        return false;
    }
    
    neighborhoods[index].size = 0;
    
    if(parameters.connection_type == RAD_CON){
        // Radially connected neighbors
        float dim = 2.0;
        if(parameters.Z_limit > 0)
            dim = 3.0;
        // Coefficient from equation (3) of FMT* paper (Janson, Schmerling, et al.)
        // Not sure what xi is in their equation (tuning parameter??)
        float coeff = (parameters.X_limit*parameters.Y_limit*log(parameters.num_pts))/(dim*parameters.num_pts*1.1);
        parameters.connection_param = 2.2*powf(coeff, 1/dim);
        // Already set somewhere else
        for(int j = 1; j < parameters.num_pts; j++){
            if(dist(points[index],points[j]) < parameters.connection_param)
            {
                if(!collision(points[index], points[j], map)){
                    // Add to neighbor list

                    if(neighborhoods[index].size < max_neighborhood_size){
                        neighborhoods[index].indices[neighborhoods[index].size] = j;
                        neighborhoods[index].costs[neighborhoods[index].size] = dist(points[index],points[j]);
                        neighborhoods[index].cache[neighborhoods[index].size] = CACHE_FREE;
                        neighborhoods[index].size += 1;
                    }
                    
                    // Since r-connected graphs are bidirectional, add to other neighborhood
                    if(neighborhoods[j].size < max_neighborhood_size){
                        neighborhoods[j].indices[neighborhoods[j].size] = index;
                        neighborhoods[j].costs[neighborhoods[j].size]   = dist(points[index],points[j]);
                        neighborhoods[j].cache[neighborhoods[j].size]   = CACHE_FREE;
                        neighborhoods[j].size += 1;

                    }else{
                    }
                }
            }
        }
        
        // TODO: Sort the neighborhood by distance (? not sure if it will improve performance at all)
        return true;
    }else if(parameters.connection_type == KNN_CON){
        /*
         // K nearest neighbors connected (this is slow to compute)
         int K = int(parameters.connection_param) + 1; // +1 compensates for including self in neighborhood
         // Array with nearest distances:
         float near_dists[K];
         int near_indices[K];
         for(int k = 0; k < K; k++)
         near_dists[k] = MAXFLOAT;
         
         // Compute distances to all points and track the closest K
         float d_ij = 0;
         for(int j = 0; j < parameters.num_pts; j++){
         d_ij = dist(points[index],points[j]);
         // If smaller than current best, insert into best vectors.
         for(int k = 0; k < K; k++){
         if(d_ij < near_dists[k])
         {
         //shift & insert at k
         for(int kk = K; kk > k; kk--){
         near_dists[kk] = near_dists[kk-1];
         near_indices[kk] = near_indices[kk-1];
         }
         near_dists[k] = d_ij; near_indices[k] = j;
         }
         break;
         }
         }
         // Insert K nearest neighbors into neighborhood:
         neighborhoods[index].indices.reserve(K);
         neighborhoods[index].costs.reserve(K);
         for(int k = 0; k < K; k++){
         neighborhoods[index].indices[k] = near_indices[k];
         neighborhoods[index].costs[k]   = near_dists[k];
         }
         
         // Now need to check all nearby neighbors. That's going to suck!
         // TODO: Do this.
         
         return true;
         */
        printf("Error: KNN not implemented\n");
        return false;
    }else{
        printf("Error (1): Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }
}


// Computes a single neighborhood
bool FMT::compute_neighborhood(Cluster* clust, int index, Map* map){
    
    if(!is_initialized){
        printf("Call to compute_neighborhood before initialization is complete. \n");
        return false;
    }
    
    neighborhoods[index].size = 0;
    
    if(parameters.connection_type == RAD_CON){
        int step = int(clust->size/32)+1; // IF the neighborhood is huge, approximate it
        
        for(int j = 1; j < parameters.num_pts; j++){
            bool neighbor = false;
            for(int i = 0; i < clust->size; i+=step){
            
                Point pt;
                map->num2ind(&pt, clust->members[i]);
            
                // Don't add if it's no better than what's already known
                if(neighbor && (dist(pt, points[j])> neighborhoods[j].costs[neighborhoods[j].size-1]))
                   continue;
                if(dist(pt,points[j]) < parameters.connection_param)
                {
                    if(!collision(pt, points[j], map)){
                        if(neighbor){
                            //Forget about the last element in neighborhoods[j]
                            neighborhoods[j].size--;
                            neighborhoods[index].size--;
                        }else{
                            neighbor = true;
                        }
                        // Add to neighbor list
                        if(neighborhoods[index].size < max_neighborhood_size){
                            neighborhoods[index].indices[neighborhoods[index].size]= j;
                            neighborhoods[index].costs[neighborhoods[index].size] = dist(pt,points[j]);
                            neighborhoods[index].cache[neighborhoods[index].size] = CACHE_FREE;
                            neighborhoods[index].size += 1;
                        }

                        // Since r-connected graphs are bidirectional, add to other neighborhood
                        if(neighborhoods[j].size < max_neighborhood_size){
                            neighborhoods[j].indices[neighborhoods[j].size] = index;
                            neighborhoods[j].costs[neighborhoods[j].size]   = dist(pt,points[j]);
                            neighborhoods[j].cache[neighborhoods[j].size]   = CACHE_FREE;
                            neighborhoods[j].size += 1;
                        }else{
    #ifdef FMT_DEBUG
                            printf("Couldn't add %d as neighbor to %d because neighborhood full\n", index,j);
    #endif
 
                        }
                    }
                }
            }
        }
#ifdef FMT_DEBUG
        printf("neighborhood has %d elements\n", neighborhoods[index].size);
#endif
        
        // TODO: Sort the neighborhood by distance (? not sure if it will improve performance at all)
        return true;
    }else if(parameters.connection_type == KNN_CON){
        /*
         // K nearest neighbors connected (this is slow to compute)
         int K = int(parameters.connection_param) + 1; // +1 compensates for including self in neighborhood
         // Array with nearest distances:
         float near_dists[K];
         int near_indices[K];
         for(int k = 0; k < K; k++)
         near_dists[k] = MAXFLOAT;
         
         // Compute distances to all points and track the closest K
         float d_ij = 0;
         for(int j = 0; j < parameters.num_pts; j++){
         d_ij = dist(points[index],points[j]);
         // If smaller than current best, insert into best vectors.
         for(int k = 0; k < K; k++){
         if(d_ij < near_dists[k])
         {
         //shift & insert at k
         for(int kk = K; kk > k; kk--){
         near_dists[kk] = near_dists[kk-1];
         near_indices[kk] = near_indices[kk-1];
         }
         near_dists[k] = d_ij; near_indices[k] = j;
         }
         break;
         }
         }
         // Insert K nearest neighbors into neighborhood:
         neighborhoods[index].indices.reserve(K);
         neighborhoods[index].costs.reserve(K);
         for(int k = 0; k < K; k++){
         neighborhoods[index].indices[k] = near_indices[k];
         neighborhoods[index].costs[k]   = near_dists[k];
         }
         
         // Now need to check all nearby neighbors. That's going to suck!
         // TODO: Do this.
         
         return true;
         */
        printf("Error: KNN not implemented\n");
        return false;
    }else{
        printf("Error (1): Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }
}



// Resets neighborhoods. Meant for use only by start_ind and goal_ind
bool FMT::reset_neighborhood(int index){
#ifdef FMT_DEBUG
    printf("reset_neighborhood\n");
#endif
    
    if(index != start_pt_ind && index != goal_pt_ind){
        printf("Error: reset neighborhoods only meant for start and goal points");
        return false;
    }
    
    if(parameters.connection_type == RAD_CON){
        // clear cache
        // Remove from neighborhood of every neighbor. Because these are added at the end, only need to check last two:
        for(int i_n = 0; i_n < neighborhoods[index].size; i_n++){
//            neighborhoods[index].cache[i_n] = CACHE_UNK;
            // check neighbor's end points for target:
            int n = neighborhoods[index].indices[i_n];
            int n_size = neighborhoods[n].size;
            if(n_size >= 1 && neighborhoods[n].indices[n_size-1] == index){
                //Remove end elements by forgetting they exist:
                neighborhoods[n].size -=1;
            }else if(n_size >= 2 && neighborhoods[n].indices[n_size-2] == index){
                //Shift end elements and remove end:
                neighborhoods[n].indices[n_size-2]=neighborhoods[n].indices[n_size-1];
                neighborhoods[n].costs[n_size-2]=neighborhoods[n].costs[n_size-1];
                neighborhoods[n].size -=1;
            }
        }
        // Reset neighborhood at index
        neighborhoods[index].size = 0;
        
    }else if(parameters.connection_type == KNN_CON){
        printf("Error: KNN not fully implemented yet!\n");
        return false;
    }else{
        printf("Error (2) : Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }
    
    return true;
}


// Computes neighborhoods for sampled points (NOT initial, final point)
bool FMT::compute_neighborhoods(){
    
    neighborhoods = new Neighborhood[(parameters.num_pts+2)];
    
    //Initialize neighborhoods
    for(int i = 0; i < parameters.num_pts+2; i++){
        neighborhoods[i].size = 0;
        neighborhoods[i].indices = new int[max_neighborhood_size];
        neighborhoods[i].costs   = new float[max_neighborhood_size];
        neighborhoods[i].cache   = new uint8_t[max_neighborhood_size];
        // Likely unnecessary.
        for(int j = 0; j < max_neighborhood_size; j++){
            neighborhoods[i].cache[j] = CACHE_UNK;
        }
    }
    if(parameters.connection_type == RAD_CON){
        // Radially connected neighbors
        float dim = 2.0;
        if(parameters.Z_limit > 0)
            dim = 3.0;
        
        // Coefficient from equation (3) of FMT* paper (Janson, Schmerling, et al.)
        // Not sure what xi is in their equation (tuning parameter??)
        float coeff = (parameters.X_limit*parameters.Y_limit*log(parameters.num_pts))/(dim*parameters.num_pts*1.1);
        parameters.connection_param = 2.2*powf(coeff, 1/dim);
        
        // Radially connected neighbors
        for(int i = 0; i < parameters.num_pts; i++){
            for(int j = 0; j < parameters.num_pts; j++){
                if(dist(points[i],points[j]) < parameters.connection_param)
                {
                    // Add to neighbor list
                    neighborhoods[i].indices[neighborhoods[i].size] = j;
                    neighborhoods[i].costs[neighborhoods[i].size] = dist(points[i],points[j]);
                    neighborhoods[i].cache[neighborhoods[i].size] = CACHE_UNK;
                    neighborhoods[i].size += 1;
                    if(neighborhoods[i].size >= max_neighborhood_size){
#ifdef FMT_DEBUG
                        printf("Warning: Neighborhood larger than specified maximum size\n");
#endif
                        break; // This is a hack
                    }
                }
            }
        }
        // Everything is sorted by index here. Could also add a remap for sorting by cost?
        return true;
    }else if(parameters.connection_type == KNN_CON){
        /*
         // K nearest neighbors connected (this is slow to compute)
         int K = int(parameters.connection_param) + 1; // +1 compensates for including self in neighborhood
         for(int i = 0; i < parameters.num_pts; i++){
         // Array with nearest distances:
         float near_dists[K];
         int near_indices[K];
         for(int k = 0; k < K; k++)
         near_dists[k] = MAXFLOAT;
         
         // Compute distances to all points and track the closest K
         float d_ij = 0;
         for(int j = 0; j < parameters.num_pts; j++){
         d_ij = dist(points[i],points[j]);
         // If smaller than current best, insert into best vectors.
         for(int k = 0; k < K; k++){
         if(d_ij < near_dists[k])
         {
         //shift & insert at k
         for(int kk = K; kk > k; kk--){
         near_dists[kk] = near_dists[kk-1];
         near_indices[kk] = near_indices[kk-1];
         }
         near_dists[k] = d_ij; near_indices[k] = j;
         }
         break;
         }
         }
         // Insert K nearest neighbors into neighborhood:
         neighborhoods[i].indices.reserve(K);
         neighborhoods[i].costs.reserve(K);
         for(int k = 0; k < K; k++){
         neighborhoods[i].indices[k] = near_indices[k];
         neighborhoods[i].costs[k]   = near_dists[k];
         }
         }
         */
        printf("Error: KNN not implemented\n");
        return false;
    }else{
        printf("Error (3): Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }
}


void FMT::clear_cache(){
    for(int i = 0; i < parameters.num_pts; i++){
        for(int j = 0; j < neighborhoods[i].size; j++){
            if(neighborhoods[i].cache[j] == CACHE_FREE_U || neighborhoods[i].cache[j] == CACHE_OCC_U){
                neighborhoods[i].cache[j] = CACHE_UNK;
            }
        }
    }
}

// Checks for collision using cache to speed up.
bool FMT::collision_inds(int ind1, int ind2, Map *map){
    num_total_checks++;
    // ====  First find neighbor indices: ==== //
    int neighb_1 = -1; // Index of ind2 in ind1's neighbor list
    int neighb_2 = -1; // Index of ind1 in ind2's neighbor list
    
    neighb_1 = find_element(neighborhoods[ind1].indices, neighborhoods[ind1].size, ind2);
    neighb_2 = find_element(neighborhoods[ind2].indices, neighborhoods[ind2].size, ind1);
    
    // ==== Check for collisions unless an answer is already cached ==== //
    // Case 1: Not in each others' neighborhoods so cache is irrelevant
    if(neighb_1 == -1 && neighb_2 == -1){
        return collision(points[ind1], points[ind2], map);
    }
    
    // If one has collision cached, treat as collision.
    if( (neighb_1 != -1 && neighborhoods[ind1].cache[neighb_1] >= CACHE_OCC) || (neighb_2 != -1 && neighborhoods[ind2].cache[neighb_2] >= CACHE_OCC)){
        num_skipped_checks++;
        return true;
    }
    
    // If one has free (we already know neither has collision), then return false
    if( (neighb_1 != -1 && (neighborhoods[ind1].cache[neighb_1] == CACHE_FREE || neighborhoods[ind1].cache[neighb_1] == CACHE_FREE_U)) || (neighb_2 != -1 && (neighborhoods[ind2].cache[neighb_2] == CACHE_FREE || neighborhoods[ind2].cache[neighb_2] == CACHE_FREE_U))){
        num_skipped_checks++;
        return false;
    }

    // If we're here, then the cache is "CACHE_UNK"
    // Compute the collision value and update the cache
    bool confirmed_value = false;
    bool collision_value = collision(points[ind1], points[ind2], map, &confirmed_value);
    if(collision_value == true){
        if(confirmed_value){
            if(neighb_1 != -1)
                neighborhoods[ind1].cache[neighb_1] = CACHE_OCC;
            if(neighb_2 != -1)
                neighborhoods[ind2].cache[neighb_2] = CACHE_OCC;
        }else{
            if(neighb_1 != -1)
                neighborhoods[ind1].cache[neighb_1] = CACHE_OCC_U;
            if(neighb_2 != -1)
                neighborhoods[ind2].cache[neighb_2] = CACHE_OCC_U;
        }
    }else{
        if(confirmed_value){
            if(neighb_1 != -1)
                neighborhoods[ind1].cache[neighb_1] = CACHE_FREE;
            if(neighb_2 != -1)
                neighborhoods[ind2].cache[neighb_2] = CACHE_FREE;
        }else{
            if(neighb_1 != -1)
                neighborhoods[ind1].cache[neighb_1] = CACHE_FREE_U;
            if(neighb_2 != -1)
                neighborhoods[ind2].cache[neighb_2] = CACHE_FREE_U;
        }
    }
    
    return collision_value;
}



void FMT::print_point(Point pt){
    printf("(%f, %f, %f)",pt.x,pt.y,pt.z);
}

void FMT::print_neighborhood(Neighborhood n){
    printf("[%d",n.indices[0]);
    for(int i = 1; i < n.size; i++)
        printf(", %d",n.indices[i]);
    printf("]");
}

void FMT::print_parameters(){
    printf("X_limit: %d\n", parameters.X_limit);
    printf("Y_limit: %d\n", parameters.Y_limit);
    printf("Z_limit: %d\n", parameters.Z_limit);
    
    printf("num_pts: %d\n", parameters.num_pts);
    printf("connection_type: %d\n", parameters.connection_type);
    printf("connection_param: %d\n", parameters.connection_type);
    
}

void FMT::print_array(int* vec, int size){
    printf("[ ");
    for(int i = 0; i < size; i++)
        printf("%d ", vec[i]);
    printf("]");
}

// Used for plotting in Julia
void FMT::print_datafile(Point *path, int path_length, Map* map){
    FILE* datafile = fopen("Datafile.jl", "w");
    
    // Print points
    fprintf(datafile, "points = [");
    for(int i = 0; i < parameters.num_pts+2; i++)
        fprintf(datafile, "%f %f\n", points[i].x, points[i].y);
    fprintf(datafile, "]\n");
    
    //Print neighborhoods
    fprintf(datafile, "adjacency = [");
    for(int i = 0; i < parameters.num_pts+2; i++){
        for(int j = 0; j < parameters.num_pts+2; j++){
            bool connected = false;
            for(int k = 0; k < neighborhoods[i].size; k++){
                if(neighborhoods[i].indices[k] == j){
                    connected = true;
                    break;
                }
            }
            fprintf(datafile, "%d ", connected);
        }
        if(i==parameters.num_pts+1){
            fprintf(datafile,"]");
        }
        fprintf(datafile, "\n");
    }

    
    //Print path
    fprintf(datafile, "path_length = %f\n", path[0].x);
    fprintf(datafile, "cost  = %f\n", path[0].y);
    fprintf(datafile, "path = [");
    for(int i = 1; i < path_length+1; i++)
        fprintf(datafile, "%f %f\n", path[i].x, path[i].y);
    fprintf(datafile, "]\n");
    fprintf(datafile, "start_pt = [%f %f]\n", points[start_pt_ind].x, points[start_pt_ind].y);
    fprintf(datafile, "goal_pt  = [%f %f]\n", points[goal_pt_ind].x, points[goal_pt_ind].y);
    
    fprintf(datafile, "map = [");
    for(int i = 0; i < parameters.X_limit; i++){
        for(int j = 0; j < parameters.Y_limit; j++){
            fprintf(datafile, "%d ", (*map)(i,j,0));
        }
        if(!(i+1 < parameters.X_limit))
            fprintf(datafile, "]");
        fprintf(datafile, "\n");
    }
    
    fprintf(datafile, "map_scan = [");
    Point pt1; pt1.z = 0.0;
    Point pt2 = points[goal_pt_ind];
    
    int coll_cnt = 0;
    for(int i = 0; i < parameters.X_limit; i++){
        for(int j = 0; j < parameters.Y_limit; j++){
            bool col = false;
                pt1.x = float(i); pt1.y = float(j);
                col = collision(pt1, pt2, map);
                if(col)
                    coll_cnt++;
                fprintf(datafile, "%d ", col);
        }
        if(!(i+1 < parameters.X_limit))
            fprintf(datafile, "]");
        fprintf(datafile, "\n");
    }
    printf("Collisions %d\n", coll_cnt);
    
}


void FMT::print_stats(){
    printf("Total collision checks: %d. skipped checks: %d. Ratio: %f\n", num_total_checks, num_skipped_checks, float(num_skipped_checks)/float(num_total_checks));
}