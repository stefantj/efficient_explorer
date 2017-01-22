//
//  main.cpp
//  FMT
//
//  Created by Stefan Jorgensen
//  MIT licence
//

#include "fmt.h"


// Temporary function for printing the neighborhood of an agent:
void FMT::plot_neighborhood(int i){
    Neighborhood* n = &(neighborhoods[i]);
    Point* pt = points[i].state;
    printf("# Begin neighborhood plot\n");
    printf("scatter3D(%f,%f,%f,marker=\"o\",color=:green);\n", pt[0].x,pt[0].y,pt[0].z); // Mark current location
    printf("plot3D( [%f,%f],[%f,%f],[%f,%f],marker=\"v\",color=:red);", pt[0].x, pt[0].x+pt[1].x, pt[0].y, pt[0].y+pt[1].y, pt[0].z, pt[0].z+pt[1].z);

    printf("plot3D( [%f,%f],[%f,%f],[%f,%f],color=:black);\n", pt[0].x, pt[0].x+pt[2].x, pt[0].y, pt[0].y+pt[2].y, pt[0].z, pt[0].z+pt[2].z);

    
    PolyState* p;
    Point* init = pt;
    Point* final;
    // Print the polynomial for each point
    for(int j = 0; j < n->size; j++){
        p = &(n->paths)[j];
        
        final = points[n->indices[j]].state;
        
        // Load coefficients for polynomial:
        
        int num_init = FMT_ORD;
        int num_final= FMT_ORD;
        
        printf("x_coeffs = [");
        for(int k= 0; k < p->order-1; k++)
            printf("%f,",p->coefficients_x[k]);
        printf("%f]\n", p->coefficients_x[p->order-1]);
        
        printf("y_coeffs = [");
        for(int k= 0; k < p->order-1; k++)
            printf("%f,",p->coefficients_y[k]);
        printf("%f]\n", p->coefficients_y[p->order-1]);
        
        printf("z_coeffs = [");
        for(int k= 0; k < p->order-1; k++)
            printf("%f,",p->coefficients_z[k]);
        printf("%f]\n", p->coefficients_z[p->order-1]);
        
        printf("t = %f\n", p->duration);
        printf("init = [");
        for(int k = 0; k < num_init-1; k++)
            printf("%f %f %f;",init[k].x,init[k].y,init[k].z);
        printf("%f %f %f]\n", init[num_init-1].x,init[num_init-1].y,init[num_init-1].z);
        printf("final = [");
        for(int k = 0; k < num_final-1; k++)
            printf("%f %f %f;",final[k].x,final[k].y,final[k].z);
        printf("%f %f %f]\n", final[num_final-1].x,final[num_final-1].y,final[num_final-1].z);
        
        
        // Plot polynomial path
        printf("times = linspace(0,t,1000);\n");
        printf("px = zeros(times); py = zeros(times); pz = zeros(times);\n");
        printf("for k=1:size(x_coeffs,1)\n");
        printf("    px += x_coeffs[k]*times.^(k-1)\n");
        printf("    py += y_coeffs[k]*times.^(k-1)\n");
        printf("    pz += z_coeffs[k]*times.^(k-1)\n");
        printf("end\n");
        printf("plot3D(px,py,pz);\n");
        printf("waitforbuttonpress()\n");
    }
    

    
    // Print the initial point, with barbs for where it is:
    
}

// Explicit constructor method
FMT::FMT(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param, Map* map_structure){
    // set parameters
    set_parameters(xlimit, ylimit, zlimit, num_pts, connection_type, connection_param);
    // set velocities manually - should rework interface here.
    parameters.X_vel_limit = 5.0; // meters/second
    parameters.X_acc_limit = 2.0; // meters/second^2
    parameters.Y_vel_limit = 5.0; // meters/second
    parameters.Y_acc_limit = 2.0; // meters/second^2
    parameters.Z_vel_limit = 0.5; // meters/second
    parameters.Z_acc_limit = 0.3; // meters/second^2
    
    if(connection_type == RAD_CON){
        max_neighborhood_size = int(num_pts/5);
    }else if(connection_type == KNN_CON){
        max_neighborhood_size = int(connection_param+10);
    }
    
    smoother = new PolynomialSmoother();

    // Preallocate vectors (+2 because we add start and end points at the end)
    points = new FMT_Point[(parameters.num_pts+2)];
    sample_points();
    
    // Compute neighborhoods (return if error)
    if(!compute_neighborhoods(map_structure))
        return;
#ifdef FMT_DEBUG
    print_parameters();
    printf("compute neighborhoods successful\n");
//    for(int i = 0; i < parameters.num_pts; i++)
//        plot_neighborhood(i);

#endif
    
    
    int max_num_cells = 0;
    float mean_cells = 0;
    float cells_count = 0;
    for(int i = 0; i < parameters.num_pts+2; i++){
        for(int j = 0; j < neighborhoods[i].size; j++){
            if(neighborhoods[i].paths[j].num_cells > max_num_cells)
                max_num_cells = neighborhoods[i].paths[j].num_cells;
            mean_cells = mean_cells*(cells_count/(cells_count+1)) + (float)neighborhoods[i].paths[j].num_cells/(cells_count+1);
            cells_count+=1;
        }
    }
    
    printf(" Maximum number of cells used: %d. Maximum allocated: %d\n", max_num_cells, MAX_POLY_CELLS);
    printf(" Average number of cells used: %f. Wasted allocation: %f\n", mean_cells, (MAX_POLY_CELLS-mean_cells)*cells_count);
    
    
    
    
    // Initialize planning containers: // all of these could be made smaller for memory savings
    Windex = new int[parameters.num_pts+1];
    Parents = new int[parameters.num_pts+2];
    Costs = new float[parameters.num_pts+2];
    Active = new int[parameters.num_pts+2];
    H_new = new int[parameters.num_pts+2];

    x_near.costs = new float[parameters.num_pts+2];
    x_near.indices = new int[parameters.num_pts+2];
    y_near.costs = new float[parameters.num_pts+2];
    y_near.indices = new int[parameters.num_pts+2];
    
    
    init_vel.x = 0.0;
    init_vel.y = 0.0;
    init_vel.z = 0.0;
    
    init_acc.x = 0.0;
    init_acc.y = 0.0;
    init_acc.z = 0.0;
    

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
            for(int j = 0; j < neighborhoods[i].size; j++)
                delete [] neighborhoods[i].paths[j].cells;
            delete [] neighborhoods[i].paths;
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

int FMT::fmtstar(Point start, Point goal, Map* map, PolyState* path){
    
    //push start, goal onto points list
    points[start_pt_ind].state[0].x = start.x;
    points[start_pt_ind].state[0].y = start.y;
    points[start_pt_ind].state[0].z = start.z;
    
    points[goal_pt_ind].state[0].x  = goal.x;
    points[goal_pt_ind].state[0].y  = goal.y;
    points[goal_pt_ind].state[0].z  = goal.z;

    // Error check (this is incomplete)
    if(isnan(start.x) || isnan(start.y) || isnan(start.z) || isnan(goal.x) || isnan(goal.y) || isnan(goal.z)){
        path[0].cost = -1;
        return 0;
    }

    // Start by constructing neighborhoods for start/finish points
    if(! compute_neighborhood(goal, goal_pt_ind, map) || !compute_neighborhood(start, start_pt_ind, map) ){
        path[0].cost = -1;
        return 0;
    }
    
    
#ifdef FMT_DEBUG
    printf("Start index has %d neighbors, end index has %d. \n", neighborhoods[start_pt_ind].size, neighborhoods[goal_pt_ind].size);
    
    int n_have_start = 0;
    int n_have_end = 0;
    for(int i = 0; i < parameters.num_pts; i++){
        if(find_element(neighborhoods[i].indices, neighborhoods[i].size, start_pt_ind) != -1)
            n_have_start++;
        
        if(find_element(neighborhoods[i].indices, neighborhoods[i].size, goal_pt_ind) != -1)
            n_have_end++;
    }
    
    printf("%d other points have the start point as a neighbor, %d have the end point. \n", n_have_start, n_have_end);
    
    bool feasible = false;
    for(int i = 0; i < neighborhoods[start_pt_ind].size; i++){
        if(!collision_inds(start_pt_ind, neighborhoods[start_pt_ind].indices[i], map))
            feasible = true;
    }
    if(!feasible){
        printf("No connections to starting point. ");
    }
    feasible = false;
    for(int i = 0; i < neighborhoods[goal_pt_ind].size; i++){
        if(!collision_inds(goal_pt_ind, neighborhoods[goal_pt_ind].indices[i], map))
            feasible = true;
    }
    if(!feasible){
        printf("No connections to end point.\n");
    }
    
#endif
    
    return plan_path(map, path);
}

// Plans the path. Need to change the path representation.
int FMT::plan_path(Map* map, PolyState* path){
    
    // Initialize containers for search
    int Windex_size = parameters.num_pts+2; // includes end point but not start point.
    for(int i = 0; i < Windex_size; i++)
        Windex[i] = i;
    
    int z = start_pt_ind;
    
    int Active_size = 1;
    Active[0] = z;
    
    for(int i = 0; i < parameters.num_pts+2; i++){
        Costs[i] = 0.0;
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

#ifdef FMT_DEBUG
    int num_points_checked = 0;
#endif

    while(z!=goal_pt_ind){
        // Pull out current element
        Windex_size = remove_element(Windex, Windex_size, z);

        // Reset containers:
        H_new_size = 0;
        
        //Construct neighborhood of z
        filter(&x_near, z, Windex, Windex_size, -1);
        
#ifdef FMT_DEBUG
        int cnt = 0;
        bool sorted = true;
        for(int i = 0; i < neighborhoods[z].size; i++){
            if(find_element(Windex, Windex_size, neighborhoods[z].indices[i]) >= 0)
                cnt++;
            for(int j = 0; j < i; j++){
                if(neighborhoods[z].indices[i] <= neighborhoods[z].indices[j])
                    sorted = false;
            }
        }
        
        printf("X_near (%d) has %d points. Should have %d points. Sorted? %d\n", z, x_near.size, cnt, sorted);
        printf("Windex has %d points\n", Windex_size);
        printf("Z_neighborhood has %d points\n", neighborhoods[z].size);
#endif

        int ind_x = 0;
        for(int i_x = 0; i_x < x_near.size; i_x++){
            ind_x = x_near.indices[i_x];
            
            
            /*
            // I don't think this transfers quite right to the KNN setting. Problem is this searches for connections _from_ x points _to_ active, when we want the vice-versa connections.
            filter(&y_near, ind_x, Active, Active_size, ind_x);
            
            if(y_near.size > 0){
                // Find minimum parent:
                float cmin = MAXFLOAT;
                int y_min = -1;
                for(int i_y = 0; i_y < y_near.size; i_y++){
                    if(ind_x == y_near.indices[i_y]) // Don't connect to self!
                        continue;
                    
                    // TODO: Update
                    float c = Costs[y_near.indices[i_y]]+y_near.costs[i_y];
                    if(c < cmin){
                        cmin  = c;
                        y_min = i_y;
                    }
                }
             */
            int min_parent = -1;
            float c_min = MAXFLOAT;
            // Search each element in Active to see whether ind_x is in it. This is _much_ slower!! Re-implement a backwards neighborhood if this fixes things.
            for(int a = 0; a < Active_size; a++){
                if(Active[a] == ind_x)
                    continue;
                
                int active_ind = find_element(neighborhoods[Active[a]].indices, neighborhoods[Active[a]].size, ind_x);
                if(active_ind != -1){
                    float c = Costs[Active[a]] + neighborhoods[Active[a]].costs[active_ind];
                    if(c < c_min){
                        c_min = c;
                        min_parent = Active[a];
                    }
                }
            }
        
            if(min_parent != -1){
#ifdef FMT_DEBUG
                printf("Found parent for point %d -> %d\n", ind_x, min_parent);
#endif
                if(min_parent != -1 && !collision_inds(ind_x, min_parent, map)){
#ifdef FMT_DEBUG

#endif
                    // Add parent
                    Parents[ind_x] = min_parent;
                    
                    // Update cost to go
                    Costs[ind_x] = c_min;
                    
                    // Mark as visited
                    H_new[H_new_size] = ind_x;
                    H_new_size += 1;

                    Windex_size = remove_element(Windex, Windex_size, ind_x);
                    // Add to heap
                    Heap_pt new_pt;
                    new_pt.index = ind_x;
                    new_pt.cost  = c_min;
                    HHeap.push(new_pt);
                } else{
#ifdef FMT_DEBUG
                    printf("collision? %d. Otherwise y_near is empty.\n", collision_inds(ind_x, min_parent, map));
#endif
                }
            }else{
                printf("y_near is empty (no connection to %d active points)\n", Active_size);
            }
        }

        // Update active list
        Active_size = join_sets(Active, Active_size, H_new, H_new_size);

        // Remove z
        Active_size = remove_element(Active,Active_size, z);
#ifdef FMT_DEBUG
        printf("Active size: %d\n", Active_size);
        for(int a = 0; a < Active_size; a++){
            int goal_connection = find_element(neighborhoods[Active[a]].indices, neighborhoods[Active[a]].size, goal_pt_ind);
            if(goal_connection >= 0){
                printf("Goal is reachable by an element of Active (%d)\n", Active[a]);
                break;
            }
        }
#endif

        if(HHeap.size() > 0){
            z = HHeap.top().index;
            HHeap.pop();
#ifdef FMT_DEBUG
            num_points_checked++;
#endif
        }else{
#ifdef FMT_DEBUG
            printf("Heap is empty - failure. Checked %d points\n", num_points_checked);
#endif
            path[0].cost = -1;
            return 0;
        }
    }
    
    // Construct path by assembling parents from goal to start:
    // Going to end up backwards this way
    int path_length = 0;    // Number of segments
    int curr_ind = 0;
    int start_ind = 0;
    path[0].cost = 0;
    // Connect to goal if LOS was used
    if(z!= goal_pt_ind){
        curr_ind = goal_pt_ind;
        start_ind = z;
        copy_polystate(&(path[parameters.num_pts+1]), (get_neighborhood_poly(start_ind, curr_ind)));
#ifdef FMT_DEBUG
        printf("(LOS) Placed polynomial at %d\n", parameters.num_pts+1);
#endif
        path_length = 1;
    }
    
// Doesn't seem right... DEFINITE bugs here.
    curr_ind = z;
    // Insert points backwards from the back
    while(Parents[curr_ind] >= 0){
        // Error checking:
        if(path_length >= parameters.num_pts+2) // Error in path - too many elements
        {path_length = -1; break;}
        if(curr_ind == Parents[curr_ind]) // Error in path - loop
        {path_length = -1; break;}
        
        //Insert appropriate polynomial in appropriate place of path:
        start_ind = Parents[curr_ind];
        copy_polystate(&(path[(parameters.num_pts+1)- path_length]), get_neighborhood_poly(start_ind, curr_ind));
#ifdef FMT_DEBUG
        printf("Placed polynomial at %d\n", parameters.num_pts+1-path_length);
#endif
        // Update costs
        if(path[(parameters.num_pts+1)-path_length].cost == -1){
            path_length = -1; break; // New segment is infeasible or not found
        }else{
            path[0].cost += path[(parameters.num_pts+1)-path_length].cost;
        }
        curr_ind = start_ind;
        path_length++;
    }

    if(path_length >= 0){
        if(start_ind != start_pt_ind){
            copy_polystate(&(path[1]), get_neighborhood_poly(start_ind, curr_ind));
            // Update costs
            if(path[1].cost == -1){
                path[0].cost = -1;
            }else{
                path[0].cost += path[(parameters.num_pts+1)-path_length].cost;
            }
            path_length++;
        }

        // Shift to front
        int shift_ind = parameters.num_pts+1;
        for(int i = 0; i < path_length; i++){
#ifdef FMT_DEBUG
            printf("Shifting from %d to %d\n", shift_ind - i, i+1);
#endif
            copy_polystate(&(path[i+1]), &(path[shift_ind-i]));
            path[i+1] = path[shift_ind - i];
        }
        
    }else{
        // Path generation failed - return badpath.
        path[0].cost = -1;
        return 0;
    }

#ifdef FMT_DEBUG
    printf("Polynomial has %d segments with coefficients :\n",path_length);
    for(int i = 0; i <=parameters.num_pts+1; i++){
        if(path[i].order != 0){
            printf("Seg %d: order %d (rev= %d)\n", i, path[i].order, path[i].reverse);
            printf("x = [ ");
            for(int k = 0; k < path[i].order; k++)
                printf(" %f ", path[i].coefficients_x[k]);
            printf("]\ny = [ ");
            for(int k = 0; k < path[i].order; k++)
                printf(" %f ", path[i].coefficients_y[k]);
            printf("]\n");
            printf("cells: ");
            for(int k = 0; k < path[i].num_cells; k++)
                printf("%d, ", path[i].cells[k]);
            printf("\n");
        }
    }
    
    printf("Path length: %d, path cost: %f\n", path_length, path[0].cost);
#endif
    
    // Clean up planner containers
    reset_neighborhood(start_pt_ind);
    reset_neighborhood(goal_pt_ind);

    return path_length;
}

PolyState* FMT::get_neighborhood_poly(int start_ind, int end_ind){
    // See whether in start_ind's neighborhood:
    for( int i_s=0; i_s < neighborhoods[start_ind].size; i_s++){
        if(neighborhoods[start_ind].indices[i_s] == end_ind){
            return &(neighborhoods[start_ind].paths[i_s]);
        }
    }
    // Otherwise return bad polynomial
    return NULL;
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
        points[i].state[0].x = rand()*x_scaling;
        points[i].state[0].y = rand()*y_scaling;
    }
    
    // Only call rand() if we're planning in Z
    if(parameters.Z_limit == 0){
        for(int i = 0; i < parameters.num_pts; i++)
            points[i].state[0].z = 0.0;
    }else{
        float z_scaling = parameters.Z_limit/RAND_MAX;
        for(int i = 0; i < parameters.num_pts; i++)
            points[i].state[0].z = rand()*z_scaling;
    }
    for(int order = 1; order < FMT_ORD; order++){
        // Generate random numbers
        if(order ==1){
            x_scaling = float(2.0*parameters.X_vel_limit)/RAND_MAX;
            y_scaling = float(2.0*parameters.Y_vel_limit)/RAND_MAX;
        }else{
            x_scaling = float(2.0*parameters.X_acc_limit)/RAND_MAX;
            y_scaling = float(2.0*parameters.Y_acc_limit)/RAND_MAX;
        }
        for(int i = 0; i < parameters.num_pts; i++)
        {
            points[i].state[order].x = rand()*x_scaling-x_scaling/2.0;
            points[i].state[order].y = rand()*y_scaling-y_scaling/2.0;
        }
        
        // Only call rand() if we're planning in Z
        if(parameters.Z_limit == 0){
            for(int i = 0; i < parameters.num_pts; i++)
                points[i].state[order].z = 0.0;
        }else{
            float z_scaling = float(2.0*parameters.Z_vel_limit)/RAND_MAX;
            for(int i = 0; i < parameters.num_pts; i++)
                points[i].state[order].z = rand()*z_scaling-z_scaling/2.0;
        }
    }
}


// Resets neighborhoods. Meant for use only by start_ind and goal_ind
bool FMT::reset_neighborhood(int index){
#ifdef FMT_DEBUG
//    printf("reset_neighborhood\n");
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
            
            // Check last item
            if(n_size >= 1 && neighborhoods[n].indices[n_size-1] == index){
                // Delete the path:
                //                neighborhoods[n].paths[neighborhoods[n].size].num_cells = 0;
                //Remove end elements by forgetting they exist:
                neighborhoods[n].size -=1;
            }else if(n_size >= 2 && neighborhoods[n].indices[n_size-2] == index){
                //Shift end elements and remove end:
                neighborhoods[n].indices[n_size-2]=neighborhoods[n].indices[n_size-1];
                neighborhoods[n].costs[n_size-2]=neighborhoods[n].costs[n_size-1];
                //                neighborhoods[n].paths[neighborhoods[n].size].num_cells= 0;
                neighborhoods[n].size -=1;
            }
            //            neighborhoods[index].paths[i_n].num_cells = 0;
        }
        // Reset neighborhood at index
        neighborhoods[index].size = 0;
    }else if(parameters.connection_type == KNN_CON){
            // clear cache
            // Remove from neighborhood of every neighbor. Because these are added at the end, only need to check last two:
            for(int i_n = 0; i_n < neighborhoods[index].size; i_n++){
                //            neighborhoods[index].cache[i_n] = CACHE_UNK;
                // check neighbor's end points for target:
                int n = neighborhoods[index].indices[i_n];
                int n_size = neighborhoods[n].size;
                
                // Check last item
                if(n_size >= 1 && neighborhoods[n].indices[n_size-1] == index){
                    // Delete the path:
                    //                neighborhoods[n].paths[neighborhoods[n].size].num_cells = 0;
                    //Remove end elements by forgetting they exist:
                    neighborhoods[n].size -=1;
                }else if(n_size >= 2 && neighborhoods[n].indices[n_size-2] == index){
                    //Shift end elements and remove end:
                    neighborhoods[n].indices[n_size-2]=neighborhoods[n].indices[n_size-1];
                    neighborhoods[n].costs[n_size-2]=neighborhoods[n].costs[n_size-1];
                    //                neighborhoods[n].paths[neighborhoods[n].size].num_cells= 0;
                    neighborhoods[n].cache[n_size-2] = neighborhoods[n].cache[n_size-1];
                    copy_polystate(&(neighborhoods[n].paths[n_size-2]),&(neighborhoods[n].paths[n_size-1]));
                    neighborhoods[n].size -=1;
                }
                //            neighborhoods[index].paths[i_n].num_cells = 0;
            }
            // Reset neighborhood at index
            neighborhoods[index].size = 0;
    }else{
        printf("Error (2) : Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }
    
    return true;
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
        // Not sure what xi is in their equation (tuning parameter?)
        float coeff = (parameters.X_limit*parameters.Y_limit*log(parameters.num_pts))/(dim*parameters.num_pts*1.1);
        parameters.connection_param = 2.2*powf(coeff, 1/dim);
        
        for(int j = 0; j < parameters.num_pts && neighborhoods[index].size < max_neighborhood_size; j++){
            if(j==index)
                continue;
            
            PolyState* p = &(neighborhoods[index].paths[neighborhoods[index].size]);
            
            float d = smoother->fit_polynomial(p, points[index].state, POLY_ORD, points[j].state, POLY_ORD, map);

            if(true || d < parameters.connection_param)
            {

                bool collision_value = false;
                for(int i = 0; i < (neighborhoods[index].paths[neighborhoods[index].size]).num_cells; i++){
                    if(!(map->is_free((neighborhoods[index].paths[neighborhoods[index].size]).cells[i]))){
                        collision_value = true;
                        break;
                    }
                }
                collision_value=false;
                
                if(!collision_value){
                    // Add to neighbor list

                    if(neighborhoods[index].size < max_neighborhood_size){
                        neighborhoods[index].indices[neighborhoods[index].size] = j;
                        neighborhoods[index].costs[neighborhoods[index].size] = (neighborhoods[index].paths[neighborhoods[index].size]).cost;
                        neighborhoods[index].cache[neighborhoods[index].size] = CACHE_FREE;
                        neighborhoods[index].size += 1;
                    }
                    
                    // Since r-connected graphs are bidirectional, add to other neighborhood
                    if(neighborhoods[j].size < max_neighborhood_size){
                        printf("adding to neighbor\n");
                        printf("This does not work right (RADCON form neighborhoods)");
                        neighborhoods[j].indices[neighborhoods[j].size] = index;
                        neighborhoods[j].costs[neighborhoods[j].size]   = (neighborhoods[index].paths[neighborhoods[index].size]).cost;
                        neighborhoods[j].cache[neighborhoods[j].size]   = CACHE_FREE;
                        // Need to copy polystate...
                        PolyState* p_j = &(neighborhoods[j].paths[neighborhoods[j].size]);
                        for(int o = 0; o < p->order; o++){
                            p_j->coefficients_x[o] = p->coefficients_x[o];
                            p_j->coefficients_y[o] = p->coefficients_y[o];
                            p_j->coefficients_z[o] = p->coefficients_z[o];
                            p_j->coefficients_p[o] = p->coefficients_p[o];
                        }
                        p_j->duration = p->duration;
                        p_j->order = p->order;
                        p_j->num_cells = p->num_cells;
                        if(p_j->cells == NULL)
                            p_j->cells = new size_t[MAX_POLY_CELLS];
                        for(int c = 0; c < p->num_cells; c++)
                            p_j->cells[c] = p->cells[c];
                        p_j->reverse = !p->reverse;
                        p_j->cost = p->cost;
                        neighborhoods[j].size += 1;
                    }else{
//                        printf("Neighbor full\n");
                    }
                }else{
                    
                }

            }
        }
#ifdef FMT_DEBUG
        if(neighborhoods[index].size <= 10)
            printf("Warning: neighborhood for %d only has %d elements\n", index, neighborhoods[index].size);
#endif
        // TODO: Sort the neighborhood by distance (? not sure if it will improve performance at all)
        return true;

    }else if(parameters.connection_type==KNN_CON){
        
        // Placeholder polynomial
        PolyState p;
        p.cells = new size_t[MAX_POLY_CELLS];
        
        // Set bound to inf
        float bound = MAXFLOAT;
        // Index of the furthest neighbor
        int bound_index = -1;
        
        for(int j = 0; j < parameters.num_pts; j++){
            if(index==j) continue;
            
            // If we either have too few neighbors, or point j might be a nearest neighbor:
            if(neighborhoods[index].size + 1 < parameters.connection_param){
                // Compute the polynomial:
                float d = smoother->fit_polynomial(&p, points[index].state, POLY_ORD, points[j].state, POLY_ORD, map);
                if(d < 0)
                    continue;
                if(neighborhoods[index].size + 1 < parameters.connection_param || d < bound){
                    if(neighborhoods[index].size+1 > parameters.connection_param){
                        if(bound_index == -1){
                            printf("Error - trying to remove -1 indexed element. Size: %d, bound %f\n", neighborhoods[index].size,bound);
                            continue;
                        }
                        // Remove the furthest neighbor:
                        // Since we want to preserve ordering in terms of index, remove by copying up.
                        for(int k = bound_index+1; k < neighborhoods[index].size; k++){
                            // Index of neighbor:
                            neighborhoods[index].indices[k-1] = neighborhoods[index].indices[k];
                            // Cost
                            neighborhoods[index].costs[k-1] = neighborhoods[index].costs[k];
                            // Cache (though all should be unknown anyway)
                            neighborhoods[index].cache[k-1] = neighborhoods[index].cache[k];
                            // Path - might need to implement copy constructor?
                            if(neighborhoods[index].paths[k-1].cells == NULL)
                                printf("%d, %d has null cells in path\n", index, k-1);
                            if(neighborhoods[index].paths[k].cells == NULL)
                                printf("%d, %d has null cells in path\n", index, k);
                            
                            copy_polystate(&(neighborhoods[index].paths[k-1]), &(neighborhoods[index].paths[k]));
                        }
                    }
                    // Add to neighbor list
                    neighborhoods[index].indices[neighborhoods[index].size] = j;
                    neighborhoods[index].costs[neighborhoods[index].size] = p.cost;
                    neighborhoods[index].cache[neighborhoods[index].size] = CACHE_UNK;
                    
                    // Make sure this copies.
                    copy_polystate(&(neighborhoods[index].paths[neighborhoods[index].size]),&p);

                    
                    // Add self to neighbor's list:
                    if(neighborhoods[j].size < max_neighborhood_size){
                        
                        // Need to recompute the polynomial, since reverse flag does not work right.
                        float d = smoother->fit_polynomial(&p, points[j].state, POLY_ORD, points[index].state, POLY_ORD, map);
                        
                        if(d == -1) // Means polynomial is somehow invalid.
                            continue;

                        neighborhoods[j].indices[neighborhoods[j].size] = index;
                        neighborhoods[j].costs[neighborhoods[j].size] = p.cost;
                        neighborhoods[j].cache[neighborhoods[j].size] = CACHE_UNK;
                        copy_polystate(&(neighborhoods[j].paths[neighborhoods[j].size]), &p);
                        // Increment size:
                        neighborhoods[j].size += 1;
                    }
                    
                    
                    
                    // If you did not replace a neighbor, increment size
                    if(neighborhoods[index].size+1 <= parameters.connection_param){
                        neighborhoods[index].size++;
                    }
                    
                    // Now recompute the bound. Could do this somewhat more efficiently, but I don't think it is worth it.
                    bound = 0;
                    for(int k = 0; k < neighborhoods[index].size; k++){
                        if(neighborhoods[index].costs[k] >= bound){
                            bound = neighborhoods[index].costs[k];
                            bound_index = k;
                        }
                    }
                    //printf("Bound set to %f\n", bound);
                }
            }
        }

        delete p.cells;
        
        return true;

    }else{
        printf("Error (1): Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }
}

// Computes neighborhoods for sampled points (NOT initial, final point)
bool FMT::compute_neighborhoods(Map* map_structure){
    
    neighborhoods = new Neighborhood[(parameters.num_pts+2)];
    
    //Initialize neighborhoods
    for(int i = 0; i < parameters.num_pts+2; i++){
        neighborhoods[i].size = 0;
        neighborhoods[i].indices = new int[max_neighborhood_size];
        neighborhoods[i].costs   = new float[max_neighborhood_size];
        neighborhoods[i].cache   = new uint8_t[max_neighborhood_size];
        neighborhoods[i].paths   = new PolyState[max_neighborhood_size];
        // Likely unnecessary.
        for(int j = 0; j < max_neighborhood_size; j++){
            neighborhoods[i].cache[j] = CACHE_UNK;
            neighborhoods[i].paths[j].cells = new size_t[MAX_POLY_CELLS];
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
        
        PolyState* p;
        
        // Radially connected neighbors
        for(int i = 0; i < parameters.num_pts; i++){
            for(int j = 0; j < parameters.num_pts; j++){
                if(i==j) continue;
                
                // Fit a polynomial:
                p = &(neighborhoods[i].paths[neighborhoods[i].size]);
                float d = smoother->fit_polynomial(p, points[i].state, POLY_ORD, points[j].state, POLY_ORD, map_structure);
                if(d < parameters.connection_param)
                {
                    // Add to neighbor list
                    neighborhoods[i].indices[neighborhoods[i].size] = j;
                    neighborhoods[i].costs[neighborhoods[i].size] = p->cost;
                    neighborhoods[i].cache[neighborhoods[i].size] = CACHE_UNK;
                    // Path already copied above.
                    neighborhoods[i].size += 1;
                    if(neighborhoods[i].size >= max_neighborhood_size){
#ifdef FMT_DEBUG
                        printf("Warning: Neighborhood larger than specified maximum size\n");
#endif
                        break; // This is a hack
                    }
                }
            }
            
#ifdef FMT_DEBUG
            if(neighborhoods[i].size <= 10)
                printf("Warning: neighborhood for %d only has %d elements\n", i, neighborhoods[i].size);
#endif
        }
        return true;
    }else if(parameters.connection_type==KNN_CON){
        
        // Placeholder polynomial
        PolyState p;
        p.cells = new size_t[MAX_POLY_CELLS];
        
        for(int i = 0; i < parameters.num_pts; i++){
            // Set bound to inf
            float bound = MAXFLOAT;
            // Index of the furthest neighbor
            int bound_index = -1;

            for(int j = 0; j < parameters.num_pts; j++){
                if(i==j) continue;
                
                // If we either have too few neighbors, or point j might be a nearest neighbor:
                if(neighborhoods[i].size + 1 < parameters.connection_param){
                    // Compute the polynomial:
                    float d = smoother->fit_polynomial(&p, points[i].state, POLY_ORD, points[j].state, POLY_ORD, map_structure);
                    if(d == -1) // Means polynomial is somehow invalid.
                        continue;
                    if(neighborhoods[i].size + 1 <= parameters.connection_param ||  d < bound){
                        if(neighborhoods[i].size+1 > parameters.connection_param){
                            // Remove the furthest neighbor:
                            // Since we want to preserve ordering in terms of index, remove by copying up.
                            for(int k = bound_index+1; k < neighborhoods[i].size; k++){
                                // Index of neighbor:
                                neighborhoods[i].indices[k-1] = neighborhoods[i].indices[k];
                                // Cost
                                neighborhoods[i].costs[k-1] = neighborhoods[i].costs[k];
                                // Cache (though all should be unknown anyway)
                                neighborhoods[i].cache[k-1] = neighborhoods[i].cache[k];
                                // Path - might need to implement copy constructor?
                                copy_polystate(&(neighborhoods[i].paths[k-1]), &(neighborhoods[i].paths[k]));
                            }
                        }
                        
                        // Add to neighbor list
                        neighborhoods[i].indices[neighborhoods[i].size] = j;
                        neighborhoods[i].costs[neighborhoods[i].size] = p.cost;
                        neighborhoods[i].cache[neighborhoods[i].size] = CACHE_UNK;
                        // Make sure this copies.
                        copy_polystate(&(neighborhoods[i].paths[neighborhoods[i].size]),&p);
                        
                        // If you did not replace a neighbor, increment size
                        if(neighborhoods[i].size+1 <= parameters.connection_param){
                            neighborhoods[i].size++;
                        }
                        
                        // Now recompute the bound. Could do this somewhat more efficiently, but I don't think it is worth it.
                        bound = 0;
                        for(int k = 0; k < neighborhoods[i].size; k++){
                            if(neighborhoods[i].costs[k] > bound){
                                bound = neighborhoods[i].costs[k];
                                bound_index = k;
                            }
                        }
                    }else{

                    }
                }
            }
        }
        
        delete p.cells;
        
        return true;
    }else{
        printf("Error (3): Unsupported connection type %d!\n", parameters.connection_type);
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
        // Not in each other's neighborhoods, which is not good.
#ifdef FMT_DEBUG
        printf("Warning - skipping check because not in neighborhoods\n");
#endif
        return true;
    }
    
    // If one has collision cached, treat as collision.
    if( (neighb_1 != -1 && neighborhoods[ind1].cache[neighb_1] >= CACHE_OCC) || (neighb_2 != -1 && neighborhoods[ind2].cache[neighb_2] >= CACHE_OCC)){
//        printf("Cached collision\n");
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
    
    // Pull out the appropriate polynomial:
    
    PolyState* p;
    if(neighb_1 != -1){
        p = &(neighborhoods[ind1].paths[neighb_1]);
    }else{
        p = &(neighborhoods[ind2].paths[neighb_2]);
    }
    
    bool collision_value = false; // Need a collision method which checks indices of the polynomial.
    if(p == nullptr || p->cells == nullptr ){
        collision_value = true;
    }else{
        
        for(int i = 0; i < p->num_cells; i++){
            if(!(map->is_free(p->cells[i]))){
                collision_value = true;
                break;
            }
        }
    }
    
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
    printf("connection_param: %f\n", parameters.connection_param);
    
}

void FMT::print_array(int* vec, int size){
    printf("[ ");
    for(int i = 0; i < size; i++)
        printf("%d ", vec[i]);
    printf("]");
}

// Used for plotting in Julia
void FMT::print_datafile(Point *path, int path_length, Map* map){
    /*
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
 */
}

void FMT::set_initial_state(Point *vel, Point *acc){
    init_vel.x = vel->x;
    init_vel.y = vel->y;
    init_vel.z = vel->z;
    
    init_acc.x = acc->x;
    init_acc.y = acc->y;
    init_acc.z = acc->z;
#if FMT_ORD > 1
    points[start_pt_ind].state[1].x = vel->x;
    points[start_pt_ind].state[1].y = vel->y;
    points[start_pt_ind].state[1].z = vel->z;
#endif
#if FMT_ORD > 2
    points[start_pt_ind].state[2].x = acc->x;
    points[start_pt_ind].state[2].y = acc->y;
    points[start_pt_ind].state[2].z = acc->z;
#endif
}

void FMT::print_stats(){
    printf("Total collision checks: %d. skipped checks: %d. Ratio: %f\n", num_total_checks, num_skipped_checks, float(num_skipped_checks)/float(num_total_checks));
}