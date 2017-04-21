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
    Neighborhood* n = &(neighborhoodsF[i]);
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
        max_neighborhood_size = int(connection_param*1.5);
    }
    
    smoother = new PolynomialSmoother();

    // Preallocate vectors (+2 because we add start and end points at the end)
    points = new FMT_Point[(parameters.num_pts+2)];
    sample_points();
    
    // Compute neighborhoods (return if error)
    if(!initialize_neighborhoods(map_structure))
        return;
#ifdef FMT_DEBUG
    print_parameters();
    printf("compute neighborhoods successful\n");
//    for(int i = 0; i < parameters.num_pts; i++)
//        plot_neighborhood(i);

    
    int max_num_cells = 0;
    float mean_cells = 0;
    float cells_count = 0;
    for(int i = 0; i < parameters.num_pts+2; i++){
        for(int j = 0; j < neighborhoodsF[i].size; j++){
            if(neighborhoodsF[i].paths[j].num_cells > max_num_cells)
                max_num_cells = neighborhoodsF[i].paths[j].num_cells;
            mean_cells = mean_cells*(cells_count/(cells_count+1)) + (float)neighborhoodsF[i].paths[j].num_cells/(cells_count+1);
            cells_count+=1;
        }
        for(int j = 0; j < neighborhoodsR[i].size; j++){
            if(neighborhoodsR[i].paths[j].num_cells > max_num_cells)
                max_num_cells = neighborhoodsR[i].paths[j].num_cells;
            mean_cells = mean_cells*(cells_count/(cells_count+1)) + (float)neighborhoodsR[i].paths[j].num_cells/(cells_count+1);
            cells_count+=1;
        }
    }
     printf(" Maximum number of cells used: %d. Maximum allocated: %d\n", max_num_cells, MAX_POLY_CELLS);
    printf(" Average number of cells used: %f. Wasted allocation: %f\n", mean_cells, (MAX_POLY_CELLS-mean_cells)*cells_count);
#endif
    
    
    
    // Initialize planning containers: // all of these could be made smaller for memory savings
    Windex = new int[parameters.num_pts+2];
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
            delete [] neighborhoodsF[i].indices;
            delete [] neighborhoodsF[i].costs;
            delete [] neighborhoodsF[i].cache;
            for(int j = 0; j < neighborhoodsF[i].size; j++)
                delete [] neighborhoodsF[i].paths[j].cells;
            delete [] neighborhoodsF[i].paths;
            
            delete [] neighborhoodsR[i].indices;
            delete [] neighborhoodsR[i].costs;
            delete [] neighborhoodsR[i].cache;
            for(int j = 0; j < neighborhoodsR[i].size; j++)
                delete [] neighborhoodsR[i].paths[j].cells;
            delete [] neighborhoodsR[i].paths;
        }
        delete [] neighborhoodsF;
        delete [] neighborhoodsR;
        
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
#ifdef FMT_DEBUG
    printf("\n\n");
#endif
    // Error check (this is incomplete)
    if(isnan(start.x) || isnan(start.y) || isnan(start.z) || isnan(goal.x) || isnan(goal.y) || isnan(goal.z)){
        path[0].cost = -1;
        return 0;
    }
    
    
    //push start, goal onto points list
    points[start_pt_ind].state[0].x = start.x;
    points[start_pt_ind].state[0].y = start.y;
    points[start_pt_ind].state[0].z = start.z;
    points[start_pt_ind].state[1].x = init_vel.x;
    points[start_pt_ind].state[1].y = init_vel.y;
    points[start_pt_ind].state[1].z = init_vel.z;
    points[start_pt_ind].state[2].x = init_acc.x;
    points[start_pt_ind].state[2].y = init_acc.y;
    points[start_pt_ind].state[2].z = init_acc.z;
    
    
    points[goal_pt_ind].state[0].x  = goal.x;
    points[goal_pt_ind].state[0].y  = goal.y;
    points[goal_pt_ind].state[0].z  = goal.z;
    points[goal_pt_ind].state[1].x = 0;
    points[goal_pt_ind].state[1].y = 0;
    points[goal_pt_ind].state[1].z = 0;
    points[goal_pt_ind].state[2].x = 0;
    points[goal_pt_ind].state[2].y = 0;
    points[goal_pt_ind].state[2].z = 0;
    
    
    
/*
    // Try to connect directly:
    float d = smoother->fit_polynomial(&(path[1]), points[start_pt_ind].state, POLY_ORD, points[goal_pt_ind].state, POLY_ORD, map);
    if(d >= 0){
        // Check for collision:
        bool free = true;
        for(int c = 0; c < path[1].num_cells; c++)
            free = free && map->is_free(path[1].cells[c]);
        
        if(free){
            printf("Direct shortcut!\n");
            path[0].cost = path[1].cost;
            return 1;
        }
    }


 */
    
    
    
    reset_neighborhood(start_pt_ind);
    reset_neighborhood(goal_pt_ind);
    
    // Start by constructing neighborhoods for start/finish points
    if(! compute_neighborhood(goal_pt_ind, map) || !compute_neighborhood(start_pt_ind, map) ){
        path[0].cost = -1;
        return 0;
    }
    
    return plan_path(map, path);
}


int FMT::fmtstar(Point start, Point goal, Point goal_vel, Map* map, PolyState* path){
#ifdef FMT_DEBUG
    printf("\n\n");
#endif
    // Error check (this is incomplete)
    if(isnan(start.x) || isnan(start.y) || isnan(start.z) || isnan(goal.x) || isnan(goal.y) || isnan(goal.z)){
        path[0].cost = -1;
        return 0;
    }
    

    //push start, goal onto points list
    points[start_pt_ind].state[0].x = start.x;
    points[start_pt_ind].state[0].y = start.y;
    points[start_pt_ind].state[0].z = start.z;
    points[start_pt_ind].state[1].x = init_vel.x;
    points[start_pt_ind].state[1].y = init_vel.y;
    points[start_pt_ind].state[1].z = init_vel.z;
    points[start_pt_ind].state[2].x = init_acc.x;
    points[start_pt_ind].state[2].y = init_acc.y;
    points[start_pt_ind].state[2].z = init_acc.z;
    
    
    points[goal_pt_ind].state[0].x  = goal.x;
    points[goal_pt_ind].state[0].y  = goal.y;
    points[goal_pt_ind].state[0].z  = 0*goal.z;
    points[goal_pt_ind].state[1].x = goal_vel.x;
    points[goal_pt_ind].state[1].y = goal_vel.y;
    points[goal_pt_ind].state[1].z = 0*goal_vel.z;
    points[goal_pt_ind].state[2].x = 0;
    points[goal_pt_ind].state[2].y = 0;
    points[goal_pt_ind].state[2].z = 0;
    /*
    // Try to connect directly:
    float d = smoother->fit_polynomial(&(path[1]), points[start_pt_ind].state, POLY_ORD, points[goal_pt_ind].state, POLY_ORD, map);
    if(d >= 0){
        // Check for collision:
        bool free = true;
        for(int c = 0; c < path[1].num_cells; c++)
            free = free && map->is_free(path[1].cells[c]);
        
        if(free){
            printf("Direct shortcut!\n");
            path[0].cost = path[1].cost;
            return 1;
        }
    }
     */
    
    reset_neighborhood(start_pt_ind);
    reset_neighborhood(goal_pt_ind);

    // Start by constructing neighborhoods for start/finish points
    if(! compute_neighborhood(goal_pt_ind, map) || !compute_neighborhood(start_pt_ind, map) ){
        path[0].cost = -1;
        return 0;
    }
    
    /*
    FILE* f=fopen("adjmat.jl", "w");
    
    fprintf(f, "AdjF = [");
    for(int i = 0; i < parameters.num_pts+2; i++){
        int neighb_ind = 0;
        for(int j = 0; j < parameters.num_pts+2; j++){
            if(neighb_ind < neighborhoodsF[i].size && neighborhoodsF[i].indices[neighb_ind] == j){
                    fprintf(f," 1");
                    neighb_ind+=1;
            }else{
                fprintf(f," 0");
            }
        }
        if(i != parameters.num_pts+1)
            fprintf(f,";");
    }
    fprintf(f,"];\n");
    fprintf(f,"AdjR = [");
    for(int i = 0; i < parameters.num_pts+2; i++){
        int neighb_ind = 0;
        for(int j = 0; j < parameters.num_pts+2; j++){
            if(neighb_ind < neighborhoodsR[i].size && neighborhoodsR[i].indices[neighb_ind] == j){
                fprintf(f," 1");
                neighb_ind+=1;
            }else{
                fprintf(f," 0");
            }
        }
        if(i != parameters.num_pts+1)
            fprintf(f,";");
    }
    fprintf(f,"];\n");
    
    fprintf(f,"nF = [");
    for(int i = 0; i < parameters.num_pts+2; i++){
        for(int j = 0; j < max_neighborhood_size; j++){
            if(j < neighborhoodsF[i].size){
                fprintf(f," %d", neighborhoodsF[i].indices[j]);
            }else{
                fprintf(f," 0");
            }
        }
        if(i != parameters.num_pts+1)
            fprintf(f, ";");
    }
    fprintf(f, "];\n");

    fprintf(f,"nR = [");
    for(int i = 0; i < parameters.num_pts+2; i++){
        for(int j = 0; j < max_neighborhood_size; j++){
            if(j < neighborhoodsR[i].size){
                fprintf(f," %d", neighborhoodsR[i].indices[j]);
            }else{
                fprintf(f," 0");
            }
        }
        if(i != parameters.num_pts+1)
            fprintf(f, ";");
    }
    fprintf(f, "];\n");

    throw 1;
    
     */
    return plan_path(map, path);
}

// Plans the path.
int FMT::plan_path(Map* map, PolyState* path){
    
    // Initialize containers for search
    int Windex_size = parameters.num_pts+2;
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
        filterF(&x_near, z, Windex, Windex_size, -1);
        
        /*
#ifdef FMT_DEBUG
        int cnt = 0;
        bool sorted = true;
        for(int i = 0; i < neighborhoodsF[z].size; i++){
            if(find_element(Windex, Windex_size, neighborhoodsF[z].indices[i]) >= 0)
                cnt++;
            for(int j = 0; j < i; j++){
                if(neighborhoodsF[z].indices[i] <= neighborhoodsF[z].indices[j])
                    sorted = false;
            }
        }
        bool Wsorted = true;
        for(int i = 0; i < Windex_size; i++){
            for(int j = 0; j < i; j++){
                if(Windex[i] <= Windex[j])
                    Wsorted = false;
            }
        }
        
        printf("X_nearF (%d) has %d points. Should have %d points. Sorted? %d\n", z, x_near.size, cnt, sorted);
        if(!sorted){
            for(int i =0; i < neighborhoodsF[z].size; i++){
                printf("%d ", neighborhoodsF[z].indices[i]);
            }
            printf("\n");
        }
        printf("Windex has %d points. Sorted? %d\n", Windex_size,Wsorted);
        printf("Active has %d points\n", Active_size);
#endif
         */

        int ind_x = 0;
        for(int i_x = 0; i_x < x_near.size; i_x++){
            ind_x = x_near.indices[i_x];
            filterR(&y_near, ind_x, Active, Active_size, ind_x);
            
#ifdef FMT_DEBUG
//            printf("Y_nearR (%d) has %d points.\n", ind_x, y_near.size);
#endif
            int min_parent = -1;
            float cmin = MAXFLOAT;
            if(y_near.size > 0){
                // Find minimum parent:
                for(int i_y = 0; i_y < y_near.size; i_y++){
                    if(ind_x == y_near.indices[i_y]) // Don't connect to self!
                        continue;
                    
                    float c = Costs[y_near.indices[i_y]]+y_near.costs[i_y];
                    if(c < cmin){
                        cmin  = c;
                        min_parent = y_near.indices[i_y];
                    }
                }
            }
            
/*
            int min_parent = -1;
            float c_min = MAXFLOAT;
            // Search each element in Active to see whether ind_x is in it. This is _much_ slower!! Re-implement a backwards neighborhood if this fixes things.
            for(int a = 0; a < Active_size; a++){
                if(Active[a] == ind_x)
                    continue;
                
                int active_ind = find_element(neighborhoodsF[Active[a]].indices, neighborhoodsF[Active[a]].size, ind_x);
                if(active_ind != -1){
                    float c = Costs[Active[a]] + neighborhoodsF[Active[a]].costs[active_ind];
                    if(c < c_min){
                        c_min = c;
                        min_parent = Active[a];
                    }
                }
            }
*/
            if(min_parent != -1){
                if(min_parent != -1 && !collision_inds(ind_x,min_parent, map)){
#ifdef FMT_DEBUG
                    // Make sure the cache is updated appropriately:

#endif
                    // Add parent
                    Parents[ind_x] = min_parent;
                    
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
                    HHeap.push(new_pt);
                }
            }else{
#ifdef FMT_DEBUG
//                printf("y_near is empty (no connection to %d active points)\n", Active_size);
#endif
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
    
    /*
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
     */

#ifdef FMT_DEBUG
    Point tmp;
    printf("Goal point location (preassembly): (%f,%f)\n", points[z].state[0].x,points[z].state[0].y);
    get_poly_der(*(get_neighborhood_poly(Parents[z], z)), (get_neighborhood_poly(Parents[z], z))->duration, &tmp, 0);
    printf("Final point of polynomial (preassembly): (%f,%f)\n", tmp.x,tmp.y);
#endif
    
    curr_ind = z;
    // Insert points backwards from the back
    while(Parents[curr_ind] >= 0){ // Will be false when we reach the start point.
#ifdef FMT_DEBUG
        printf("%d (%f,%f) <- ", curr_ind, points[curr_ind].state[0].x,points[curr_ind].state[0].y);
#endif
        // Error checking:
        if(path_length >= parameters.num_pts+2) // Error in path - too many elements
        {path_length = -1; break;}
        if(curr_ind == Parents[curr_ind]) // Error in path - loop
        {path_length = -1; break;}
        
        // See if we can shortcut from the starting node. If so, skip the rest of the path generation.
        if(Parents[curr_ind] != start_pt_ind){
            float d = smoother->fit_polynomial(&(path[(parameters.num_pts+1)-path_length]), points[start_pt_ind].state, POLY_ORD, points[curr_ind].state, POLY_ORD, map);
            if(d >= 0){
                // Check for collision:
                bool free = true;
                for(int c = 0; c < path[(parameters.num_pts+1)-path_length].num_cells; c++)
                    free = free && map->is_free(path[(parameters.num_pts+1)-path_length].cells[c]);
                
                if(free){
                    printf("shortcut!\n");
                    path_length++;
                    start_ind = start_pt_ind;
                    break;
                }
            }
        }
        
        //Insert appropriate polynomial in appropriate place of path:
        start_ind = Parents[curr_ind];
        copy_polystate(&(path[(parameters.num_pts+1)- path_length]), get_neighborhood_poly(start_ind, curr_ind));
        
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
        // Shift to front - no need to flip.
        int shift_ind = parameters.num_pts+2 - path_length;
        int offset = 1;
        if(start_ind != start_pt_ind){
            copy_polystate(&(path[1]), get_neighborhood_poly(start_ind, curr_ind));
            // Update costs
            if(path[1].cost == -1){
                path[0].cost = -1;
            }else{
                path[0].cost += path[(parameters.num_pts+1)-path_length].cost;
            }
            path_length++;
            offset += 1;
        }

        for(int i = 0; i < path_length; i++){
            copy_polystate(&(path[i+offset]), &(path[shift_ind+i]));
        }
        
    }else{
        // Path generation failed - return badpath.
        path[0].cost = -1;
        return 0;
    }

#ifdef FMT_DEBUG
    printf("Path length: %d, path cost: %f\n", path_length, path[0].cost);
    get_poly_der(path[parameters.num_pts+1], path[parameters.num_pts+1].duration, &tmp, 0);
    printf("Final point of polynomial: (%f,%f)\n", tmp.x,tmp.y);
    get_poly_der(path[path_length], path[path_length].duration, &tmp, 0);
    printf("Final point of flipped polynomial: (%f,%f)\n", tmp.x,tmp.y);
    printf("Goal point location: (%f,%f)\n", points[goal_pt_ind].state[0].x,points[goal_pt_ind].state[0].y);
#endif
    
    // Clean up planner containers
    reset_neighborhood(start_pt_ind);
    reset_neighborhood(goal_pt_ind);
    
#ifdef FMT_WARNING
    // Double check that the path is free:
    for(int seg = 1; seg < path_length; seg++){
        for(int c=0; c < path[seg].num_cells; c++){
            if(!map->is_free(path[seg].cells[c])){
                printf("**********************************************************************************************\n");
                printf("Error: Collision on path being returned. Occurs on segment %d, while planning to (%f,%f,%f)\n", seg,points[goal_pt_ind].state[0].x,points[goal_pt_ind].state[0].y,points[goal_pt_ind].state[0].z);
                printf("**********************************************************************************************\n");
                path[0].cost = -1;
                return 0;
            }
        }
    }
#endif

    return path_length;
}

// Seems like a waste. Should use methods for search.
PolyState* FMT::get_neighborhood_poly(int start_ind, int end_ind){

    // See whether in start_ind's FORWARD neighborhood:
    int ind = find_element(neighborhoodsF[start_ind].indices, neighborhoodsF[start_ind].size, end_ind);
    if(ind != -1)
        return &(neighborhoodsF[start_ind].paths[ind]);
    
    // See whether in end_ind's reverse neighborhood:
    ind = find_element(neighborhoodsR[end_ind].indices, neighborhoodsR[end_ind].size, start_ind);
    if(ind != -1)
        return &(neighborhoodsR[end_ind].paths[ind]);
    
    // Otherwise return bad polynomial
    return NULL;
}


// Returns the forward neighborhood of index filtered by filter
// Since the filter is sorted, can speed up
// Maintains same sorting as original neighborhood
// Assumes that neighborhoods[index].indices and filter_vec are sorted, smallest first.
void FMT::filterF(Neighborhood* filtered_neighborhood, int index, int* filter_vec, int filter_size, int exclude){
    int i_n = 0; int i_f = 0; int ind_n = 0;
    filtered_neighborhood->size = 0;
    float c = 0;
    
    while(i_n < neighborhoodsF[index].size && i_f < filter_size){
        ind_n = neighborhoodsF[index].indices[i_n];
        if(ind_n == filter_vec[i_f]){
            if(ind_n != exclude){
                c = neighborhoodsF[index].costs[i_n];
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

// Returns the reverse neighborhood of index filtered by filter
// Since the filter is sorted, can speed up
// Maintains same sorting as original neighborhood
// Assumes that neighborhoods[index].indices and filter_vec are sorted, smallest first.
void FMT::filterR(Neighborhood* filtered_neighborhood, int index, int* filter_vec, int filter_size, int exclude){
    int i_n = 0; int i_f = 0; int ind_n = 0;
    filtered_neighborhood->size = 0;
    float c = 0;
    
    while(i_n < neighborhoodsR[index].size && i_f < filter_size){
        ind_n = neighborhoodsR[index].indices[i_n];
        if(ind_n == filter_vec[i_f]){
            if(ind_n != exclude){
                c = neighborhoodsR[index].costs[i_n];
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
// Could be made more efficient if we cared to, but this is an infinitesmal piece of the runtime.
bool FMT::reset_neighborhood(int index){
#ifdef FMT_DEBUG
//    printf("reset_neighborhood\n");
#endif
    
    if(index != start_pt_ind && index != goal_pt_ind){
        printf("Error: reset neighborhoods only meant for start and goal points");
        return false;
    }
    
    if(parameters.connection_type == RAD_CON || parameters.connection_type == KNN_CON){
        // clear cache
        // Remove from neighborhood of every neighbor. Because these are added at the end, only need to check last two:
        for(int n = 0; n < parameters.num_pts; n++){
            // check neighbor's end points for target:
            int n_size = neighborhoodsF[n].size;
            
            // Check last item
            if(n_size >= 1 && neighborhoodsF[n].indices[n_size-1] == index){
                //Remove end elements by forgetting they exist:
                neighborhoodsF[n].size -=1;
            }else if(n_size >= 2 && neighborhoodsF[n].indices[n_size-2] == index){
                //Shift end elements and remove end:
                neighborhoodsF[n].indices[n_size-2]=neighborhoodsF[n].indices[n_size-1];
                neighborhoodsF[n].costs[n_size-2]=neighborhoodsF[n].costs[n_size-1];
                neighborhoodsF[n].cache[n_size-2]=neighborhoodsF[n].cache[n_size-1];
                copy_polystate(&(neighborhoodsF[n].paths[n_size-2]),&(neighborhoodsF[n].paths[n_size-1]));
                neighborhoodsF[n].size -=1;
            }
#ifdef FMT_WARNING
            if(find_element(neighborhoodsF[n].indices, neighborhoodsF[n].size, index)!= -1)
                printf("******************************************* Error: neighborhoodF not reset properly\n");
            //            neighborhoods[index].paths[i_n].num_cells = 0;
#endif
        }
        // Reset neighborhood at index
        neighborhoodsF[index].size = 0;
        
        // Now reset reverse neighborhood:
        // Remove from neighborhood of every neighbor. Because these are added at the end, only need to check last two:
        for(int n = 0; n < parameters.num_pts; n++){
            // check neighbor's end points for target:
            int n_size = neighborhoodsR[n].size;
            
            // Check last item
            if(n_size >= 1 && neighborhoodsR[n].indices[n_size-1] == index){
                //Remove end elements by forgetting they exist:
                neighborhoodsR[n].size -=1;
            }else if(n_size >= 2 && neighborhoodsR[n].indices[n_size-2] == index){
                //Shift end elements and remove end:
                neighborhoodsR[n].indices[n_size-2]=neighborhoodsR[n].indices[n_size-1];
                neighborhoodsR[n].costs[n_size-2]=neighborhoodsR[n].costs[n_size-1];
                neighborhoodsR[n].cache[n_size-2]=neighborhoodsR[n].cache[n_size-1];
                copy_polystate(&(neighborhoodsR[n].paths[n_size-2]),&(neighborhoodsR[n].paths[n_size-1]));
                neighborhoodsR[n].size -=1;
            }
#ifdef FMT_WARNING
            if(find_element(neighborhoodsR[n].indices, neighborhoodsR[n].size, index)!= -1)
                printf("******************************************* Error: neighborhoodF not reset properly\n");
            //            neighborhoods[index].paths[i_n].num_cells = 0;
#endif
        }
        // Reset neighborhood at index
        neighborhoodsR[index].size = 0;
        
    }else{
        printf("Error (2) : Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }
    
    return true;
}

// Computes a single neighborhood
bool FMT::compute_neighborhood(int index, Map* map, bool rev){
    if(&(neighborhoodsF[index]) == NULL || &(neighborhoodsR[index])==NULL){
        printf("Call to FMT with uninitialized neighborhood for point %d\n", index);
        return false;
    }

    Neighborhood* nF = neighborhoodsF;
    Neighborhood* nR = neighborhoodsR;
    // In nearly all circumstances, we look forward. The exception is the goal point, or when we specify.
    if( index == goal_pt_ind || rev){
        rev = true;
        nF = neighborhoodsR;
        nR = neighborhoodsF;
    }
    nF[index].size = 0;
    
    if(parameters.connection_type == RAD_CON){
        // Radially connected neighbors
        float dim = 2.0;
        if(parameters.Z_limit > 0)
            dim = 3.0;
        // Coefficient from equation (3) of FMT* paper (Janson, Schmerling, et al.)
        // Not sure what xi is in their equation (tuning parameter?)
        float coeff = (parameters.X_limit*parameters.Y_limit*log(parameters.num_pts))/(dim*parameters.num_pts*1.1);
        parameters.connection_param = 2.2*powf(coeff, 1/dim);
        
        // Compute forward neighborhood, fill in neighbor's reverse neighborhoods
        for(int j = 0; j < parameters.num_pts && nF[index].size < max_neighborhood_size; j++){
            if(j==index)
                continue;
            
            PolyState* p = &(nF[index].paths[nF[index].size]);
            Point tmp_point;
            
            float d = 0;
            
            if(rev || index == goal_pt_ind){
                d = smoother->fit_polynomial(p, points[j].state, POLY_ORD, points[index].state, POLY_ORD, map);
#ifdef FMT_WARNING
                get_poly_der(*p, p->duration, &tmp_point, 0);
                if(d > 0 && dist(tmp_point, points[index].state[0]) > 10.0){
                    
                    printf("WARNING: Polynomial does not reach goal (%f,%f,%f), (%f,%f,%f)!\n",tmp_point.x,tmp_point.y,tmp_point.z, points[index].state[0].x,points[index].state[0].y,points[index].state[0].z);
                }
#endif
            }else{
                d = smoother->fit_polynomial(p, points[index].state, POLY_ORD, points[j].state, POLY_ORD, map);
#ifdef FMT_WARNING
                get_poly_der(*p, p->duration, &tmp_point, 0);
                if(d > 0 && dist(tmp_point, points[j].state[0]) > 10.0){
                    printf("WARNING: Polynomial does not reach goal (%f,%f,%f), (%f,%f,%f)!\n",tmp_point.x,tmp_point.y,tmp_point.z, points[index].state[0].x,points[index].state[0].y,points[index].state[0].z);
                }
#endif
            }
            
            
            if(d < 0)
                continue;
            
            if(true || d < parameters.connection_param)
            {
                // Add to own neighbor list
                if(nF[index].size < max_neighborhood_size){
                    nF[index].indices[nF[index].size] = j;
                    nF[index].costs[nF[index].size] = (nF[index].paths[nF[index].size]).cost;
                    nF[index].cache[nF[index].size] = CACHE_FREE;
                    nF[index].size += 1;
                    
                    // Since r-connected graphs are bidirectional, add to neighbor's neighborhood:
                    if(nR[j].size < max_neighborhood_size){
                        nR[j].indices[nR[j].size] = index;
                        nR[j].costs[nR[j].size]   = nF[index].costs[nF[index].size-1];
                        nR[j].cache[nR[j].size]   = CACHE_FREE;
                        copy_polystate(&(nR[j].paths[nR[j].size]), &(nF[index].paths[nF[index].size-1]));
                        nR[j].size += 1;
                    }
                }

            }
        }
        
    }else if(parameters.connection_type==KNN_CON){
        
        // Placeholder polynomial
        PolyState p;
        Point tmp_point;
        p.cells = new size_t[MAX_POLY_CELLS];
        
        // Set bound to inf
        float bound = MAXFLOAT;
        // Index of the furthest neighbor
        int bound_index = -1;
        
        for(int j = 0; j < parameters.num_pts+2; j++){
            if(index==j) continue;
            if((index!=start_pt_ind || index!= goal_pt_ind) && (j >= parameters.num_pts) ) continue;
                
            
            // If we either have too few neighbors, or point j might be a nearest neighbor:
            // Compute the polynomial:
            float d = 0;
            if(rev || index == goal_pt_ind){
                d = smoother->fit_polynomial(&p, points[j].state, POLY_ORD, points[index].state, POLY_ORD, map);
#ifdef FMT_WARNING
                get_poly_der(p, p.duration, &tmp_point, 0);
                if(d > 0 && dist(tmp_point, points[index].state[0]) > 10.0){
                    printf("Error: Polynomial does not reach goal (%f,%f,%f), (%f,%f,%f)!\n",tmp_point.x,tmp_point.y,tmp_point.z, points[index].state[0].x,points[index].state[0].y,points[index].state[0].z);
                }
#endif
            }else{
                d = smoother->fit_polynomial(&p, points[index].state, POLY_ORD, points[j].state, POLY_ORD, map);
#ifdef FMT_WARNING
                get_poly_der(p, p.duration, &tmp_point, 0);
                if(d > 0 && dist(tmp_point, points[j].state[0]) > 10.0){
                    printf("Error: Polynomial does not reach goal (%f,%f,%f), (%f,%f,%f)!\n",tmp_point.x,tmp_point.y,tmp_point.z, points[j].state[0].x,points[j].state[0].y,points[j].state[0].z);
                }
#endif
            }
            

            if(d < 0)
                continue;
            
            
            if(nF[index].size + 1 < parameters.connection_param || d < bound || (j==goal_pt_ind&&index==start_pt_ind)){
                if(nF[index].size+1 > parameters.connection_param){
                    if(bound_index == -1){
                        printf("Error - trying to remove -1 indexed element. Size: %d, bound %f\n", nF[index].size,bound);
                        continue;
                    }
                    // Remove the furthest neighbor:
                    // Since we want to preserve ordering in terms of index, remove by copying up.
                    
                    for(int k = bound_index; k < nF[index].size; k++){
                        // Index of neighbor:
                        nF[index].indices[k] = nF[index].indices[k+1];
                        // Cost
                        nF[index].costs[k] = nF[index].costs[k+1];
                        // Cache (though all should be unknown anyway)
                        nF[index].cache[k] = nF[index].cache[k+1];
                        // Path
                        if(nF[index].paths[k].cells == NULL)
                            printf("%d, %d has null cells in path\n", index, k);
                        if(nF[index].paths[k+1].cells == NULL)
                            printf("%d, %d has null cells in path\n", index, k+1);
                        
                        copy_polystate(&(nF[index].paths[k]), &(nF[index].paths[k+1]));
                    }
                    // Now forget this element exists
                    nF[index].size -= 1;
                }
                
                // Add to neighbor list
                nF[index].indices[nF[index].size] = j;
                nF[index].costs[nF[index].size] = p.cost;
                nF[index].cache[nF[index].size] = CACHE_UNK;
                
                // Make sure this copies.
                copy_polystate(&(nF[index].paths[nF[index].size]),&p);

                // Add self to neighbor's list if the start or goal point:
                if( (index == start_pt_ind || index == goal_pt_ind) && nR[j].size < max_neighborhood_size){
                    if(find_element(nR[j].indices, nR[j].size, index)==-1){
                        nR[j].indices[nR[j].size] = index;
                        nR[j].costs[nR[j].size] = p.cost;
                        nR[j].cache[nR[j].size] = CACHE_UNK;
                        copy_polystate(&(nR[j].paths[nR[j].size]), &p);
                        // Increment size:
                        nR[j].size += 1;
                    }
                }
                
                // If you did not replace a neighbor, increment size
                if(nF[index].size+1 <= parameters.connection_param){
                    nF[index].size++;
                }
                
                // Now recompute the bound. Could do this somewhat more efficiently, but I don't think it is worth it.
                bound = 0;
                for(int k = 0; k < nF[index].size; k++){
                    if(nF[index].costs[k] >= bound){
                        bound = nF[index].costs[k];
                        bound_index = k;
                    }
                }
            }
        }

        delete [] p.cells;
        

    }else{
        printf("Error (1): Unknown connection type %d!\n", parameters.connection_type);
        return false;
    }

    
    return true;

}

// Computes neighborhoods for sampled points (NOT initial, final point)
bool FMT::initialize_neighborhoods(Map* map_structure){
    
    neighborhoodsF = new Neighborhood[(parameters.num_pts+2)];
    neighborhoodsR = new Neighborhood[(parameters.num_pts+2)];
    
    //Initialize neighborhoods
    for(int i = 0; i < parameters.num_pts+2; i++){
        neighborhoodsF[i].size = 0;
        neighborhoodsF[i].indices = new int[max_neighborhood_size];
        neighborhoodsF[i].costs   = new float[max_neighborhood_size];
        neighborhoodsF[i].cache   = new uint8_t[max_neighborhood_size];
        neighborhoodsF[i].paths   = new PolyState[max_neighborhood_size];
        // Likely unnecessary.
        for(int j = 0; j < max_neighborhood_size; j++){
            neighborhoodsF[i].cache[j] = CACHE_UNK;
            neighborhoodsF[i].paths[j].cells = new size_t[MAX_POLY_CELLS];
        }
        
        neighborhoodsR[i].size = 0;
        neighborhoodsR[i].indices = new int[max_neighborhood_size];
        neighborhoodsR[i].costs   = new float[max_neighborhood_size];
        neighborhoodsR[i].cache   = new uint8_t[max_neighborhood_size];
        neighborhoodsR[i].paths   = new PolyState[max_neighborhood_size];
        // Likely unnecessary.
        for(int j = 0; j < max_neighborhood_size; j++){
            neighborhoodsR[i].cache[j] = CACHE_UNK;
            neighborhoodsR[i].paths[j].cells = new size_t[MAX_POLY_CELLS];
        }
    }
    
    for(int i = 0; i < parameters.num_pts; i++){
        if(!compute_neighborhood(i, map_structure))
            return false;
    }
    
    for(int i = 0; i < parameters.num_pts; i++){
        if(neighborhoodsR[i].size==0){
            if(!compute_neighborhood(i, map_structure,true))
                return false;
        }
    }
    return true;
}


void FMT::clear_cache(){
    for(int i = 0; i < parameters.num_pts; i++){
        for(int j = 0; j < neighborhoodsF[i].size; j++){
            if(neighborhoodsF[i].cache[j] == CACHE_FREE_U || neighborhoodsF[i].cache[j] == CACHE_OCC_U){
                neighborhoodsF[i].cache[j] = CACHE_UNK;
            }
        }
    }
    for(int i = 0; i < parameters.num_pts; i++){
        for(int j = 0; j < neighborhoodsR[i].size; j++){
            if(neighborhoodsR[i].cache[j] == CACHE_FREE_U || neighborhoodsR[i].cache[j] == CACHE_OCC_U){
                neighborhoodsR[i].cache[j] = CACHE_UNK;
            }
        }
    }
}



// Checks for collision using cache to speed up.
// ind1 is assumed to be forward from ind2, meaning ind2 should be in neighborhoodR[ind1], and ind1 in neighborhoodF[ind2]

bool FMT::collision_inds(int ind1, int ind2, Map *map){
    num_total_checks++;
    
    // ====  First find neighbor indices: ==== //
    int neighb_1 = -1; // Index of ind2 in ind1's neighbor list
    int neighb_2 = -1; // Index of ind1 in ind2's neighbor list
    
    neighb_1 = find_element(neighborhoodsR[ind1].indices, neighborhoodsR[ind1].size, ind2);
    neighb_2 = find_element(neighborhoodsF[ind2].indices, neighborhoodsF[ind2].size, ind1);
    Neighborhood* n2 = &(neighborhoodsF[ind2]);
    Neighborhood* n1 = &(neighborhoodsR[ind1]);
    
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
    if( (neighb_1 != -1 && n1->cache[neighb_1] >= CACHE_OCC) || (neighb_2 != -1 && n2->cache[neighb_2] >= CACHE_OCC)){
//        printf("Cached collision\n");
        num_skipped_checks++;
        return true;
    }
    
    // If one has free (we already know neither has collision), then return false
    if( (neighb_1 != -1 && (n1->cache[neighb_1] == CACHE_FREE || n1->cache[neighb_1] == CACHE_FREE_U)) || (neighb_2 != -1 && (n2->cache[neighb_2] == CACHE_FREE || n2->cache[neighb_2] == CACHE_FREE_U))){
        num_skipped_checks++;
        return false;
    }

    // If we're here, then the cache is "CACHE_UNK"
    // Compute the collision value and update the cache
    bool confirmed_value = false;
    
    // Pull out the appropriate polynomial:
    
    PolyState* p;
    if(neighb_1 != -1){
        p = &(n1->paths[neighb_1]);
    }else{
        p = &(n2->paths[neighb_2]);
    }
    
    bool collision_value = false; // Need a collision method which checks indices of the polynomial.
    if(p == nullptr || p->cells == nullptr ){
        collision_value = true;
    }else{
        for(int i = 0; i < p->num_cells; i++){
            if(!(map->is_free(p->cells[i]))){
                collision_value = true;
                confirmed_value = true;
                break;
            }
        }
    }
    
    if(collision_value == true){
        if(confirmed_value){
            if(neighb_1 != -1)
                n1->cache[neighb_1] = CACHE_OCC;
            if(neighb_2 != -1)
                n2->cache[neighb_2] = CACHE_OCC;
        }else{
            if(neighb_1 != -1)
                n1->cache[neighb_1] = CACHE_OCC_U;
            if(neighb_2 != -1)
                n2->cache[neighb_2] = CACHE_OCC_U;
        }
    }else{
        if(confirmed_value){
            if(neighb_1 != -1)
                n1->cache[neighb_1] = CACHE_FREE;
            if(neighb_2 != -1)
                n2->cache[neighb_2] = CACHE_FREE;
        }else{
            if(neighb_1 != -1)
                n1->cache[neighb_1] = CACHE_FREE_U;
            if(neighb_2 != -1)
                n2->cache[neighb_2] = CACHE_FREE_U;
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