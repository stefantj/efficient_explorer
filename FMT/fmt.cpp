//
//  main.cpp
//  FMT
//
//  Created by Megamind on 7/18/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#include "fmt.h"


#ifdef FMT_TIMING
typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()
#endif


// Explicit constructor method
FMT::FMT(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param, float goal_radius){
    // set parameters
    set_parameters(xlimit, ylimit, zlimit, num_pts, connection_type, connection_param, goal_radius);
    sample_points();
    
    // Compute neighborhoods (return if error)
    if(!compute_neighborhoods())
        return;
#ifdef FMT_DEBUG
    
    printf("compute neighborhoods successful\n");
#endif

    //TODO: Save precomputed neighborhoods to a parameter file
    
    is_initialized = true;
}

void FMT::set_parameters(int xlimit, int ylimit, int zlimit, int num_pts, int connection_type, float connection_param, float goal_radius){
    parameters.X_limit = xlimit;
    parameters.Y_limit = ylimit;
    parameters.Z_limit = zlimit;
    parameters.num_pts = num_pts;
    parameters.connection_type = connection_type;
    parameters.connection_param = connection_param;
    parameters.goal_radius = goal_radius;
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

std::vector<int> FMT::fmtstar(Point start, Point goal){
#ifdef FMT_DEBUG
    printf("fmtstar\n");
    printf("Planning from (%f %f %f) to (%f %f %f)\n", start.x,start.y,start.z,goal.x,goal.y,goal.z);
#endif
    
#ifdef FMT_TIMING
    TimeVar t1 = timeNow();
    double initialization_time = 0;
#endif
    std::vector<int> path;
    //push start, goal onto points list
    points[start_pt_ind].x = start.x;
    points[start_pt_ind].y = start.y;
    points[start_pt_ind].z = start.z;

    points[goal_pt_ind].x  = goal.x;
    points[goal_pt_ind].y  = goal.y;
    points[goal_pt_ind].z  = goal.z;
    
    
/*
    // First, check whether start and goal are trivially connected (check line of sight)
    if( !collision(start, goal)){
        path.push_back(start_pt_ind);
        path.push_back(goal_pt_ind);
        return path;
    }
*/
    // Start by constructing neighborhoods for start/finish points
    if(! compute_neighborhood(start, start_pt_ind))
        return path;
    if(! compute_neighborhood(goal, goal_pt_ind))
        return path;
    
    // Initialize containers for search
    std::vector<int> Vindex (parameters.num_pts+2);
    std::vector<int> Windex (parameters.num_pts+1);
    for(int i = 0; i < parameters.num_pts+1; i++){
        Vindex[i] = i;
        Windex[i] = i;
    }
    Vindex[parameters.num_pts+1] = parameters.num_pts+1;

    // Containers are [index cost] tuples where cost is the shortest path found from start to points[index]
    std::priority_queue<Heap_pt> HHeap;
    int z = start_pt_ind;
    std::vector<int> Active;
    Active.push_back(z);
    std::vector<int> Parents(parameters.num_pts+2);
    std::vector<float> Costs(parameters.num_pts+2,0);

#ifdef FMT_TIMING
    initialization_time = duration(timeNow() - t1);
    TimeVar loop_init = timeNow();
    double filter_time = 0;
    int filter_count = 0;
    int z_count = 0;
#endif
    std::vector<int> H_new;
    while(! is_goal_pt(z)){
#ifdef FMT_TIMING
        z_count++;
#endif
        remove_element(Windex, z);
        H_new.clear();
        Neighborhood x_near = filter(z, Windex, -1);
        Heap_pt new_pt;
        for(std::vector<int>::iterator it_x = x_near.indices.begin(); it_x != x_near.indices.end(); it_x++){
#ifdef FMT_TIMING
            t1 = timeNow();
            filter_count++;
#endif
            Neighborhood y_near = filter(*it_x, Windex, z);
#ifdef FMT_TIMING
            filter_time += duration(timeNow()-t1);
#endif
            if(y_near.indices.size() > 0){
                // Find minimum parent:
                float cmin = MAXFLOAT;
                int y_min = -1;
                for(int i_y = 0; i_y < y_near.size; i_y++){
                    float c = Costs[y_near.indices[i_y]]+y_near.costs[i_y]; //TODO: CHeck this
                    if(c < cmin){
                        cmin  = c;
                        y_min = i_y;
                    }
                }
                if(!collision(points[*it_x], points[y_near.indices[y_min]])){
                    // Add parent
                    Parents[*it_x] = y_near.indices[y_min];
                    // Update cost to go
                    Costs[*it_x] = cmin;
                    // Mark as visited
#ifdef FMT_DEBUG
                    std::cout << "Adding "<< *it_x << " to heap\n";
#endif
                    H_new.push_back(*it_x);
                    // Add to heap
                    new_pt.index = *it_x;
                    new_pt.cost  = cmin;
                    HHeap.push(new_pt); // TODO: CHECK THIS!! Is it copied in, or passed by reference?
                    // Remove from candidates list
                    remove_element(Windex, *it_x);
                }else{
#ifdef FMT_DEBUG
                    std::cout << "Collision found between points (" << points[*it_x].x << ", " << points[*it_x].y << ", " << points[*it_x].z << ") and (";
                    std::cout << points[y_near.indices[y_min]].x << ", " <<  points[y_near.indices[y_min]].y << ", " <<  points[y_near.indices[y_min]].z << ")\n";
#endif
                }
            }
        }
        // Remove z
        printf("Removing %d from "); print_vec(Active);printf("\n");
        Active.erase(std::remove(Active.begin(), Active.end(), z));
        // Update active list
        for(std::vector<int>::iterator it_H = H_new.begin(); it_H != H_new.end(); it_H++)
            Active.push_back(*it_H);
        
        if(HHeap.size() > 0){
            z = HHeap.top().index;
            HHeap.pop();
            printf("Z is "); print_point(points[z]); printf("/n");
            printf("Neighbors are "); print_neighborhood(neighborhoods[z]); printf("/n");

        }else{
            // Try to connect to goal in case this is a sparse sampling issue
            if(!collision(points[z],points[goal_pt_ind])){
                Parents[goal_pt_ind] = z;
                Costs[goal_pt_ind] = Costs[z] + dist(points[z], goal);
                z = goal_pt_ind;
            }else{
                Costs[z] = MAXFLOAT;
            }
            std::cout << "Heap is empty - failure.\n";
            break;
        }
    }
#ifdef FMT_TIMING
    double loop_time = duration(timeNow() - loop_init);
    t1 = timeNow();
#endif
    
    // Construct path by assembling parents from goal to start:
    path.push_back(z);
    while(path[1] != start_pt_ind)
        path.insert(path.begin(), Parents[*path.begin()]);
    
    // Clean up planner containers
    reset_neighborhood(start_pt_ind);
    reset_neighborhood(goal_pt_ind);
    points.pop_back();
    points.pop_back();
/*
    // Smooth the path and tighten cost estimate:
    int i = 1;
    while(i < path.size()-2){
        if(!collision(points[path[i-1]],points[path[i+1]])){
            // remove element i (irrelevant)
            path.erase(path.begin()+i);
        }else{
            i+= 1;
        }
    }
 */
#ifdef FMT_TIMING
    double end_time = duration(timeNow() - t1);
    double total_time = initialization_time+loop_time + end_time;
    
    std::cout<< "FMT timing breakdown: \n" << " Initialization " << initialization_time/1000000 << " ms\n Loop " << loop_time/1000000 << " ms\n Cleanup " << end_time/1000000 << "ms\n";
    std::cout<< "Time spent filtering neighborhoods: " << filter_time/1000000  << " ms ("<< filter_time/total_time<<")\n";
    std::cout << "There were " << filter_count << " calls to filter and "<< z_count << " outer loops \n";
    
    std::cout<< "Number of unexplored sites: " << Windex.size() << "\n";
    
    std::cout << "Finished at goal point (" << points[z].x << ", " << points[z].y << ", " << points[z].z << ")\n";
    
    std::vector<int>::iterator it;
    for(it = path.begin(); it!=path.end(); it++)
        std::cout<< " " << *it;
    
#endif
    
#ifdef FMT_DEBUG
    std::cout<< "Percent of sites explored: " << float(parameters.num_pts - Windex.size())/parameters.num_pts << "\n";
#endif
    
    return path;
}





// Checks whether point is close to goal
bool FMT::is_goal_pt(int index){
    if( dist(points[index], points[goal_pt_ind]) < parameters.goal_radius){
        if(!collision(points[goal_pt_ind], points[index])){
            return true;
        }
    }
    return false;
}

// Return true if collision between points
bool FMT::collision(Point pt1, Point pt2){
    return false;
}

// Removes element from _sorted_ collection. Assumes that first element is smallest
// hopefully this is fast
void FMT::remove_element(std::vector<int> collection, int element){
    std::vector<int>::iterator index = std::lower_bound(collection.begin(), collection.end(), element);
    if(*index == element){
        collection.erase(index);
    }
}

// Returns the neighborhood of index filtered by filter
// Since the filter is sorted, can speed up
// Maintains same sorting as original neighborhood
// Assumes that neighborhoods[index].indices and filter_vec are sorted, smallest first.
Neighborhood FMT::filter(int index, std::vector<int> filter_vec, int exclude){
    Neighborhood filtered_neighborhood;
    filtered_neighborhood.size = 0;
    
    int L1 = neighborhoods[index].size;
    int L2 = filter_vec.size();
    
    int i_n = 0;
    int i_f = 0;
    int ind_n = 0;
    float c = 0;
    while(i_n < L1 && i_f < L2){
        ind_n = neighborhoods[index].indices[i_n];
        if(ind_n == filter_vec[i_f]){
            c = neighborhoods[index].costs[i_n];
            filtered_neighborhood.indices.push_back(ind_n);
            filtered_neighborhood.costs.push_back(c);
            filtered_neighborhood.size+=1;
            i_n+=1; i_f += 1;
        }else if( ind_n < filter_vec[i_f]){
            i_n += 1;
        }else{
            i_f += 1;
        }
    }
    
    return filtered_neighborhood;
}



// Remove neighbor from neighborhood
void FMT::remove_neighbor(int index, int neighbor_index){
    for(int i = 0; i < neighborhoods[index].size; i++){
        if(neighborhoods[index].indices[i] == neighbor_index){
            neighborhoods[index].indices.erase(neighborhoods[index].indices.begin()+i);
            neighborhoods[index].costs.erase(neighborhoods[index].costs.begin()+i);
            neighborhoods[index].size -= 1;
            return;
        }
    }
}


/*** Private methods ***/

void FMT::sample_points(){
#ifdef FMT_DEBUG
    printf("Sampling points\n");
#endif
    
    //Seed RNG
    srand(time(NULL));
    
    // Preallocate vectors (+2 because we add start and end points at the end)
    points.reserve(parameters.num_pts+2);
    neighborhoods.reserve(parameters.num_pts+2);
    start_pt_ind = parameters.num_pts;
    goal_pt_ind  = parameters.num_pts+1;
    
    // Generate random numbers
    float x_scaling = float(parameters.X_limit)/RAND_MAX;
    float y_scaling = float(parameters.Y_limit)/RAND_MAX;
    printf("x scaling is %f y scaling is %f\n", x_scaling,y_scaling);
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
bool FMT::compute_neighborhood(Point pt, int index){
    
    if(!is_initialized){
        printf("Call to compute_neighborhood before initialization is complete. \n");
        return false;
    }
    
    neighborhoods[index].size = 0;
    
    if(parameters.connection_type == RAD_CON){
        // Radially connected neighbors
        for(int j = 1; j < parameters.num_pts; j++){
            if(dist(points[index],points[j]) < parameters.connection_param)
            {
                // Add to neighbor list
                neighborhoods[index].indices.push_back(j);
                neighborhoods[index].costs.push_back(dist(points[index],points[j]));
                neighborhoods[index].size += 1;
                
                // Since r-connected graphs are bidirectional, add to other neighborhood
                neighborhoods[j].indices.push_back(index);
                neighborhoods[j].costs.push_back(dist(points[index],points[j]));
                neighborhoods[j].size += 1;
            }
        }
#ifdef FMT_DEBUG
        printf("neighborhood has %d elements\n", neighborhoods[index].size);
#endif
        
        // TODO: Sort the neighborhood by distance (? not sure if it will improve performance at all)
        return true;
    }else if(parameters.connection_type == KNN_CON){
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
    }else{
        printf("Error: Unknown connection type!\n");
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
        // Remove from neighborhood of every neighbor:
        for(std::vector<int>::iterator it_n = neighborhoods[index].indices.begin(); it_n != neighborhoods[index].indices.end(); it_n++){
            if(*it_n != index){
                remove_neighbor(*it_n, index);
            }
        }
        // Reset neighborhood at index
        neighborhoods[index].size = 0;
        neighborhoods[index].indices.clear();
        neighborhoods[index].costs.clear();

        
    }else if(parameters.connection_type == KNN_CON){
        printf("Error: KNN not fully implemented yet!\n");
        return false;
    }else{
        printf("Error: Unknown connection type!\n");
        return false;
    }
    
    return true;
}


// Computes neighborhoods for sampled points (NOT initial, final point)
bool FMT::compute_neighborhoods(){

    //Initialize neighborhood sizes
    for(int i = 0; i < parameters.num_pts+2; i++)
        neighborhoods[i].size = 0;
    
    if(parameters.connection_type == RAD_CON){
        // Radially connected neighbors
        for(int i = 0; i < parameters.num_pts; i++){
            for(int j = i+1; j < parameters.num_pts; j++){
                if(dist(points[i],points[j]) < parameters.connection_param)
                {
                    // Add to neighbor list
                    neighborhoods[i].indices.push_back(j);
                    neighborhoods[i].costs.push_back(dist(points[i],points[j]));
                    neighborhoods[i].size += 1;
                    
                    // Since r-connected graphs are bidirectional, add to other neighborhood
                    neighborhoods[j].indices.push_back(i);
                    neighborhoods[j].costs.push_back(dist(points[i],points[j]));
                    neighborhoods[j].size += 1;
                }
            }
        }
        // TODO: Sort the neighborhood by distance (? not sure if it will improve performance at all)
        return true;
    }else if(parameters.connection_type == KNN_CON){
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
        return true;
    }else{
        printf("Error: Unknown connection type!\n");
        return false;
    }
}

void FMT::print_point(Point pt){
    printf("(%f, %f, %f)",pt.x,pt.y,pt.z);
}

void FMT::print_neighborhood(Neighborhood n){
    std::cout<<"[ " << n.indices[0];
    for(int i = 1; i < n.size; i++){
        printf(", %d",n.indices[i]);
    }
    printf("]");
}

void FMT::print_vec(std::vector<int> v){
    if(v.size() == 0){
        printf("<>");
        return;
    }
    std::cout<<"< " << v[0];
    for(int i = 1; i < v.size(); i++){
        printf(", %d",v[i]);
    }
    printf(">");
}

