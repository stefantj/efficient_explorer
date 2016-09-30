//
//  utilities.cpp
//  FMT
//
//  Created by Megamind on 7/25/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#include "utilities.h"
const int MAX_SHUFFLE_SIZE = 4096;
int SHUFFLE_ARRAY[MAX_SHUFFLE_SIZE];


float nice_angle(float theta){
    while(theta < 0)
        theta += M_2_PI;
    while(theta >= M_2_PI)
        theta -= M_2_PI;
    return theta;
}

// Checks the permutations of the point. Returns true if map is occupied
inline bool check_map(const float* pt, Map* map, bool* confirmed){
    
    if(confirmed != NULL)
        (*confirmed)=true;
    
    if(pt == NULL){
#ifdef FMT_DEBUG
        printf("Error: Null pointer!\n");
#endif
        return true;
    }
    // int casting is the same as floor
    int x=int(pt[0]);
    int y=int(pt[1]);
    int z=int(pt[2]);
#ifdef FMT_DEBUG
    if(x < 0 || y < 0 || z < 0 || x > parameters.X_limit || y > parameters.Y_limit || z > parameters.Z_limit)
    {printf("Error: checking map at invalid index: %f %f %f\n", pt[0], pt[1], pt[2]);
        return false;
    }
#endif
    
    // - - -
    if( (*map)(x,y,z) > MAP_OCC_THRESH){
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&((*map)(x,y,z) == MAP_OCCU);
        }
        return true;
    }else{
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&( ((*map)(x,y,z) <= MAP_FREE));
        }
    }
    // - + -
    if(y+1 < (*map).Y_dim && (*map)(x,y+1,z)>MAP_OCC_THRESH){
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&((*map)(x,y+1,z) == MAP_OCCU);
        }
        return true;
    }else{
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&( ((*map)(x,y+1,z) <= MAP_FREE));
        }
    }
    // + - -
    if(x+1 < (*map).X_dim && (*map)(x+1,y,z)>MAP_OCC_THRESH){
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&((*map)(x+1,y,z) == MAP_OCCU);
        }
        return true;
    }else{
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&( ((*map)(x+1,y,z) <= MAP_FREE));
        }
    }
    // + + -
    if(x+1 < (*map).X_dim && y+1 < (*map).Y_dim && (*map)(x+1,y+1,z)> MAP_OCC_THRESH){
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&((*map)(x+1,y+1,z) == MAP_OCCU);
        }
        return true;
    }else{
        if(confirmed!= NULL){
            (*confirmed) = (*confirmed)&&( ((*map)(x+1,y+1,z) <= MAP_FREE));
        }
    }
/*
    if(z != 0 && z+1 < (*map).Z_dim){ // Fix this warning by simply eliminating this if statement.
        printf("Warning: Collision checker does not expect nonzero z values\n");
        // - - +
        if( (*map)(x,y,z+1) == 1)
            return true;
        // - + +
        if(y+1 < (*map).Y_dim && (*map)(x,y+1,z+1)> MAP_OCC_THRESH)
            return true;
        // + - +
        if(x+1 < (*map).X_dim && (*map)(x+1,y,z+1)> MAP_OCC_THRESH)
            return true;
        // + + +
        if(x+1 < (*map).X_dim && y+1 < (*map).Y_dim && (*map)(x+1,y+1,z+1)> MAP_OCC_THRESH)
            return true;
    }
 */
    return false;
}

// Return true if collision between points
bool collision(const Point pt1, const Point pt2, Map* map, bool* confirmed){
    
    if(dist(pt1,pt2) < 0.001){
        if(confirmed != NULL){
            (*confirmed) = true;
        }
        return false;
    }
    int p1[3];
    int p2[3];
    // Rounded to get onto grid
    p1[0] = int(pt1.x);
    p1[1] = int(pt1.y);
    p1[2] = int(pt1.z);
    
    p2[0] = int(pt2.x);
    p2[1] = int(pt2.y);
    p2[2] = int(pt2.z);
    
#ifdef FMT_DEBUG
    if(p2[2] != 0 || p1[2] != 0)
        printf("Warning: z-dimension ignored for collision checking!\n");
#endif
    if( (*map)(p1[0],p1[1],p1[2]) > MAP_OCC_THRESH || (*map)(p2[0],p2[1],p2[2]) > MAP_OCC_THRESH){
        if(confirmed != NULL){
            (*confirmed) = ((*map)(p1[0],p1[1],p1[2]) == MAP_OCCU || (*map)(p2[0],p2[1],p2[2]) == MAP_OCCU);
        }
        return true;
    }
    
    // Normalized deltas for tracing the line
    float dz = (pt2.z-pt1.z);
    float dy = (pt2.y-pt1.y);
    float dx = (pt2.x-pt1.x);
    
    float scaling = dz;
    if(scaling < 0) scaling=scaling*-1;
    if(dy >= 0 && dy > scaling){
        scaling = dy;
    }else if(-1*dy >= 0 && -1*dy > scaling){
        scaling = -1*dy;
    }
    if(dx >= 0 && dx > scaling){
        scaling = dx;
    }else if(-1*dx >= 0 && -1*dx > scaling){
        scaling = -1*dx;
    }
    
    dz /= scaling;
    dy /= scaling;
    dx /= scaling;
    
    float pt[3];
    
    // Run line starting from pt1:
    pt[0] = pt1.x; pt[1]=pt1.y; pt[2] = pt1.z;
    bool keep_going = false;
    
    while(1){ //TODO: Add 3D element here
        // Move test point along search direction and round
        pt[0] += dx; pt[1] += dy; pt[2] += dz;
        
        // Check whether we've passed point 2:
        if(dy >= 0){
            if(dx >= 0){
                // Traveling in +x +y direction
                keep_going = (pt[0] <= p2[0] && pt[1] <= p2[1]);
            }else{
                // Traveling -x +y direction
                keep_going = (pt[0] >= p2[0] && pt[1] <= p2[1]);
            }
        }else{
            if(dx >= 0){
                // Traveling in the -y +x direction
                keep_going = (pt[0] <= p2[0] && pt[1] >= p2[1]);
            }else{
                // Traveling in the -y -x direction
                keep_going = (pt[0] >= p2[0] && pt[1] >= p2[1]);
                
            }
        }
        if(keep_going == false)
            break;
        
        if(check_map(pt, map, confirmed)){
            return true;
        }
    }
    
    // Run line starting from pt2:
    pt[0] = pt2.x; pt[1]=pt2.y; pt[2] = pt2.z;
    dx = -1*dx; dy = -1*dy; dz = -1*dz;
    keep_going = false;
    
    while(1){ //TODO: Add 3D element here
        // Move test point along search direction and round
        pt[0] += dx; pt[1] += dy; pt[2] += dz;
        
        // Check whether we've passed point 2:
        if(dy >= 0){
            if(dx >= 0){
                // Traveling in +x +y direction
                keep_going = (pt[0] <= p1[0] && pt[1] <= p1[1]);
            }else{
                // Traveling -x +y direction
                keep_going = (pt[0] >= p1[0] && pt[1] <= p1[1]);
            }
        }else{
            if(dx >= 0){
                // Traveling in the -y +x direction
                keep_going = (pt[0] <= p1[0] && pt[1] >= p1[1]);
            }else{
                // Traveling in the -y -x direction
                keep_going = (pt[0] >= p1[0] && pt[1] >= p1[1]);
                
            }
        }
        if(keep_going == false)
            break;
        
        if(check_map(pt, map, confirmed)){
            return true;
        }
    }
    
    return false;
    
}

int find_element(int* collection, int collection_size, int element){

    int curr = collection_size/2;
    // If small collection, linear search is faster
    if(collection_size < 50){
        curr = 0;
        for(int i = 0; i < collection_size; i++){
            if(collection[i]==element){
                return i;
            }
        }
        return -1;
    }
    
    // Otherwise, binary search
    
    
    int LB = 0;
    int UB = collection_size-1;
    bool found = false;
    
    // Check if in range:
    if(collection[LB] > element || collection[UB] < element){
        found = false;
        return collection_size;
    }
    
    
    //Check ends first:
    if(collection[LB] == element){
        found = true; curr = LB;
    }else if(collection[UB] == element){
        found = true; curr = UB; }
    
    while( (!found) && (LB+1 < UB) ){
        // Check current place and cut:
        if(collection[curr] == element){
            found = true;
            break;
        }else if(collection[curr] < element){
            LB = curr + 1;
            curr = (UB-LB)/2+LB;
            // Check LB just in case:
            if(collection[LB] == element){
                found = true; curr = LB; break;}
        }else{
            UB = curr - 1;
            curr = (UB-LB)/2+LB;
            if(collection[UB] == element){
                found = true; curr = UB; break;}
        }
    }
    if(found){
        return curr;
    }else{
        return -1;
    }
    
}

int remove_element(int* collection, int collection_size, int element)
{
    int loc = find_element(collection, collection_size, element);
    
    if(loc!=-1){
        // Now shift out to remove element
        for(int iter = loc+1; iter < collection_size; iter++)
            collection[iter-1] = collection[iter];
            return collection_size-1;
    }else{
        return collection_size;
    }
}

int join_sets(int* set1, int set1_size, int* set2, int set2_size){
    int ind1 = 0; int ind2 = 0;
    int new_set_size = 0;
    // Shuffle together arrays (60% time)
    while(ind1 < set1_size || ind2 < set2_size){
        if(new_set_size >= MAX_SHUFFLE_SIZE){
        }
        if(ind1 < set1_size && ind2 < set2_size){
            if(set1[ind1] <= set2[ind2]){
                SHUFFLE_ARRAY[new_set_size] = set1[ind1];
                new_set_size++;
                if(set1[ind1] == set2[ind2])
                    ind2++;
                ind1++;
            }else{
                SHUFFLE_ARRAY[new_set_size] = set2[ind2];
                new_set_size++;
                ind2++;
            }
        }else if(ind1 < set1_size){
            SHUFFLE_ARRAY[new_set_size] = set1[ind1];
            new_set_size++; ind1++;
        }else if(ind2 < set2_size){
            SHUFFLE_ARRAY[new_set_size] = set2[ind2];
            new_set_size++; ind2++;
        }
    }
    
    // Push data to new set: (30% time)
    if(new_set_size > MAX_SHUFFLE_SIZE){
        printf("Error: merging arrays that are too large!");
    }else{
    for(int i = 0; i < new_set_size; i++)
        set1[i] = SHUFFLE_ARRAY[i];
    }
    
    return new_set_size;
}

// Place the contents of both clusters into the smalleter numbered one.
void merge_clusters(Cluster* C1, Cluster* C2){
    Cluster* target = C1;
    Cluster* other  = C2;
    if(C1->number > C2->number){
        target = C2;
        other =  C1;
    }
    
    // Merge member sets
    target->size = join_sets(target->members, target->size, other->members, other->size);
    
    //Mark other as empty
    other->size = 0;
    
}

void save_julia_var(FILE* f, std::string var_name, float* variable, int num_vars){
    fprintf(f, "%s = [", var_name.c_str());
    for(int i = 0; i < num_vars; i++){
        if(variable[i]!=variable[i]){
            fprintf(f, " NaN ");
        }else{
            fprintf(f, " %f ", variable[i]);
        }
    }
    fprintf(f,"]\n");
    fflush(f);
}

void save_julia_var(FILE* f, std::string var_name, Point* variable, int num_vars){
    fprintf(f, "%s = [", var_name.c_str());
    for(int i = 0; i < num_vars; i++){
        if(variable[i].x != variable[i].x || variable[i].y != variable[i].y || variable[i].z != variable[i].z){
            fprintf(f, "NaN NaN NaN");
        }else{
            fprintf(f, "%f %f %f", variable[i].x, variable[i].y, variable[i].z);
        }
        if(i+1 < num_vars)
            fprintf(f,"; ");
    }
    fprintf(f,"]\n");
    fflush(f);
}

void save_julia_var(FILE* f, std::string var_name, Map* variable){
    fprintf(f, "%s = [", var_name.c_str());
    int code = 0;
    for(int i = 0; i < variable->X_dim; i++){
        for(int j = 0; j < variable->Y_dim; j++){
            if((*variable)(i,j,0) == MAP_UNKN){
                code = fprintf(f, "%.1f ", 0.5);
            }else if((*variable)(i,j,0) >= MAP_OCC_THRESH){
                code = fprintf(f, "%d ", 1);
            }else{
                code = fprintf(f, "%d ",0);
            }
        }
        if(i+1 < variable->X_dim)
            fprintf(f,";\n");
    }
    fprintf(f,"]\n");
    fflush(f);
}


