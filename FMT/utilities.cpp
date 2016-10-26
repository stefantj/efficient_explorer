//
//  utilities.cpp
//  FMT
//
//  Created by Stefan Jorgensen
//  MIT licence

#include "utilities.h"
const int MAX_SHUFFLE_SIZE = 4096;
int SHUFFLE_ARRAY[MAX_SHUFFLE_SIZE];

// Helper for nice angles in [0,2PI]
float nice_angle(float theta){
    while(theta < 0)
        theta += M_2_PI;
    while(theta >= M_2_PI)
        theta -= M_2_PI;
    return theta;
}


// Cheap implementation.
// Returns true if a collision found
bool collision(const Point pt1, const Point pt2,  Map* map, bool* confirmed){
    return true;
}

// Finds integer element in a collection.
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

// Removes integer element from a collection
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

// Joins sets. If the sets are sorted smallest to largest, it preserves the sorting.
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

// Helpers for saving data in a julia friendly manner

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

/*
void save_julia_var(FILE* f, std::string var_name, PolyState* variable, int num_vars){
    
    // Print x coefficients:
    fprintf(f, "%s_cx_%d = [", var_name.c_str());
    for(int i = 0; i < )
    
    
    
    
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
*/

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


