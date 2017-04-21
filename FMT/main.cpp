//
//  main.cpp
//  FMT
//
//  Created by Stefan Jorgensen
//  MIT licence
//

#include <stdio.h>
#include <chrono>
#include "utilities.h"
#include "fmt.h"
#include "Explorer.h"

#include "smoother.h"


int main(){
    srand(time(NULL));
    
    
    
//    PathSmoother sm = PathSmoother();
//    sm.test_smoother();
    
    
    Simulator S = Simulator();
    S.run_simulator(5000, 1);
    
    return 0;
}




/* Old tests:
 
 
 
 
 void print_map(FILE* filename, Map* map_vals){
 for(int i = 0; i < map_vals->X_dim; i++){
 for(int j = 0; j < map_vals->Y_dim; j++){
 if((*map_vals)(i,j,0)==1){
 fprintf(filename, "%s","‡ ");
 }else if((*map_vals)(i,j,0) == 0){
 fprintf(filename, "%s","  ");
 }else if((*map_vals)(i,j,0) == 2){
 fprintf(filename, "%s","∫∫");
 }
 }
 fprintf(filename, "\n");
 }
 }
 
 void make_map(Map* map_vals){
 map_vals->X_dim = 400;
 map_vals->Y_dim = 300;
 map_vals->Z_dim = 1;
 map_vals->data = new char[map_vals->X_dim*map_vals->Y_dim*map_vals->Z_dim];
 
 for(int i = 0; i < map_vals->X_dim; i++){
 for(int j = 0; j < map_vals->Y_dim; j++){
 if(map[0][i][j] == 0){
 map_vals->set(i,j,0, MAP_FREE);
 }else if(map[0][i][j] == 1){
 map_vals->set(i,j,0, MAP_OCCU);
 }else{
 map_vals->set(i,j,0, MAP_UNKN);
 }
 }
 }
 FILE* filename = fopen("acsii_plots.txt", "a");
 print_map(filename, map_vals);
 print_map(filename, map_vals);
 fclose(filename);
 }
 
 
 void test_multi(){
 Point start;
 Point goal;
 start.x = 40.0;
 start.y = 40.0;
 start.z = 0.0;
 
 goal.x  = 120;
 goal.y  = 120.0;
 goal.z  = 0.0;
 
 Map map_vals;
 make_map(&map_vals);
 
 FMT F = FMT(200,150,0,1250, RAD_CON, 3);
 int num_success = 0;
 int num_tests  = 200;
 for(int j = 0; j < num_tests; j++){
 printf("Trial %d:\n",j);
 Point* path = F.fmtstar(start, goal, &map_vals);
 if(path[0].x == -1){
 }else{ // Success!
 //            printf("Successful plan! Path has length %f and points\n",path[0].y);
 //            for(int i = 1; i < path[0].x; i++)
 //                printf("(%f,%f,%f)\n", path[i].x,path[i].y,path[i].z);
 num_success++;
 }
 }
 printf("%d/%d success rate\n", num_success, num_tests);
 }
 
 void test_clustering(){
 Map mapvals;
 make_map(&mapvals);
 Explorer exp = Explorer(0, mapvals.X_dim, mapvals.Y_dim);
 exp.force_map_update(&mapvals);
 TimeVar t1 = timeNow();
 exp.cluster_frontiers();
 float time = duration(timeNow()-t1);
 
 
 printf("Time to compute the frontiers is %f ms\n", time/1000000);
 t1 = timeNow();
 exp.compute_costs();
 time = duration(timeNow()-t1);
 printf("Time to compute the costs is %f ms\n", time/1000000);
 //    exp.plot_clusters();
 exp.print_costs();
 
 }
 
 
 // Test passed.
 void test_indexing(){
 Map mapvals;
 make_map(&mapvals);
 Point loc;
 for(int i = 0; i < mapvals.X_dim; i++){
 for(int j = 0; j < mapvals.Y_dim; j++){
 mapvals.num2ind(&loc, mapvals.ind2num(i, j, 0));
 if(i != int(loc.x) || j != int(loc.y) || 0 != int(loc.z)){
 printf("(%d,%d,%d) -> (%f,%f,%f)\n", i,j,0,loc.x,loc.y, loc.z );
 }
 }
 }
 }


*/