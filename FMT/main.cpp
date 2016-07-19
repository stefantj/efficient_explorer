//
//  main.cpp
//  FMT
//
//  Created by Megamind on 7/18/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#include <stdio.h>
#include <chrono>
#include "fmt.h"

#ifndef FMT_TIMING
typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

#endif

int main(){
//    TimeVar t1 = timeNow();
    FMT F = FMT(100,100,0,5000,RAD_CON, 10, 5);
//    double time = duration(timeNow()-t1);
//    std::cout<<"Initialization time: "<< time/1000000000 <<"\n";

    Point start;
    Point goal;
    start.x = 0.0;
    start.y = 0.0;
    start.z = 0.0;
    
    goal.x  = 100.0;
    goal.y  = 100.0;
    goal.z  = 0.0;
    printf("Starting FMT timing test\n");
//    t1 = timeNow();
    std::vector<int> path = F.fmtstar(start,goal);
//    time = duration(timeNow()-t1);
//    std::cout<<"Collision-free time: "<< time/1000000000 <<"\n";
    
    
    return 0;
}