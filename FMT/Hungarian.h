//
//  Hungarian.cpp
//  FMT
//
//  Code is or-tools code from Google with minor interface change:
//    github.com/google/or-tools/blob/master/src/algorithms/hungarian.h
//

// Copyright 2010-2014 Google
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// See: //depot/google3/java/com/google/wireless/genie/frontend
//       /mixer/matching/HungarianOptimizer.java

#ifndef __FMT__Hungarian__
#define __FMT__Hungarian__

#include <iostream>
#include <vector>

void assign_max(const std::vector<std::vector<double>>& cost, int* assignments, int num_agents);
void assign_min(const std::vector<std::vector<double>>& cost, int* assignments, int num_agents);

#endif /* defined(__FMT__Hungarian__) */
