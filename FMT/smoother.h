//
//  smoother.h
//  FMT
//
//  Created by Megamind on 8/25/16.
//  Copyright (c) 2016 ASL. All rights reserved.
//

#ifndef __FMT__smoother__
#define __FMT__smoother__
#include "utilities.h"
#include <iostream>

#define MAX_DEGREE 30   // Biggest polynomial we'll try
#define DOF        0    // Probably shouldn't be a #define, but the degrees of freedom for planning

class PathSmoother{
public:
    PathSmoother();

    // Takes in the point path, initial derivatives (and initial bearing?)
    // Returns the degree of the polynomial
    int smooth_path(Point* path, Point* v_init, Point* a_init, Point* v_final, Point* a_final);

    // Compute length of last polynomial path
    float get_length();
    
    // Returns value of polynomial at parameter t
    void get_point(float t, Point* loc_at_t);
    void get_der(float t, int d, Point* vel_at_t);
    
    void test_smoother();
    void test_householder();
    
    // Times each point is visited
    double T_split[MAX_DEGREE];
    int num_points;
    // Degree of path
    int degree;
    double x_coeffs[MAX_DEGREE];
    double y_coeffs[MAX_DEGREE];
    double z_coeffs[MAX_DEGREE];
    double p_coeffs[MAX_DEGREE];
private:
    // Initializes gradient descent variables
    void initialize_gradient_descent();
    
    // Forms hessian given end time
    void form_Q(double time);
    
    // Computes A^-1B using householder QR decomposition. Not tested for non-square A
    // Result is copied in place to B
    void householder_Ldiv(double A[][MAX_DEGREE], int ai_off, int aj_off, int a_sizei, int a_sizej, double B[][MAX_DEGREE], int bi_off, int bj_off, int b_sizej);

    // Multiplies A*B and saves result in mult_res[][]
    void mat_mult(double A[][MAX_DEGREE], int ai_off, int aj_off, int a_sizei,int a_sizej, double B[][MAX_DEGREE], int bi_off, int bj_off, int b_sizei, int b_sizej, bool a_trans=false, bool b_trans=false);
    void slow_transpose(double A[][MAX_DEGREE]);
    
    // Derivatives
    double D_x[MAX_DEGREE];
    double D_y[MAX_DEGREE];
    double D_z[MAX_DEGREE];
    double D_p[MAX_DEGREE];
    
    
    
    // Optimization cost weights
    double Q_coeffs[MAX_DEGREE];
    // Hessian
    double Q_mat[MAX_DEGREE][MAX_DEGREE];
    // Constraint matrix
    double A_mat[MAX_DEGREE][MAX_DEGREE];
    double A_inv[MAX_DEGREE][MAX_DEGREE];
    double mult_res[MAX_DEGREE][MAX_DEGREE];
    
    // Gradient descent variables
    int max_search_steps;
    double step_size;
    double t_diff[MAX_DEGREE];
    double J_curr[MAX_DEGREE];
    double J_last[MAX_DEGREE];
    double best_J;
    double best_x[MAX_DEGREE];
    double best_y[MAX_DEGREE];
    double best_z[MAX_DEGREE];
    double best_p[MAX_DEGREE];
    double best_T[MAX_DEGREE];
    
    
    
    // used for swapping
    double tmp_vec[MAX_DEGREE];
};

#endif /* defined(__FMT__smoother__) */
