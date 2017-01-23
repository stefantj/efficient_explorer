//
//  polynomial.h
//  FMT
//
//  Created by Megamind on 10/6/16.
//  MIT Licence

#ifndef __FMT__polynomial__
#define __FMT__polynomial__
#define DOF 0           // Degrees of freedom for optimization

#include <iostream>
#include "utilities.h"

//#define POLY_DEBUG

#define XMIN 0
#define XMAX 400
#define YMIN 0
#define YMAX 300


void print_poly_state(PolyState p);
void print_poly_state_julia(PolyState p, Point* init, int num_init, Point* final, int num_final);
void get_poly_der(PolyState p, double dt, Point* res, int order);


// Holds a polynomial
class PolynomialSmoother{
public:
    PolynomialSmoother();
    
    // Fits a polynomial object between two state space points:
    float fit_polynomial(PolyState* p, Point* init, int num_init, Point* final, int num_final, Map* map);

    // Sets cost coefficients
    void set_q(float* q_coeffs, int num_q);
    void set_kt(float kt);
    
private:
    
    // Linear algebra routines //

    // Forms Q matrix
    void form_Q(double time);
    
    // Computes A^-1B using householder QR decomposition. Not tested for non-square A
    // Result is copied in place to B
    void householder_Ldiv(double A[][MAX_ORDER], int ai_off, int aj_off, int a_sizei, int a_sizej, double B[][MAX_ORDER], int bi_off, int bj_off, int b_sizej);
    
    // Multiplies A*B and saves result in mult_res[][]
    void mat_mult(double A[][MAX_ORDER], int ai_off, int aj_off, int a_sizei,int a_sizej, double B[][MAX_ORDER], int bi_off, int bj_off, int b_sizei, int b_sizej, bool a_trans=false, bool b_trans=false);

    //Degree of polynomial
    int degree;
    
    // Optimization cost weights
    float kT;
    double Q_coeffs[MAX_ORDER];
    
    // Hessian
    double Q_mat[MAX_ORDER][MAX_ORDER];
    // Constraint matrix
    double A_mat[MAX_ORDER][MAX_ORDER];
    double A_inv[MAX_ORDER][MAX_ORDER];
    double mult_res[MAX_ORDER][MAX_ORDER];
};







#endif /* defined(__FMT__polynomial__) */
