//
//  smoother.cpp
//  FMT
//
//  Created by Stefan Jorgensen
//  MIT licence
//

#include "smoother.h"

PathSmoother::PathSmoother(){

//    printf("Warning: Path smoother only plans in 3d space (skips yaw).\n");
    
    
    // Cost coefficients. Index corresponds to order, and number corresponds to priority.
    Q_coeffs[0] = 0.01;
    Q_coeffs[1] = 0.01;
    Q_coeffs[2] = 0.01;
    Q_coeffs[3] = 1;
    Q_coeffs[4] = 0.01;
    
    // Number of search steps for gradient descent
    max_search_steps = 1;
    step_size = 0.1;
    
    num_points = 0;
    
    
}

void PathSmoother::initialize_gradient_descent(){
    for(int k=0; k < MAX_DEGREE; k++){
        t_diff[k] = 0.0;
        J_curr[k] = 0.0;
        J_last[k] = 0.0;
        
        best_x[k] = 0.0;
        best_y[k] = 0.0;
        best_z[k] = 0.0;
        best_p[k] = 0.0;

        best_T[k] = 0.0;
    }
    best_J = 10000000; // Should replace this by MAX_double or something.
}


// Just a single segment
int PathSmoother::smooth_path(Point* path, Point* v_init, Point* a_init, Point* v_final, Point* a_final){
//    TimeVar T;
//    T = timeNow();
//    double clock;
    // Extract information about problem size:
    num_points = (int)path[0].x;            // Path metadata is stored in first element.
    int constr_dim = 4;                     // Assume both initial and final derivatives (vel/acc) are constrained
    if(v_final == NULL)
        constr_dim = 2;                     // Only initial derivatives are constrained.
    degree = num_points + constr_dim + DOF; // Degree is number of constraints + DOF

    // Initialize containers for tracking best path
    initialize_gradient_descent();
 
    // Initialize coefficients/derivatives
    for(int k=0; k < MAX_DEGREE; k++){
        x_coeffs[k] = 0.0;
        y_coeffs[k] = 0.0;
        z_coeffs[k] = 0.0;
        p_coeffs[k] = 0.0;
        
        D_x[k] = 0.0;
        D_y[k] = 0.0;
        D_z[k] = 0.0;
        D_p[k] = 0.0;
        
        T_split[k] = 0.0;
    }
    
    // Form constraint vectors
    // Point constraints (points 2:end)
    int con_ind = 0;
    for(con_ind = 0; con_ind < num_points-1; con_ind++){
        D_x[con_ind] = path[con_ind+2].x; // path points are indexed starting from 1
        D_y[con_ind] = path[con_ind+2].y; // path points are indexed starting from 1
        D_z[con_ind] = path[con_ind+2].z; // path points are indexed starting from 1
        D_p[con_ind] = 0.0;
    }
    if(v_final != NULL){
        D_x[con_ind] = v_final->x;
        D_y[con_ind] = v_final->y;
        D_z[con_ind] = v_final->z;
        D_p[con_ind] = 0;
        con_ind+=1;
        D_x[con_ind] = a_final->x;
        D_y[con_ind] = a_final->y;
        D_z[con_ind] = a_final->z;
        D_p[con_ind] = 0;
    }
    // Add in initial constraints:
    con_ind+=1;
    D_x[con_ind] = path[1].x;
    D_y[con_ind] = path[1].y;
    D_z[con_ind] = path[1].z;
    D_p[con_ind] = 0;
    
    con_ind+=1;
    D_x[con_ind] = v_init->x;
    D_y[con_ind] = v_init->y;
    D_z[con_ind] = v_init->z;
    D_p[con_ind] = 0;
    
    con_ind+=1;
    D_x[con_ind] = a_init->x;
    D_y[con_ind] = a_init->y;
    D_z[con_ind] = a_init->z;
    D_p[con_ind] = 0;
    
    
    // Initial time hueristic:
    double nom_speed = 2.0*36.0/5.0; // TODO: Make this generic
    double T_nom = 1.2*(path[0].y)/nom_speed; // Nominal distance/nominal velocity
//    T_nom = 15.0; // can't be too large or bad stuff

    if(T_nom > 200)
        T_nom = 200;

//    printf("T_nom = %f\n", T_nom);
    
//    printf("T_split = ");
    for(int k = 0; k < num_points-1; k++){
        T_split[k+1] = T_nom/((float)num_points-1) + T_split[k];
//        printf("%f ", T_split[k+1]);
    }
//    printf("\n");
    
    for(int j=0; j < degree; j++){
        for(int k=0; k < degree; k++){
            A_mat[j][k] = 0.0;
            A_inv[j][k] = 0.0;
//            if(DOF > 0){
//                Q_mat[j][k] = 0.0;
//            }
        }
    }
//    clock = duration(timeNow() - T);
//    printf("T_init: %f\n", clock/1000.0);
//    T = timeNow();
    

    // END INITIALIZATION PHASE //

        if(DOF > 0){
        // Form Hessian for given time:
            form_Q(T_split[num_points-1]);
        }
    
    
    // Form A matrix: //
        double constraint_times[degree];
        int constraint_order[degree];
//        con_ind = 0;
// TODO: Double check this
        for(int k =0; k < num_points-1; k++){
            constraint_times[k] = T_split[k+1];
            constraint_order[k] = 0;
        }
        con_ind= num_points-1;
        if(v_final != NULL){ // Add final derivative constraints
            constraint_times[con_ind] = T_split[num_points-1];
            constraint_order[con_ind] = 1;
            con_ind+=1;
            constraint_times[con_ind] = T_split[num_points-1];
            constraint_order[con_ind] = 2;
            con_ind+=1;
        }
        // Initial constraints:
        for(int k = 0; k <= DOF+2; k++){
            constraint_order[con_ind+k] = k;
            constraint_times[con_ind+k] = 0.0;
        }
        
        for(int k=0; k < degree; k++){
            if(constraint_order[k] >= degree)
                continue;
            
            int order = constraint_order[k];
            double time = constraint_times[k];
            
            for(int n= order; n < degree; n++){
                double coeff = 1;
                for(int j =0; j < order; j++){
                    coeff*=(n-j);
                }
                for(int pwr = 1; pwr <= n-order; pwr++)
                    coeff *= time;
                A_mat[k][n] = coeff;
            }
        } // End A matrix formation
//        clock = duration(timeNow() - T);
//        printf("T_A_mat: %f\n", clock/1000.0);
//        T = timeNow();

        // Compute inverse of A
        int n_A = 3+DOF;
//        if(v_final==NULL)
//            n_A -= 2;
        
        // Compute A_inv using block matrix identity:
        // Amat = [A B     Amat_inv = [0    Cinv
        //         C 0],               Binv -Binv*A*Cinv

        // Useful shortcuts for referencing the blocks of A_mat:
        int aB_ioff = 0;
        int aB_joff = n_A;
        int aC_ioff = degree-n_A;
        int aC_joff = 0;

        // Useful shortcuts for referencing the blocks of Amat_inv:
        int aIB_ioff = 0;
        int aIB_joff = degree-n_A;
        int aIC_ioff = n_A;
        int aIC_joff = 0;
        int aID_ioff = n_A;
        int aID_joff = degree-n_A;
        
        // A_inv 1,1 block is already zero
        // Fill in A_inv 1,2 block:
        // A_c is diagonal with size (degree-n_A)
        for(int k = 0; k < n_A; k++){
            A_inv[aIB_ioff+k][aIB_joff+ k] = 1/(A_mat[aC_ioff+k][aC_joff+k]);
        }
        // Fill in A_inv 2,1 block:
        // Compute B_inverse
        for(int k =0; k < degree-n_A; k++)
            A_inv[aIC_ioff+k][aIC_joff+k] = 1.0;
        // Put A_B\I into 2,1 block
        householder_Ldiv(A_mat, aB_ioff, aB_joff, degree-n_A, degree-n_A, A_inv, aIC_ioff, aIC_joff, degree-n_A);

        // Fill in A_inv 2,2 block:
        // Binv*A
        mat_mult(A_inv, aIC_ioff, aIC_joff, degree-n_A, degree-n_A, A_mat, 0, 0, degree-n_A, n_A);
        
        // -(Binv*A)*Cinv
        for(int i = 0; i < degree-n_A; i++)
            for(int j =0; j < n_A; j++)
                A_inv[aID_ioff+i][aID_joff+j] = -mult_res[i][j]/A_mat[aC_ioff+j][aC_joff+j];
        
//        clock = duration(timeNow() - T);
//        printf("T_A_inv: %f\n", clock/1000.0);
//        T = timeNow();

        // IF there are degrees of freedom, solve optimization problem: (untested)
        if(false && DOF > 0){
            
            // Compute Ainv'QA_inv
            mat_mult(Q_mat, 0, 0, degree, degree, A_inv, 0, 0, degree, degree);
            double R_mat[MAX_DEGREE][MAX_DEGREE];
            for(int i =0; i < degree; i++)
                for(int j =0; j < degree; j++)
                    R_mat[i][j] = mult_res[i][j];
            mat_mult(A_inv, 0, 0, degree, degree, R_mat, 0, 0, degree, degree, true);
            for(int i =0; i < degree; i++)
                for(int j =0; j < degree; j++)
                    R_mat[i][j] = mult_res[i][j];
            int mod = 4;
            if(v_final==NULL)
                mod = 2;
            
            double opt_mat[MAX_DEGREE][MAX_DEGREE];
            // Fill in -Rfp' to optmat
            for(int i = 0; i < degree-(num_points+mod); i++)
                for(int j=0; j < num_points+mod; j++)
                    opt_mat[i][j] = -R_mat[j][(num_points+mod)+i];
            // Compute -Rfp'\ Rpp
            householder_Ldiv(R_mat, num_points+mod, num_points+mod, degree, degree, opt_mat, 0, 0,num_points+mod );

            // Fill in free variables:
            for(int i = 0; i < (degree-num_points-mod); i++){
                for(int k=0; k < num_points+mod; k++){
                   D_x[num_points+mod+i] += opt_mat[i][k]*D_x[k];
                   D_y[num_points+mod+i] += opt_mat[i][k]*D_y[k];
                   D_z[num_points+mod+i] += opt_mat[i][k]*D_z[k];
                   D_p[num_points+mod+i] += opt_mat[i][k]*D_p[k];
                }
            }
        } // End DOF calculations
        
        // Fill in polynomials:
        for(int i =0; i < degree; i++){
            for(int k = 0; k < degree; k++){
                x_coeffs[i] += A_inv[i][k]*D_x[k];
                y_coeffs[i] += A_inv[i][k]*D_y[k];
                z_coeffs[i] += A_inv[i][k]*D_z[k];
                p_coeffs[i] += A_inv[i][k]*D_p[k];
            }
        }
        
        // Now we compute the gradient:
        
        
//        clock = duration(timeNow() - T);
//        printf("T_end: %f\n", clock/1000.0);

    
    return degree;
}


float PathSmoother::get_length(){
    int num_t_steps = 100;
    float dt = T_split[num_points-1]/((float)num_t_steps);
    float length = 0;
    float t = 0;
    // This is a very bad way of integrating but is fine for our purposes.
    for(int n = 0; n < num_t_steps; n++){
        float x_val = 0;
        float y_val = 0;
        float z_val = 0;
        t += dt;
        float t_pow = 1;
        // First derivative
        for(int k = 1; k < degree; k++){
            x_val += x_coeffs[k]*(k)*t_pow;
            y_val += y_coeffs[k]*(k)*t_pow;
            z_val += z_coeffs[k]*(k)*t_pow;
            t_pow *= t;
        }
        length += sqrtf( x_val*x_val + y_val*y_val + z_val*z_val )*dt;
    }
    return length;
}


void PathSmoother::form_Q(double time){
    for(int k = 0; k < degree; k++){
        if(Q_coeffs[k] != 0){
            for(int i = 0; i < degree; i++){
                for(int l = 0; l < degree; l++){
                    if( (i >= l || i+1 >= k ) && l+1 >= k ){
                        double c_k = Q_coeffs[k]*2;
                        for(int m = 0; m < k; m++)
                            c_k *= (i+1-m)*(l+1-m);
                        Q_mat[i][l] += c_k * (powf(time, i+l+3-2*k)/(i+l+3-2*k));
                    }
                }
            }
        }
    }
}

void PathSmoother::test_smoother(){
    // Make a point path:
    int num_pts=4;
    Point pp[num_pts+1];
    pp[0].x = double(num_pts);
    
    pp[1].x = 0.0;
    pp[1].y = 0.0;
    pp[1].z = 0.0;

    pp[2].x = 1.0;
    pp[2].y = 5.0;
    pp[2].z = 0.0;

    pp[3].x = 2.0;
    pp[3].y = 5.0;
    pp[3].z = 0.0;
    
    pp[4].x = 3.0;
    pp[4].y = 5.0;
    pp[4].z = 0.0;

    
    Point v_i;
    v_i.x = 0.0; v_i.y = 0.0; v_i.z = 0.0;
    Point v_f;
    v_f.x = 0.0; v_f.y = 0.0; v_f.z = 0.0;
    Point a_i;
    a_i.x = 0.0; a_i.y = 0.0; a_i.z = 0.0;
    Point a_f;
    a_f.x = 0.0; a_f.y = 0.0; a_f.z = 0.0;
    
    double time;
    TimeVar T = timeNow();
    smooth_path(pp, &v_i, &a_i, &v_f, &a_f);
    time = duration(timeNow()-T);
    printf("time is %f microseconds\n", time/1000.0);
    // Print output to make julia plot:
    printf("coeff_x=[");
    for(int i =0; i < degree; i++){
        printf("%f", x_coeffs[i]);
        if(i+1==degree){
            printf("]\n");
        }else{
            printf("; ");
        }
    }
    printf("coeff_y=[");
    for(int i =0; i < degree; i++){
        printf("%f", y_coeffs[i]);
        if(i+1==degree){
            printf("]\n");
        }else{
            printf("; ");
        }
    }
    printf("path=[");
    for(int i =0; i < num_pts; i++){
        printf("%f %f 0.0 0.0", pp[i+1].x, pp[i+1].y);
        if(i+1==num_pts){
            printf("]\n");
        }else{
            printf("; ");
        }
    }
    printf("plot_poly(coeff_x,coeff_y,path,%f,3)", T_split[num_pts-1]);
}



void PathSmoother::test_householder(){
    // Form a random matrix
    double A[MAX_DEGREE][MAX_DEGREE];
    double B[MAX_DEGREE][MAX_DEGREE];
    double B_copy[MAX_DEGREE][MAX_DEGREE];
    for(int k =0; k < MAX_DEGREE; k++){
        for(int j=0; j<MAX_DEGREE; j++){
            B[k][j] = k+j*j;
            B_copy[k][j] = B[k][j];

            if(k==j){
                A[k][j] = 1.0;
            }else{
                A[k][j] = 0.0;
            }
        }
    }
    printf("A is identity, B is random. Try Ainv B.");
    householder_Ldiv(A, 0, 0, MAX_DEGREE-20, MAX_DEGREE-20, B, 0, 0, MAX_DEGREE-20);
    double error = 0;
    
    for(int k =0; k < MAX_DEGREE; k++){
        for(int j=0; j<MAX_DEGREE; j++){
            error += abs(B[k][j]-B_copy[k][j]);
        }
    }
    printf(" Error is: %f\n", error );
}


void PathSmoother::householder_Ldiv(double A[][MAX_DEGREE], int ai_off, int aj_off, int a_sizei, int a_sizej, double B[][MAX_DEGREE], int bi_off, int bj_off, int b_sizej){
    double R[a_sizei][a_sizej];
    for(int i = 0; i < a_sizei; i++){
        for(int j = 0; j < a_sizej; j++){
            R[i][j] = A[ai_off+i][aj_off+j];
        }
    }
    
    double v_k[a_sizei];
    for(int k = 0; k < a_sizei; k++)
        v_k[k] = 0.0;
    
    for(int k = 0; k < a_sizej; k++){
        // Form v_k:
        int v_size = a_sizei-k;
        double normv= 0;
        for(int i=0; i<v_size;i++){
            v_k[i] = R[k+i][k];
            if(i>0) // Skip first element to microoptimize..
                normv += (v_k[i])*(v_k[i]);
        }
//        v_k[0] += (1-2*(v_k[0]<0))*sqrtf(normv + v_k[0]*v_k[0]);

        v_k[0] += sign(v_k[0])*sqrtf(normv + v_k[0]*v_k[0]);
        normv = sqrtf(normv + v_k[0]*v_k[0]);
        for(int i =0; i < v_size; i++)
            v_k[i] /= normv;
        
        // Update R:
        for(int col=k; col < a_sizej; col++){
            // Compute v_k'R[:,col]
            double rowsum = 0;
            for(int c = 0; c < v_size; c++){
                rowsum += v_k[c]*R[k+c][col];
            }
            // Update R
            for(int r=0; r < v_size; r++){
                R[r+k][col] -= 2*v_k[r]*rowsum;
            }
        }
        
        // Update Q'B in place
        for(int col = 0; col < b_sizej; col++){
            // Compute V_k'B[:,col]
            double rowsum = 0;
            for(int c = 0; c < v_size; c++)
                rowsum += v_k[c]*B[k+c+bi_off][col+bj_off];
            for(int r = 0; r < v_size; r++)
                B[bi_off+k+r][bj_off+col] -= 2*rowsum*v_k[r];

        }
    }
    
    // B now holds Q'*B. We solve R\(Q'B) using backsubstitution:
    for(int col = 0; col < a_sizej; col++){
        for(int row = a_sizej-1; row >=0; row--){
            double rowsum = 0;
            for(int k = a_sizei-1; k > row; k--){
                rowsum += B[bi_off+k][bj_off+col]*R[row][k];
            }
            B[bi_off+row][bj_off+col] = (B[bi_off+row][bj_off+col]-rowsum)/R[row][row];
        }
    }
    
}

// Multiplies A*B and saves result in mult_res[][]
void PathSmoother::mat_mult(double A[][MAX_DEGREE], int ai_off, int aj_off, int a_sizei, int a_sizej, double B[][MAX_DEGREE], int bi_off, int bj_off, int b_sizei, int b_sizej, bool a_trans, bool b_trans){
    int size1 = a_sizei;
    int size_common = a_sizej;
    int size2 = b_sizej;
    int b_inside = b_sizei;
    if(a_trans){
        size1 = a_sizej;
        size_common = a_sizei;
    }
    
    if(b_trans){
        size2 = b_sizei;
        b_inside = b_sizej;
    }
    if(b_inside != size_common){
        printf("ERROR: Matrix dimensions are not compatible! (%dx%d)(%dx%d)\n", size1, size_common, b_inside, size2);
    }
    for(int i = 0; i < size1; i++){
        for(int j =0; j < size2; j++){
            double sum = 0;
            for(int k=0; k < size_common; k++){
                int a_i = i+ai_off;
                int a_j = k+aj_off;
                int b_i = k+bi_off;
                int b_j = j+bj_off;
                if(a_trans){
                    a_i = k + ai_off; a_j = i+aj_off;
                }
                if(b_trans){
                    b_i = j+bi_off; b_j = k+bj_off;
                }
                sum += A[a_i][a_j]*B[b_i][b_j];
            }
            mult_res[i][j] = sum;
        }
    }
}

void PathSmoother::get_point(float t, Point* loc_at_t){
    loc_at_t->x = 0.0;
    loc_at_t->y = 0.0;
    loc_at_t->z = 0.0;
    
    double coeff = 1;
    for(int k = 0; k < degree; k++){
        loc_at_t->x += x_coeffs[k]*coeff;
        loc_at_t->y += y_coeffs[k]*coeff;
        loc_at_t->z += z_coeffs[k]*coeff;
        coeff *= t;
    }    
}

void PathSmoother::get_der(float t, int d, Point* vel_at_t){
    vel_at_t->x = 0.0;
    vel_at_t->y = 0.0;
    vel_at_t->z = 0.0;
    
    double coeff = 1;
    for(int k = d; k < degree; k++){
        vel_at_t->x += x_coeffs[k]*coeff;
        vel_at_t->y += y_coeffs[k]*coeff;
        vel_at_t->z += z_coeffs[k]*coeff;
        coeff *= t;
    }
}
