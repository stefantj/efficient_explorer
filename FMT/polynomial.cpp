//
//  polynomial.cpp
//  FMT
//
//  Created by Megamind on 10/6/16.
//  MIT Licence
//

#include "polynomial.h"

PolynomialSmoother::PolynomialSmoother(){

    for(int k = 0; k < MAX_ORDER; k++){
        Q_coeffs[k] = 0.0;
    }
    kT = 0;
    
    degree = 0;
    
}


// Fits a polynomials in the state space between two state space points
// Inputs:
// p            Container for the state space polynomial
// init         Array containing initial constraints
// num_init     number of initial constraints
// final        Array containing final constraints
// num_final    number of final constraints
//
// Returns the cost? of the polynomial
float PolynomialSmoother::fit_polynomial(PolyState* p, Point* init, int num_init, Point* final, int num_final, Map* map){
    
#ifdef POLY_DEBUG
//    printf("Planning polynomial from (%.2f,%.2f,%.2f) to (%.2f,%.2f,%.2f)\n", init[0].x,init[0].y,init[0].z, final[0].x,final[0].y,final[0].z);
#endif
    
    if(p == NULL || init == NULL || final == NULL){
        printf("Null pointer\n");
        return -1;
    }

    // For each dimension xyz:

    // Heuristic for initial time:
    double nom_speed = 2.0*36.0/5.0; // TODO: Make this generic (currently in sim-world units)
    p->duration = 1.2*dist(init[0],final[0])/nom_speed; // Nominal distance/nominal velocity
    if(p->duration > 200)
        p->duration = 200;
    
    // Degree of the polynomials:
    this->degree = num_init + num_final + DOF;
    if(this->degree > MAX_ORDER){
        return -1.0;
    }
    p->order = this->degree;

    // Constraint times and orders:
    int constraint_order[degree];
    double constraint_times[degree];
    for(int k=0;k<degree;k++){
        if(k < num_init){
            constraint_order[k] = k;
            constraint_times[k] = 1.0;
        }else if(k==num_init){
            constraint_order[k] = 0;
            constraint_times[k] = 0.0;
        }else{
            constraint_order[k] = constraint_order[k-1]+1;
            constraint_times[k] = 0.0;
        }
    }

    
    for(int j=0; j < degree; j++){
        for(int k=0; k < degree; k++){
            A_mat[j][k] = 0.0;
            A_inv[j][k] = 0.0;
        }
    }
    
    // FORM OPTIMIZATION PROBLEM //
    // Note: Gradient descent loop would start here.
    
    if(DOF > 0){
        // Form Hessian for given time:
        form_Q(p->duration);
    }
    
    // Form A matrix:
    for(int k = 0; k < degree; k++){
        int order = constraint_order[k];
        double time = constraint_times[k]*p->duration;
        
        for(int n= order; n < degree; n++){
            double coeff = 1;
            for(int j =0; j < order; j++){
                coeff*=(n-j);
            }
            for(int pwr = 1; pwr <= n-order; pwr++)
                coeff *= time;
            A_mat[k][n] = coeff;
        }
    }
    
    // Compute inverse of A
    int n_A = num_init+DOF;
    
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
    
    // Compute X coefficients:
    for(int i =0; i < degree; i++){
        p->coefficients_x[i]= 0;
        for(int k = 0; k < num_final; k++)
            p->coefficients_x[i] += A_inv[i][k]*final[k].x;
        for(int k = 0; k < num_init; k++)
            p->coefficients_x[i] += A_inv[i][num_final+k]*init[k].x;
    }
    // Compute Y coefficients:
    for(int i =0; i < degree; i++){
        p->coefficients_y[i]= 0;
        for(int k = 0; k < num_final; k++)
            p->coefficients_y[i] += A_inv[i][k]*final[k].y;
        for(int k = 0; k < num_init; k++)
            p->coefficients_y[i] += A_inv[i][num_final+k]*init[k].y;
    }

    // Compute Z coefficients:
    for(int i =0; i < degree; i++){
        p->coefficients_z[i]= 0;
        for(int k = 0; k < num_final; k++)
            p->coefficients_z[i] += A_inv[i][k]*final[k].z;
        for(int k = 0; k < num_init; k++)
            p->coefficients_z[i] += A_inv[i][num_final+k]*init[k].z;
    }
    
    /*
    //Print A and Ainv to check:
    for(int i = 0; i < degree; i++){
        for(int j = 0; j < degree; j++){
            printf(" %.2f ",A_mat[i][j]);
        }
        printf("\n");
    }
    printf("---\n");
    //Print A and Ainv to check:
    for(int i = 0; i < degree; i++){
        for(int j = 0; j < degree; j++){
            printf(" %.2f ",A_inv[i][j]);
        }
        printf("\n");
    }
    */
    
    
    
    for(int k =0; k < p->order; k++){
        if(p->coefficients_x[k] != p->coefficients_x[k])
            printf("Error: x coeff overflowed for degree %d\n",k);
        if(p->coefficients_y[k] != p->coefficients_y[k])
            printf("Error: y coeff overflowed for degree %d\n",k);
        if(p->coefficients_z[k] != p->coefficients_z[k])
            printf("Error: z coeff overflowed for degree %d\n",k);
    }
    
    
    //   Evaluate the polynomial to get the cells it intersects (replace this with trevor's code):
    double dt = 0.1;
    int num_cells = 0;
    size_t cell_list[MAX_POLY_CELLS]={0};
    for(double t = 0; t < p->duration && num_cells < MAX_POLY_CELLS; t+=dt){
        // Evaluate polynomial at time t
        Point t_loc;
        t_loc.x=p->coefficients_x[0];
        t_loc.y=p->coefficients_y[0];
        t_loc.z=p->coefficients_z[0];


        double t_pow = t;
        for(int k = 1; k < p->order; k++){
            t_loc.x += t*p->coefficients_x[k];
            t_loc.y += t*p->coefficients_y[k];
            t_loc.z += t*p->coefficients_z[k];
            t_pow *= t;
            if(t_pow != t_pow)
                printf("Error: time overflow");
        }

        // Ask map for cell location
        size_t cell_id = map->pt2num(t_loc);
//        printf("t = %.2f: Checking (%.2f,%.2f,%.2f): at %lu\n", t,t_loc.x,t_loc.y,t_loc.z,cell_id);
        
        // Needs to handle error case better.
        
        // Insert if unique:
        bool inserted = false;
        for(int k=0; k < num_cells; k++){
            if(cell_list[k] == cell_id){
                inserted=true;
                break;
            }
        }
        if(!inserted){
            cell_list[num_cells] = cell_id;
            num_cells++;
        }
    }
    
    // If this is a bottleneck, we can do a better sort.
    // Naive seems good enough because there are only a small number of items
    p->num_cells = 0;
    
    while(num_cells > 0 && p->num_cells < MAX_POLY_CELLS){
        // Find smallest item
        int ind_smallest = -1;
        size_t val_smallest = (map->X_dim)*(map->Y_dim)*(map->Z_dim)+1;
        for(int k=0; k < num_cells; k++){
            if(cell_list[k] <= val_smallest){
                ind_smallest = k;
                val_smallest = cell_list[k];
            }
        }
        // Put into p
        if(ind_smallest == -1){

        }else{
            if(ind_smallest >= num_cells || p->num_cells >= MAX_POLY_CELLS){
                printf("Smallest index out of bounds!");
            }else{
                p->cells[p->num_cells] = cell_list[ind_smallest];
                p->num_cells++;
            }
            // Remove item
            cell_list[ind_smallest] = cell_list[num_cells-1];
            num_cells--;
        }
    }
    

    //   Store coefficients, duration, list of unique cells, number of cells and set reverse = false

    
    print_poly_state_julia(*p, init, num_init, final, num_final);

    
    
    return 0.0;
}


// Sets the cost coefficients
void PolynomialSmoother::set_q(float* q_coeffs, int num_q){
    if(num_q >= MAX_ORDER)
        num_q = MAX_ORDER;
    for(int o = 0; o < num_q; o++){
        Q_coeffs[o] = q_coeffs[o];
    }
}

void PolynomialSmoother::set_kt(float kt){
    kT = kt;
}


void PolynomialSmoother::form_Q(double time){
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


void PolynomialSmoother::householder_Ldiv(double A[][MAX_ORDER], int ai_off, int aj_off, int a_sizei, int a_sizej, double B[][MAX_ORDER], int bi_off, int bj_off, int b_sizej){
    double R[a_sizei][a_sizej];
    for(int i = 0; i < a_sizei; i++){
        for(int j = 0; j < a_sizej; j++){
            R[i][j] = A[ai_off+i][aj_off+j];
        }
    }
    
    double v_k[a_sizei];
    
    for(int k = 0; k < a_sizej; k++){
        // Form v_k:
        int v_size = a_sizei-k;
        double normv= 0;
        for(int i=0; i<v_size;i++){
            v_k[i] = R[k+i][k];
            if(i>0) // Skip first element to microoptimize..
                normv += (v_k[i])*(v_k[i]);
        }
        v_k[0] += (1-2*(v_k[0]<0))*sqrtf(normv + v_k[0]*v_k[0]);
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
void PolynomialSmoother::mat_mult(double A[][MAX_ORDER], int ai_off, int aj_off, int a_sizei, int a_sizej, double B[][MAX_ORDER], int bi_off, int bj_off, int b_sizei, int b_sizej, bool a_trans, bool b_trans){
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


// TODO: Test & check
void get_poly_der(PolyState p, double dt, Point* res, int order)
{
    res->x = 0;
    res->y = 0;
    res->z = 0;
    
    if(dt > p.duration)
        return;
    if(p.reverse)
        dt = p.duration-dt;
    
    double t = 1;
    for(int k = order; k < p.order; k++){
        float coeff = 1.0;
        for(int n = order; n < k; n++)
            coeff*= n;
        res->x += coeff*t*p.coefficients_x[k];
        res->y += coeff*t*p.coefficients_y[k];
        res->z += coeff*t*p.coefficients_z[k];
        t*=dt;
    }
}

// Prints for plotting
void print_poly_state_julia(PolyState p, Point* init, int num_init, Point* final, int num_final){
    printf("\n\n");
    printf("x_coeffs = [");
    for(int k= 0; k < p.order-1; k++)
        printf("%f,",p.coefficients_x[k]);
    printf("%f]\n", p.coefficients_x[p.order-1]);

    printf("y_coeffs = [");
    for(int k= 0; k < p.order-1; k++)
        printf("%f,",p.coefficients_y[k]);
    printf("%f]\n", p.coefficients_y[p.order-1]);

    printf("z_coeffs = [");
    for(int k= 0; k < p.order-1; k++)
        printf("%f,",p.coefficients_z[k]);
    printf("%f]\n", p.coefficients_z[p.order-1]);

    printf("t = %f\n", p.duration);
    printf("init = [");
    for(int k = 0; k < num_init-1; k++)
        printf("%f %f %f;",init[k].x,init[k].y,init[k].z);
    printf("%f %f %f]\n", init[num_init-1].x,init[num_init-1].y,init[num_init-1].z);
    printf("final = [");
    for(int k = 0; k < num_final-1; k++)
        printf("%f %f %f;",final[k].x,final[k].y,final[k].z);
    printf("%f %f %f]\n", final[num_final-1].x,final[num_final-1].y,final[num_final-1].z);
    
    // Plotting data
    
    printf("times = linspace(0,t,1000);\n");
    printf("px = zeros(times); py = zeros(times); pz = zeros(times);\n");
    printf("for k=1:size(x_coeffs,1)\n");
    printf("    px += x_coeffs[k]*times.^(k-1)\n");
    printf("    py += y_coeffs[k]*times.^(k-1)\n");
    printf("    pz += z_coeffs[k]*times.^(k-1)\n");
    printf("end\n");
    printf("Using PyPlot; clf();\n");
    printf("plot3D(px,py,pz);\n");
    printf("scatter3D(final[1,1],final[1,2],final[1,3],marker=\"o\")\n");
    printf("scatter3D(init[1,1],init[1,2],init[1,3],marker=\"x\")\n");
}

void print_poly_state(PolyState p){
    printf("duration = %f\n",p.duration);
    printf("num_cells = %d \n",p.num_cells);
    printf("order = %d\n", p.order);
    printf("C_x =[");
    for(int k =0; k < p.order; k++)
        printf("%.2f ", p.coefficients_x[k]);
    printf("]\nC_y =[");
    for(int k =0; k < p.order; k++)
        printf("%.2f ", p.coefficients_y[k]);
    printf("]\nC_z =[");
    for(int k =0; k < p.order; k++)
        printf("%.2f ", p.coefficients_z[k]);
    printf("]\n");
}
