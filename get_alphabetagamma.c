#include "header.h"



void get_union(int n, double * a, double * b, int * len_a_union_b, double * a_union_b){

    int counter =0;

    int a_counter = 0;
    int b_counter = 0;


    if(a[0]<b[0]){
        a_union_b[counter++] = a[0];
        a_counter++;
    }else{
        a_union_b[counter++] = b[0];
        b_counter++;
    }

    while(( a_counter <n )|| ( b_counter < n)){
        if(a_counter>=n){
            if(a_union_b[counter-1] < b[b_counter]){
                a_union_b[counter++] = b[b_counter];
            }
            b_counter++;

        }
        else if(b_counter>=n){
            if(a_union_b[counter-1] < a[b_counter]){
                a_union_b[counter++] = a[a_counter];
            }
            a_counter++;

        }
        else if (a[a_counter]<b[b_counter]){
            if(a_union_b[counter-1] < a[a_counter]){
                a_union_b[counter++] = a[a_counter];
            }
            a_counter++;
            
        }
        else if(b[b_counter]<a[a_counter]){
            if(a_union_b[counter-1] < b[b_counter])
            {
                a_union_b[counter++] = b[b_counter];
            }
            b_counter++;
        }

        else{
            if (a_union_b[counter-1] < b[b_counter])
            {
                a_union_b[counter++] = b[b_counter];
            }
            a_counter++;
            
        }
    }
    (*len_a_union_b) = counter;

}
void get_alphabetagamma(double * x, int n, double * a, double * b, double * a_union_b, int len_a_union_b, int * alpha, 
    int * beta, int * gamma, int * k, int * m, int debug){


    // initialize alpha beta gamma to zero
    memset(alpha, 0, sizeof(int)*n);
    memset(beta, 0, sizeof(int)*n);
    memset(gamma, 0, sizeof(int)*n);
    
    *m = n;
    memset(k, 0, sizeof(int));

    

    int alpha_counter = 0;
    int curr_x_ind = 0;

    int cumsum_beta = 0;
    int cumsum_beta_prev = 0;

    int cumsum_gamma = 0;
    int cumsum_gamma_next = 0;

    for (int i = 0; i < len_a_union_b-1; ++i)
    {
        if(alpha_counter>= n){
            break;
        }
        if(debug){
            printf("a_union_b[%d] = %f\n", i, a_union_b[i]);

        }
        for (int j = 0; j < n; ++j)
        {
            if ( (a_union_b[i] <= x[j]) && (x[j]<= a_union_b[i+1]))
            {
                alpha[alpha_counter]+=1;
                curr_x_ind = j;
            }
        }
        if(debug){
            printf("alpha = %d\n", alpha[alpha_counter]);

        }
        


        if(alpha[alpha_counter]>0){

            if(alpha_counter>=1){
                cumsum_beta = 0;
                for (int z = 0; z < n; ++z)
                {
                    if(b[z] < x[curr_x_ind]){
                        cumsum_beta++;
                    }else{break;}
                }
                beta[alpha_counter-1] = cumsum_beta - cumsum_beta_prev;
            }

            cumsum_gamma_next = 0;
            for (int z = n-1; z >=0; --z)
                {
                    if(a[z] > x[curr_x_ind]){
                        cumsum_gamma_next++;
                    }else{break;}
                }

            if(alpha_counter>=1){
                gamma[alpha_counter-1] = cumsum_gamma - cumsum_gamma_next;
            }

            cumsum_gamma = cumsum_gamma_next;


            
            cumsum_beta_prev = cumsum_beta;
            alpha_counter++;

        }
        if(debug){
            printf("alphacounter = %d\n", alpha_counter);

        }
    }

    *k = alpha_counter;

    return;
}
