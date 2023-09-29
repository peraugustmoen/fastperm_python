#define cord(r,c) ((r) + (DEPTH)*(c))
#define cord_spec(r,c, D) ((r) + (D)*(c))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <Python.h>
#define PY_SSIZE_T_CLEAN
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#include "xxhash.h"
#define RATIO 10

// structs 
typedef struct {
    int x;
    int y;
} pair;

typedef struct {
	double * value_array; // log perms
	pair * array; // array with index
	int * table;
	size_t used_len; // current size of dict / number of vals
	size_t table_size; // current max size
	//SEXP arraySEXP;
	//SEXP value_arraySEXP;
} dictionary;


// Function declaration 

size_t hash(pair, int);
double Clog_sum_exp(double *, int, double);
void print_matrix(int, int , double * );
void print_int_vector(int ,  int * );
void print_float_vector(int ,  double * );

static PyObject *C_get_log_permanents(PyObject *, PyObject *);
static PyObject *log_sum_exp(PyObject *, PyObject *);
static PyObject *C_get_log_ML(PyObject *, PyObject *);
static PyObject *C_get_log_permanents_bioassay(PyObject *, PyObject *);
static PyObject *C_get_log_ML_bioassay(PyObject *, PyObject *);

double Clog_sum_exp(double * , int, double );
double Csparse_log_sum_exp(dictionary *);


int nonzero_perm(double * , double * , double * , int );

void get_union(int, double *, double * , int * , double *);


void update_dict(pair, double, dictionary *);
void print_sparse_matrix(dictionary * );

void sparse_reverse_tt(dictionary * ,  dictionary *,
						int, int,
						int * , int * , int *, double *, int, int *,
						int * );

void sparse_reverse_bt(dictionary * ,  dictionary * ,
						int , int ,
						int * , int *, int *, double *, int , int * ,
						int * );

void sparse_reverse_lm(dictionary *,  dictionary * ,
						int , int,
						int * , int * , int * , double * , int , int * ,
						int * );

void sparse_reverse_bs(dictionary * ,  dictionary * ,
						int, int ,
						int * , int *, int * , double * , int , int * ,
						int * );
void sparse_reverse_rm(dictionary * ,  dictionary * ,
						int, int ,
						int * , int *, int * , double * , int , int * ,
						int * );
void sparse_reverse_ts(dictionary * ,  dictionary * ,
						int, int ,
						int * , int *, int * , double * , int , int * ,
						int * );
void sparse_get_reduced_log_subperms(dictionary *, int * , int * , int *, 
				double * , int, int *, int *);

dictionary* sparse_reverse_reduction(dictionary * , dictionary * , int * , 
					   int * , int * , int *, int , int *, int * ,
			           int * , int * , double * );

double Csparse_log_sum_exp(dictionary *);


double lookup_logperm(dictionary *, pair );

void add_to_dictionary(dictionary * , pair, double);
dictionary * init_dictionary(size_t);

void free_dictionary(dictionary *);

void expand_dictionary(dictionary * );
void nullset_dictionary(dictionary * );

// reverse top trim
void reverse_tt(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				int * );

// Reverse bottom trim
void reverse_bt(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				int * );

// Reverse bottom split
void reverse_bs(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );
// Reverse top split
void reverse_ts(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );


// Reverse left merge
void reverse_lm(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );
// Reverse right merge
void reverse_rm(double * , double * , int , int ,
				int * , int * , int * , double * , int , int * , 
				double * , int * );

void get_reduced_log_subperms(double * , int *, int *, int * , 
				double * , int, int *, int * );

int reduction(int * , int *, int * , int *, int, int *, int *,
			   int * , int * , int);

double * reverse_reduction(double * , double * , int * , 
					   int * , int * , int *, int , int *, int * ,
			           int *, int * , double * , 
			           double *);

void get_alphabetagamma(double * , int , double * , double * , double * , int , int *, 
    int *, int *, int * , int *, int );

char check_if_reduced(int *, int *, int *, 
					  int * , int * );
