#include "header.h"

static PyObject *C_get_log_permanents(PyObject *self, PyObject *args) {

	PyArrayObject* Xo; // X (python object)
	PyArrayObject* ao; // a (python object)
	PyArrayObject* bo; // b (python object)
	int n;
	int T;
	int debug;

	if (!PyArg_ParseTuple(args, "O!O!O!iii", &PyArray_Type, &Xo,&PyArray_Type, &ao,&PyArray_Type, &bo, &n, &T,&debug)){
        return NULL;
  	}

  	if(PyArray_NDIM(Xo)!= 2){
  		if(T!=1){
  			PyErr_SetString(PyExc_ValueError, "X must be 2-dimensional whenever T> 1");
  			return NULL;
  		}
  	}
  	if(!PyArray_ISFLOAT(Xo)){
  		PyErr_SetString(PyExc_ValueError, "X must be of type float");
  		return NULL;
  	}
  	if(PyArray_TYPE(Xo)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "X must have dtype numpy.float64");
  		return NULL;
  	}
  	if(PyArray_TYPE(ao)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "a must have dtype numpy.float64");
  		return NULL;
  	}
  	if(PyArray_TYPE(bo)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "b must have dtype numpy.float64");
  		return NULL;
  	}
  	if( PyArray_NDIM(ao) != 1 || !PyArray_ISFLOAT(ao) ){
  		PyErr_SetString(PyExc_ValueError, "a must be one-dimensional and of type float");
  		return NULL;
  	}
  	if( PyArray_NDIM(bo) != 1 || !PyArray_ISFLOAT(bo) ){
  		PyErr_SetString(PyExc_ValueError, "b must be one-dimensional and of type float");
  		return NULL;
  	}

  	npy_intp * shapeX = PyArray_SHAPE(Xo);

  	if( (int)shapeX[0] != T || (int)shapeX[1] != n ){
  		PyErr_SetString(PyExc_ValueError, "X must be T x n");
  		return NULL;
  	}

  	//fprintf(stdout,"fortran style = %d\n", PyArray_IS_F_CONTIGUOUS(Xo));
  	if(PyArray_IS_F_CONTIGUOUS(Xo)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming to C style.\n");
  		}
  		Xo = (PyArrayObject *)PyArray_CastToType(Xo, PyArray_DESCR(Xo), 0);
  	}

  	//fprintf(stdout,"fortran style = %d\n", PyArray_IS_F_CONTIGUOUS(Xo));

  	if(debug){
  		fprintf(stdout,"dim(X)[0] = %d, dim(X)[1] = %d\n", (int)shapeX[0], (int)shapeX[1]);

  	}
  	

    double *X = PyArray_DATA(Xo);
    double *a = PyArray_DATA(ao);
    double *b = PyArray_DATA(bo);

	PyArray_Sort(ao,0,NPY_QUICKSORT);
	PyArray_Sort(bo,0,NPY_QUICKSORT);
	PyArray_Sort(Xo, 1, NPY_QUICKSORT);



	

	double * logperms = (double*)  malloc(sizeof(double) * T);
	memset(logperms, 0, sizeof(double)*T);


	double * a_union_b = (double*)  malloc(sizeof(double) * 2*n);
	memset(a_union_b, 0, sizeof(double)*2*n);

	int len_a_union_b =0;

	get_union(n, a, b, &len_a_union_b, a_union_b);


	
	int * alpha = (int*) malloc(sizeof(int) * n);
	int * beta = (int*) malloc(sizeof(int) * n);
	int * gamma = (int*) malloc(sizeof(int) * n);

	
	double * log_factorials =(double*) malloc(sizeof(double) * (n+1));
	int * m = (int*) malloc(sizeof(int) );
	int * k = (int*) malloc(sizeof(int) );


	dictionary * new_log_subperms = init_dictionary(n);
	dictionary * old_log_subperms = init_dictionary(n);

	
	memset(alpha, 0, sizeof(int)*n);
	memset(beta, 0, sizeof(int)*n);
	memset(gamma, 0, sizeof(int)*n);
	memset(log_factorials, 0, sizeof(double)*(n+1));
	memset(m, 0, sizeof(int));
	memset(k, 0, sizeof(int));

	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	

	int * history = (int * ) malloc(sizeof(int)*3*n);
	int * amount_history = (int * ) malloc(sizeof(int)*6*n);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);


	for (int t = 0; t < T; ++t)
	{
		//fprintf(stdout,"t = %d\n",t);
		//double * x = X + (t*n);
		double * x = &(X[t*n]);
		
		
		/*if(t ==0){
			fprintf(stdout,"x[0] = %f, x[1] = %f, x[n] = %f\n", x[0], x[1], x[n-1]);
		}*/

		//fprintf(stdout,"x[0] = %f, x[1] = %f, x[n] = %f\n", x[0], x[1], x[n-1]);
		if(!nonzero_perm(x, a,  b, n)){
			logperms[t] = -1;
			continue;
		}
		memset(alpha, 0, sizeof(int)*n);
		memset(beta, 0, sizeof(int)*n);
		memset(gamma, 0, sizeof(int)*n);
		memset(m, 0, sizeof(int));
		memset(k, 0, sizeof(int));

		get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
	    beta, gamma,  k, m, debug);


	    if(debug){
	    	fprintf(stdout,"T=%d, t=%d\n", T, t);
	    	fprintf(stdout,"len_a_union_b = %d\n", len_a_union_b);
	    	fprintf(stdout,"x:\n");
	    	print_float_vector(n,x);
	    	fprintf(stdout,"a:\n");
	    	print_float_vector(n,a);
	    	fprintf(stdout,"b:\n");
	    	print_float_vector(n,b);
	    	fprintf(stdout,"a_union_b:\n");
	    	print_float_vector(2*n,a_union_b);
	    	fprintf(stdout,"len a_union_b:%d\n", len_a_union_b);
	    	fprintf(stdout,"alpha:\n");
	    	print_int_vector(n,  alpha);
	    	fprintf(stdout,"beta:\n");
	    	print_int_vector(n,  beta);
	    	fprintf(stdout,"gamma:\n");
	    	print_int_vector(n,  gamma);
	    	fprintf(stdout,"m:%d\n", *m);
	    	fprintf(stdout,"k:%d\n", *k);
	    	
	    }

		int history_len = 0;

	
		memset(history, 0, sizeof(int)*3*n);
		memset(amount_history, 0, sizeof(int)*6*n);


		if(debug){
			fprintf(stdout,"REDUCING NOW\n");
		}
		
		//return(XSEXP);
		int result = reduction(alpha,  beta,  gamma, m, n, k, history,
				   amount_history, &history_len, debug);

		if(result != 0){

			fprintf(stdout,"Error recorded, rerunning and returning NULL");

			memset(alpha, 0, sizeof(int)*n);
			memset(beta, 0, sizeof(int)*n);
			memset(gamma, 0, sizeof(int)*n);
			memset(m, 0, sizeof(int));
			memset(k, 0, sizeof(int));
			debug = 1;
			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);


		    if(debug){
		    	fprintf(stdout,"len_a_union_b = %d\n", len_a_union_b);
		    	fprintf(stdout,"x:\n");
		    	print_float_vector(n,x);
		    	fprintf(stdout,"a:\n");
		    	print_float_vector(n,a);
		    	fprintf(stdout,"b:\n");
		    	print_float_vector(n,b);
		    	fprintf(stdout,"a_union_b:\n");
		    	print_float_vector(2*n,a_union_b);
		    	fprintf(stdout,"len a_union_b:%d\n", len_a_union_b);
		    	fprintf(stdout,"alpha:\n");
		    	print_int_vector(n,  alpha);
		    	fprintf(stdout,"beta:\n");
		    	print_int_vector(n,  beta);
		    	fprintf(stdout,"gamma:\n");
		    	print_int_vector(n,  gamma);
		    	fprintf(stdout,"m:%d\n", *m);
		    	fprintf(stdout,"k:%d\n", *k);
		    	
		    }

			int history_len = 0;

		
			memset(history, 0, sizeof(int)*3*n);
			memset(amount_history, 0, sizeof(int)*6*n);


			if(debug){
				fprintf(stdout,"REDUCING NOW\n");
			}
			
			//return(XSEXP);
			int result = reduction(alpha,  beta,  gamma, m, n, k, history,
					   amount_history, &history_len, debug);


			free_dictionary(new_log_subperms);
			free_dictionary(old_log_subperms);

			PyErr_SetString(PyExc_RuntimeError, sprintf("Failed to compute permanent of sample t=%d\n", t));

			return NULL;
		}

		if(debug){
			fprintf(stdout,"history len = %d\n", history_len);

			fprintf(stdout,"REDUCED SUBPERMS\n");
		}
		sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
						log_factorials, n,  m, k);
		//fprintf(stdout,"RESULT:\n");
		//print_matrix(n+1, n+1, new_log_subperms);

		dictionary * tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;



		if(debug){
			fprintf(stdout,"==========\nReverse reduction:\n==========\n");
		}
		//fprintf(stdout,"old = %d\n", old_log_subperms);
		//fprintf(stdout,"new = %d\n", new_log_subperms);
		dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
						   beta,  gamma, m,  n, k,  history,
				           amount_history, &history_len, log_factorials);

		


		
		double logperm =  Csparse_log_sum_exp(the_log_subperms);
		logperms[t] = logperm;
		if(debug){
			fprintf(stdout,"logperm = %f\n", logperm);

		}



	}
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);


	npy_intp dims[1];
	dims[0] = T;
	
	free(a_union_b);
	free(alpha);
	free(beta);
	free(gamma);
	free(log_factorials);
	free(m);
	free(k);
	free(history);
	free(amount_history);

	return(PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, logperms));

}



static PyObject *log_sum_exp(PyObject *self, PyObject *args) {

	PyArrayObject* arrayo;

	if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &arrayo)){
        return NULL;
  	}

  	
  	npy_intp ndim = PyArray_NDIM(arrayo);
  	npy_intp * shape = PyArray_SHAPE(arrayo);

  	npy_intp totsize=1;

  	for (int i = 0; i < ndim; ++i)
  	{
  		totsize = totsize * shape[i];
  	}

  	

    double *array = PyArray_DATA(arrayo);
    
    // find max
    
    double maxval = array[0];

    for (int i = 1; i < totsize; ++i)
    {
    	if(array[i]>maxval){
    		maxval = array[i];
    	}
    }

    double exp_result = 0;



	for (int i = 0; i < totsize; ++i)
	{
		/*if(array[i]<0){
			continue;
		}*/

		exp_result += exp(array[i] - maxval);
	}

	////printf("res = %f\n", (maxval + log(exp_result)));
	return PyFloat_FromDouble(maxval + log(exp_result));

}


static PyMethodDef fastpermMethods[] = {
  {"get_log_permanents", C_get_log_permanents, METH_VARARGS, "get_log_permanents(X, a, b,n,T,debug)\n\
\n\
Computes log permanents (i.e.importance weights)\n\
associated with a simulated data set X and \n\
data a,b.\n\
\n\
Given a matrix X of samples from a proposal distribution,\n\
and vectors a,b containing left and right censoring points,\n\
the function returns a vector of log permanents corresponding \n\
to each sample. \n\
\n\
Parameters \n\
---------- \n\
X : ndarray \n\
    A numpy array of dimension T x n, in \n\
    which each row contains a sample from \n\
    the proposal distribution. \n\
a : ndarray\n\
    A flat numpy array of length n\n\
    containing the left censoring\n\
    points of the data.\n\
b : ndarray\n\
    A flat numpy array of length n\n\
    containing the right censoring\n\
    points of the data.\n\
n : int\n\
    Sample size.\n\
T : int\n\
    Number of samples from the\n\
    proposal distribution.\n\
debug : boolean\n\
    If true, debug information\n\
    is printed to stdout.\n\
\n\
Returns \n\
------- \n\
ndarray \n\
    Numpy array of log permanents,\n\
    each element associated to \n\
    the corresponding row in X.\n"},
  {"log_sum_exp", log_sum_exp, METH_VARARGS, "log_sum_exp(array)\n\
\n\
Computes the log sum exp of an array. \n\
\n\
Given input array = [x_1, ..., x_n], function returns \n\
x_* + log(exp(x_1 - x_*) + ... + exp(x_n - x_*)), \n\
where x_* = max(x_1, ... x_n).\n\
\n\
Parameters \n\
---------- \n\
array : ndarray \n\
    data array \n\
\n\
Returns \n\
------- \n\
float \n"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fastperm = {
  PyModuleDef_HEAD_INIT,
  "fastperm",
  "Module for computing permanents of a block rectangular matrix",
  -1,
  fastpermMethods
};

PyMODINIT_FUNC PyInit_fastperm(void)
{
    import_array();
    return PyModule_Create(&fastperm);
}

int nonzero_perm(double * x, double * a, double * b, int n){

	for (int i = 0; i < n; ++i)
	{
		if(x[i]< a[i] || x[i] > b[i]){
			return 0;
		}
	}
	return 1;

}