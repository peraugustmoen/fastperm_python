#include "header.h"


static PyObject *C_get_log_permanents(PyObject *self, PyObject *args) {
//(SEXP XSEXP, SEXP aSEXP, SEXP bSEXP, SEXP nSEXP, SEXP TSEXP, SEXP debugSEXP){
	

	// X MUST ALREADY BE SORTED!!!!!!!!
	
	// REMEMBER THAT PYTHON IS COLUMN MAJOR! FIX THIS PLS


	PyArrayObject* Xo;
	PyArrayObject* ao;
	PyArrayObject* bo;
	int n;
	int T;
	int debug;

	if (!PyArg_ParseTuple(args, "O!O!O!iii", &PyArray_Type, &Xo,&PyArray_Type, &ao,&PyArray_Type, &bo, &n, &T,&debug)){
        return NULL;
  	}
  	if (PyArray_NDIM(Xo)!= 2 || PyArray_NDIM(ao) != 1 || PyArray_NDIM(bo) != 1 || !PyArray_ISFLOAT(Xo) || !PyArray_ISFLOAT(ao) 
  		|| !PyArray_ISFLOAT(bo)) {
        PyErr_SetString(PyExc_ValueError, "arrays must be one-dimensional and of type float");
    	/*printf("%d\n",PyArray_NDIM(Xo));
    	printf("%d\n",PyArray_NDIM(ao));
    	printf("%d\n",PyArray_NDIM(bo));
    	printf("%d\n",PyArray_ISFLOAT(Xo));
    	printf("%d\n",PyArray_ISFLOAT(ao));
    	printf("%d\n",PyArray_ISFLOAT(bo));*/

        return NULL;
    }


    double *X = PyArray_DATA(Xo);
    double *a = PyArray_DATA(ao);
    double *b = PyArray_DATA(bo);

	PyArray_Sort(ao,0,0);
	PyArray_Sort(bo,0,0);

	//R_qsort(a, 1, n);
	//R_qsort(b, 1, n);


	

	double * logperms = (double*)  malloc(sizeof(double) * T);
	memset(logperms, 0, sizeof(double)*T);


	double * a_union_b = (double*)  malloc(sizeof(double) * 2*n);
	memset(a_union_b, 0, sizeof(double)*2*n);

	int len_a_union_b =0;

	get_union(n, a, b, &len_a_union_b, a_union_b);

	/*SEXP alphaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP betaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP gammaSEXP = PROTECT(allocVector(INTSXP, n));
	SEXP logfactorialsSEXP = PROTECT(allocVector(REALSXP, n+1));
	SEXP mSEXP = PROTECT(allocVector(INTSXP, 1));
	SEXP k_SEXP = PROTECT(allocVector(INTSXP, 1));*/

	
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

	
	
/*	SEXP historySEXP = PROTECT(allocVector(INTSXP, 3*n));
	SEXP amount_historySEXP = PROTECT(allocVector(INTSXP, 6*n));*/

	int * history = (int * ) malloc(sizeof(int)*3*n);
	int * amount_history = (int * ) malloc(sizeof(int)*6*n);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);


	for (int t = 0; t < T; ++t)
	{
		//printf("t = %d\n",t);
		//double * x = X + (t*n);
		double * x = &(X[t*n]);
		
		
		/*if(t ==0){
			printf("x[0] = %f, x[1] = %f, x[n] = %f\n", x[0], x[1], x[n-1]);
		}*/

		//printf("x[0] = %f, x[1] = %f, x[n] = %f\n", x[0], x[1], x[n-1]);
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
	    	printf("T=%d, t=%d\n", T, t);
	    	printf("len_a_union_b = %d\n", len_a_union_b);
	    	printf("x:\n");
	    	print_float_vector(n,x);
	    	printf("a:\n");
	    	print_float_vector(n,a);
	    	printf("b:\n");
	    	print_float_vector(n,b);
	    	printf("a_union_b:\n");
	    	print_float_vector(2*n,a_union_b);
	    	printf("len a_union_b:%d\n", len_a_union_b);
	    	printf("alpha:\n");
	    	print_int_vector(n,  alpha);
	    	printf("beta:\n");
	    	print_int_vector(n,  beta);
	    	printf("gamma:\n");
	    	print_int_vector(n,  gamma);
	    	printf("m:%d\n", *m);
	    	printf("k:%d\n", *k);
	    	
	    }

		int history_len = 0;

	
		memset(history, 0, sizeof(int)*3*n);
		memset(amount_history, 0, sizeof(int)*6*n);


		if(debug){
			printf("REDUCING NOW\n");
		}
		
		//return(XSEXP);
		int result = reduction(alpha,  beta,  gamma, m, n, k, history,
				   amount_history, &history_len, debug);

		if(result != 0){

			printf("Error recorded, rerunning and returning NULL");

			memset(alpha, 0, sizeof(int)*n);
			memset(beta, 0, sizeof(int)*n);
			memset(gamma, 0, sizeof(int)*n);
			memset(m, 0, sizeof(int));
			memset(k, 0, sizeof(int));
			debug = 1;
			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);


		    if(debug){
		    	printf("len_a_union_b = %d\n", len_a_union_b);
		    	printf("x:\n");
		    	print_float_vector(n,x);
		    	printf("a:\n");
		    	print_float_vector(n,a);
		    	printf("b:\n");
		    	print_float_vector(n,b);
		    	printf("a_union_b:\n");
		    	print_float_vector(2*n,a_union_b);
		    	printf("len a_union_b:%d\n", len_a_union_b);
		    	printf("alpha:\n");
		    	print_int_vector(n,  alpha);
		    	printf("beta:\n");
		    	print_int_vector(n,  beta);
		    	printf("gamma:\n");
		    	print_int_vector(n,  gamma);
		    	printf("m:%d\n", *m);
		    	printf("k:%d\n", *k);
		    	
		    }

			int history_len = 0;

		
			memset(history, 0, sizeof(int)*3*n);
			memset(amount_history, 0, sizeof(int)*6*n);


			if(debug){
				printf("REDUCING NOW\n");
			}
			
			//return(XSEXP);
			int result = reduction(alpha,  beta,  gamma, m, n, k, history,
					   amount_history, &history_len, debug);


			free_dictionary(new_log_subperms);
			free_dictionary(old_log_subperms);

			
			return NULL;
		}

		if(debug){
			printf("history len = %d\n", history_len);

			printf("REDUCED SUBPERMS\n");
		}
		sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
						log_factorials, n,  m, k);
		//printf("RESULT:\n");
		//print_matrix(n+1, n+1, new_log_subperms);

		dictionary * tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;



		if(debug){
			printf("==========\nReverse reduction:\n==========\n");
		}
		//printf("old = %d\n", old_log_subperms);
		//printf("new = %d\n", new_log_subperms);
		dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
						   beta,  gamma, m,  n, k,  history,
				           amount_history, &history_len, log_factorials);

		


		
		double logperm =  Csparse_log_sum_exp(the_log_subperms);
		logperms[t] = logperm;
		if(debug){
			printf("logperm = %f\n", logperm);

		}



	}
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);

	// FREE MALLOCED STUFF!!!!!!!!!!!
	// 
	// 
	long dims[1];
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


static PyMethodDef fastpermMethods[] = {
  {"get_log_permanents", C_get_log_permanents, METH_VARARGS, "Function for computing log permanents"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fastperm = {
  PyModuleDef_HEAD_INIT,
  "fastperm",
  "Custom fastperm module",
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