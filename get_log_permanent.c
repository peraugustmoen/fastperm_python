#include "header.h"

static PyObject *C_get_log_perms(PyObject *self, PyObject *args) {

	PyArrayObject* Xo; // X (python object)
	PyArrayObject* to; // t (python object)
	PyArrayObject* yo; // y (python object)
	int n;
	int S;
	int debug;
	int constant_ts=1;

	if (!PyArg_ParseTuple(args, "O!O!O!i", &PyArray_Type, &Xo,&PyArray_Type, &to,&PyArray_Type, &yo, &debug)){
        return NULL;
  	}

  	S = 1;
  	n = 1;

  	if(PyArray_NDIM(Xo)== 2){
  		S = 2;
  	}else if(PyArray_NDIM(Xo)!= 1){
  		PyErr_SetString(PyExc_ValueError, "X must be either one-dimensional or two-dimensional");
  		return NULL;
  	}
  	if(!PyArray_ISFLOAT(Xo)){
  		PyErr_SetString(PyExc_ValueError, "X must be of type float");
  		return NULL;
  	}
  	if(PyArray_TYPE(Xo)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "X must have dtype numpy.float64");
  		return NULL;
  	}
  	if(PyArray_TYPE(to)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "t must have dtype numpy.float64");
  		return NULL;
  	}
  	if(PyArray_TYPE(yo)!= NPY_INT32){
  		PyErr_SetString(PyExc_ValueError, "y must have dtype numpy.int32");
  		return NULL;
  	}
  	if( PyArray_NDIM(to) != 1 ){
  		constant_ts = 0;
  		if(PyArray_NDIM(to) != 2){
  			PyErr_SetString(PyExc_ValueError, "t must be either one-dimensional or two-dimensional");
  			return NULL;
  		}
  	}
  	if( PyArray_NDIM(yo) != 1 ){
  		PyErr_SetString(PyExc_ValueError, "y must be one-dimensional");
  		return NULL;
  	}

  	npy_intp * shapeX = PyArray_SHAPE(Xo);
  	npy_intp * shapet = PyArray_SHAPE(to);
  	npy_intp * shapey = PyArray_SHAPE(yo);

  	if(S>1){
  		n = (int)shapeX[1];	
  		S = (int)shapeX[0];
	}else{
		n = (int)shapeX[0];	
	} 

  	if ((int)shapey[0] != n)
  		{
  			PyErr_SetString(PyExc_ValueError, "y must be of length n");
  			return NULL;
  		}
  	if(constant_ts){
  		if ((int)shapet[0] != n)
  		{
  			PyErr_SetString(PyExc_ValueError, "t must be of length n when one-dimensional");
  			return NULL;
  		}
  		
  	}else{
  		if( (int)shapet[0] != S || (int)shapet[1] != n ){
  			PyErr_SetString(PyExc_ValueError, "t must be S x n when two-dimensional");
  			return NULL;
  		}
  	}
  	

  	if(PyArray_IS_F_CONTIGUOUS(Xo)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming X to C style.\n");
  		}
  		Xo = (PyArrayObject *)PyArray_CastToType(Xo, PyArray_DESCR(Xo), 0);
  	}

  	if(PyArray_IS_F_CONTIGUOUS(yo)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming y to C style.\n");
  		}
  		yo = (PyArrayObject *)PyArray_CastToType(yo, PyArray_DESCR(yo), 0);
  	}

  	if(PyArray_IS_F_CONTIGUOUS(to)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming t to C style.\n");
  		}
  		to = (PyArrayObject *)PyArray_CastToType(to, PyArray_DESCR(to), 0);
  	}


  	if(debug){
  		if(S>1){
  			fprintf(stdout,"dim(X)[0] = %d, dim(X)[1] = %d\n", (int)shapeX[0], (int)shapeX[1]);
  		}else{
  			fprintf(stdout,"dim(X)[0] = %d\n", (int)shapeX[0]);
  		}
  		fprintf(stdout,"dim(y)[0] = %d\n", (int)shapey[0]);
  		if(!constant_ts){
  			fprintf(stdout,"dim(t)[0] = %d, dim(t)[1] = %d\n", (int)shapet[0], (int)shapet[1]);
  		}else{
  			fprintf(stdout,"dim(t)[0] = %d\n", (int)shapet[0]);
  			
  		}

  	}
  	

    double *X = PyArray_DATA(Xo);
    double *t = PyArray_DATA(to);
    int *y = PyArray_DATA(yo);

    double * a = (double*) PyMem_Calloc(n, sizeof(double));
    double * b = (double*) PyMem_Calloc(n, sizeof(double));

    PyArrayObject* ao;
    PyArrayObject* bo;

    npy_intp dimss[1];
	dimss[0] = n;
    if (constant_ts)
    {

	    for (int i = 0; i < n; ++i)
	    {
	    	if(y[i]==1){
	    		b[i] = t[i];
	    		a[i] = -DBL_MAX;
	    	}else{
	    		b[i] = DBL_MAX;
	    		a[i] = t[i];

	    	}
	    }

    	ao = PyArray_SimpleNewFromData(1, dimss, NPY_FLOAT64, a);
    	bo = PyArray_SimpleNewFromData(1, dimss, NPY_FLOAT64, b);
    	PyArray_Sort(ao,0,NPY_QUICKSORT);
		PyArray_Sort(bo,0,NPY_QUICKSORT);
		/*printf("a:\n");
		for (int i = 0; i < n; ++i)
		{
			printf("%f\n", a[i]);
		}
		printf("b:\n");
		for (int i = 0; i < n; ++i)
		{
			printf("%f\n", b[i]);
		}*/
    }
    

	if(S>1){
		PyArray_Sort(Xo, 1, NPY_QUICKSORT);
	}else{
		PyArray_Sort(Xo, 0, NPY_QUICKSORT);
	}



	

	double * log_perms = (double*)  PyMem_Malloc(sizeof(double) * S);
	memset(log_perms, 0, sizeof(double)*S);


	double * a_union_b = (double*)  PyMem_Malloc(sizeof(double) * 2*n);
	int len_a_union_b=0;
	if (constant_ts)
	{
		memset(a_union_b, 0, sizeof(double)*2*n);

		len_a_union_b =0;

		get_union(n, a, b, &len_a_union_b, a_union_b);
	}
	

	
	int * alpha = (int*) PyMem_Malloc(sizeof(int) * n);
	int * beta = (int*) PyMem_Malloc(sizeof(int) * n);
	int * gamma = (int*) PyMem_Malloc(sizeof(int) * n);

	
	double * log_factorials =(double*) PyMem_Malloc(sizeof(double) * (n+1));
	int * m = (int*) PyMem_Malloc(sizeof(int) );
	int * k = (int*) PyMem_Malloc(sizeof(int) );


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

	

	int * history = (int * ) PyMem_Malloc(sizeof(int)*3*n);
	int * amount_history = (int * ) PyMem_Malloc(sizeof(int)*6*n);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);


	for (int s = 0; s < S; ++s)
	{
		double * x = &(X[s*n]);
		memset(alpha, 0, sizeof(int)*n);
		memset(beta, 0, sizeof(int)*n);
		memset(gamma, 0, sizeof(int)*n);
		memset(m, 0, sizeof(int));
		memset(k, 0, sizeof(int));

		if(constant_ts){
			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);
		    if(!nonzero_perm(x, a,  b, n)){
		    	// setting log_perms[s] to python NAN:
		    	log_perms[s] = NPY_NAN;
				continue;
			}
			if((*k)==1){
				log_perms[s] = log_factorials[n];
				continue;
			}
		}else{
		    for (int i = 0; i < n; ++i)
		    {
		    	if(y[i]==1){
		    		b[i] = t[i + s*n];
		    		a[i] = -DBL_MAX;
		    	}else{
		    		b[i] = DBL_MAX;
		    		a[i] = t[i + s*n];

		    	}
		    }

	    	ao = PyArray_SimpleNewFromData(1, dimss, NPY_FLOAT64, a);
	    	bo = PyArray_SimpleNewFromData(1, dimss, NPY_FLOAT64, b);
	    	PyArray_Sort(ao,0,NPY_QUICKSORT);
			PyArray_Sort(bo,0,NPY_QUICKSORT);

			if(!nonzero_perm(x, a,  b, n)){
				log_perms[s] = NPY_NAN;
				continue;
			}
			memset(a_union_b, 0, sizeof(double)*2*n);

			len_a_union_b =0;

			get_union(n, a, b, &len_a_union_b, a_union_b);

			get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
		    beta, gamma,  k, m, debug);
		    
			if((*k)==1){
				log_perms[s] = log_factorials[n];
				continue;
			}
	   
		}

		


	    if(debug){
	    	fprintf(stdout,"S=%d, s=%d\n", S, s);
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
			
			result = reduction(alpha,  beta,  gamma, m, n, k, history,
					   amount_history, &history_len, debug);


			free_dictionary(new_log_subperms);
			free_dictionary(old_log_subperms);

			PyErr_Format(PyExc_RuntimeError,
                 "Failed to compute permanent at iteration s=%d\n",s
                 );

			return NULL;
		}

		if(debug){
			fprintf(stdout,"history len = %d\n", history_len);

			fprintf(stdout,"REDUCED SUBPERMS\n");
		}
		sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
						log_factorials, n,  m, k);

		dictionary * tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;



		if(debug){
			fprintf(stdout,"==========\nReverse reduction:\n==========\n");
		}
		dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
						   beta,  gamma, m,  n, k,  history,
				           amount_history, &history_len, log_factorials);

		


		
		double logperm =  Csparse_log_sum_exp(the_log_subperms);
		log_perms[s] = logperm;
		if(debug){
			fprintf(stdout,"logperm = %f\n", logperm);

		}



	}
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);

	npy_intp dims[1];
	dims[0] = S;
	PyMem_Free(a_union_b);
	PyMem_Free(alpha);
	PyMem_Free(beta);
	PyMem_Free(gamma);
	PyMem_Free(log_factorials);
	PyMem_Free(m);
	PyMem_Free(k);
	PyMem_Free(history);
	PyMem_Free(amount_history);
	PyMem_Free(a);
	PyMem_Free(b);
	return(PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, log_perms));
	//return Py_BuildValue("i", *log_perms);

}
static PyObject *C_get_log_perms_bioassay(PyObject *self, PyObject *args) {

	PyArrayObject* Xo; // X (python object)
	PyArrayObject* levelso; // X (python object)
	PyArrayObject* successeso; // t (python object)
	PyArrayObject* trialso; // y (python object)
	int n;
	int num_trials;
	int S;
	int debug;


	if (!PyArg_ParseTuple(args, "O!O!O!O!i", &PyArray_Type, &Xo,&PyArray_Type, &levelso,&PyArray_Type, &successeso,&PyArray_Type, 
		&trialso,&debug)){
        return NULL;
  	}

  	if(PyArray_NDIM(Xo)== 2){
  		S = 2;
  	}else if(PyArray_NDIM(Xo)!= 1){
  		PyErr_SetString(PyExc_ValueError, "X must be either one or two-dimensional");
  			return NULL;
  	}
  	if(!PyArray_ISFLOAT(Xo)){
  		PyErr_SetString(PyExc_ValueError, "X must be of type float");
  		return NULL;
  	}
  	if(PyArray_TYPE(Xo)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "X must have dtype numpy.float64");
  		return NULL;
  	}
  	if(PyArray_TYPE(levelso)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "levels must have dtype numpy.float64");
  		return NULL;
  	}
  	if(PyArray_TYPE(successeso)!= NPY_INT32){
  		PyErr_SetString(PyExc_ValueError, "successes must have dtype numpy.int32");
  		return NULL;
  	}
  	if(PyArray_TYPE(trialso)!= NPY_INT32){
  		PyErr_SetString(PyExc_ValueError, "trials must have dtype numpy.int32");
  		return NULL;
  	}


  	if( PyArray_NDIM(levelso) != 1 ){
  		PyErr_SetString(PyExc_ValueError, "levels must be one-dimensional");
  		return NULL;

  	}
  	if( PyArray_NDIM(trialso) != 1 ){
  		PyErr_SetString(PyExc_ValueError, "trials must be one-dimensional");
  		return NULL;
  		
  	}
  	if( PyArray_NDIM(successeso) != 1 ){
  		PyErr_SetString(PyExc_ValueError, "successes must be one-dimensional");
  		return NULL;
  		
  	}


  	npy_intp * shapeX = PyArray_SHAPE(Xo);
  	npy_intp * shapetrials = PyArray_SHAPE(trialso);
  	npy_intp * shapesuccesses = PyArray_SHAPE(successeso);
  	npy_intp * shapelevels = PyArray_SHAPE(levelso);

  	if(S>1){
	  	S = (int)shapeX[0];
	  	n = (int)shapeX[1];
	}else{
		S = 1;
		n = (int)shapeX[0];
	}
	num_trials = (int)shapelevels[0];

  	if((int)shapetrials[0] != num_trials){
  		PyErr_SetString(PyExc_ValueError, "trials must have same length as levels");
  		return NULL;
  	}
  	if((int)shapesuccesses[0] != num_trials){
  		PyErr_SetString(PyExc_ValueError, "successes must have same length as levels");
  		return NULL;
  	}
  	


  	if(PyArray_IS_F_CONTIGUOUS(Xo)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming X to C style.\n");
  		}
  		Xo = (PyArrayObject *)PyArray_CastToType(Xo, PyArray_DESCR(Xo), 0);
  	}

  	if(PyArray_IS_F_CONTIGUOUS(trialso)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming trials to C style.\n");
  		}
  		trialso = (PyArrayObject *)PyArray_CastToType(trialso, PyArray_DESCR(trialso), 0);
  	}

  	if(PyArray_IS_F_CONTIGUOUS(successeso)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming successes to C style.\n");
  		}
  		successeso = (PyArrayObject *)PyArray_CastToType(successeso, PyArray_DESCR(successeso), 0);
  	}
  	if(PyArray_IS_F_CONTIGUOUS(levelso)){
  		if(debug){
  			fprintf(stdout,"Fortran style memory layout detected. Transforming levels to C style.\n");
  		}
  		levelso = (PyArrayObject *)PyArray_CastToType(levelso, PyArray_DESCR(levelso), 0);
  	}


  	if(debug){
  		if(S>1){
  			fprintf(stdout,"dim(X)[0] = %d, dim(X)[1] = %d\n", (int)shapeX[0], (int)shapeX[1]);
  		}else{
  			fprintf(stdout,"dim(X)[0] = %d\n", (int)shapeX[0]);
  		}

		fprintf(stdout,"dim(successes)[0] = %d\n", (int)shapesuccesses[0]);
		fprintf(stdout,"dim(trials)[0] = %d\n", (int)shapetrials[0]);
		fprintf(stdout,"dim(levels)[0] = %d\n", (int)shapelevels[0]);
		fprintf(stdout,"n = %d, num_trials = %d\n",n, num_trials);

  		

  	}
  	

    double *X = PyArray_DATA(Xo);
    double *levels = PyArray_DATA(levelso);
    int *trials = PyArray_DATA(trialso);
    int *successes = PyArray_DATA(successeso);

    int count = 0;
    for (int j = 0; j < num_trials; ++j)
    {
    	count += trials[j];
    }
    if(count != n){
    	PyErr_SetString(PyExc_ValueError, "The total number of trials must sum to n.");
    	return NULL;
    }

    double * a = (double*) PyMem_Calloc(n, sizeof(double));
    double * b = (double*) PyMem_Calloc(n, sizeof(double));

    PyArrayObject* ao;
    PyArrayObject* bo;

    npy_intp dimss[1];
	dimss[0] = n;
	ao = PyArray_SimpleNewFromData(1, dimss, NPY_FLOAT64, a);
	bo = PyArray_SimpleNewFromData(1, dimss, NPY_FLOAT64, b);

	int succ=0;
	int trial=0;
	int totcount = 0;
    
    for (int j = 0; j < num_trials; ++j)
    {
    	succ = successes[j];
    	trial = trials[j];
    	//printf("succ = %d, trial = %d\n", succ, trial);

    	for (int i = 0; i < succ; ++i)
    	{
    		b[totcount] = levels[j];
    		a[totcount++] = -DBL_MAX;
    		//printf("%d\n", totcount);
    	}

    	for (int i = succ; i < trial; ++i)
    	{
    		a[totcount] = levels[j];
    		b[totcount++] = DBL_MAX;
    		//printf("%d\n", totcount);
    	}


    }

	PyArray_Sort(ao,0,NPY_QUICKSORT);
	PyArray_Sort(bo,0,NPY_QUICKSORT);
	/*printf("a:\n");
	for (int i = 0; i < n; ++i)
	{
		printf("%f\n", a[i]);
	}
	printf("b:\n");
	for (int i = 0; i < n; ++i)
	{
		printf("%f\n", b[i]);
	}*/

    

	
	if(S>1){
		PyArray_Sort(Xo, 1, NPY_QUICKSORT);
	}else{
		PyArray_Sort(Xo, 0, NPY_QUICKSORT);
	}



	

	double * log_perms = (double*)  PyMem_Malloc(sizeof(double) * S);
	memset(log_perms, 0, sizeof(double)*S);


	double * a_union_b = (double*)  PyMem_Malloc(sizeof(double) * 2*n);
	int len_a_union_b=0;

	memset(a_union_b, 0, sizeof(double)*2*n);

	len_a_union_b =0;

	get_union(n, a, b, &len_a_union_b, a_union_b);

	

	
	int * alpha = (int*) PyMem_Malloc(sizeof(int) * n);
	int * beta = (int*) PyMem_Malloc(sizeof(int) * n);
	int * gamma = (int*) PyMem_Malloc(sizeof(int) * n);

	
	double * log_factorials =(double*) PyMem_Malloc(sizeof(double) * (n+1));
	int * m = (int*) PyMem_Malloc(sizeof(int) );
	int * k = (int*) PyMem_Malloc(sizeof(int) );


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

	

	int * history = (int * ) PyMem_Malloc(sizeof(int)*3*n);
	int * amount_history = (int * ) PyMem_Malloc(sizeof(int)*6*n);

	memset(history, 0, sizeof(int)*3*n);
	memset(amount_history, 0, sizeof(int)*6*n);


	for (int s = 0; s < S; ++s)
	{
		double * x = &(X[s*n]);
		memset(alpha, 0, sizeof(int)*n);
		memset(beta, 0, sizeof(int)*n);
		memset(gamma, 0, sizeof(int)*n);
		memset(m, 0, sizeof(int));
		memset(k, 0, sizeof(int));

		
		get_alphabetagamma(x, n, a, b, a_union_b, len_a_union_b, alpha, 
	    beta, gamma,  k, m, debug);
	    if(!nonzero_perm(x, a,  b, n)){
			log_perms[s] = NPY_NAN;
			continue;
		}
		if((*k)==1){
				log_perms[s] = log_factorials[n];
				continue;
			}


		


	    if(debug){
	    	fprintf(stdout,"S=%d, s=%d\n", S, s);
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
			
			result = reduction(alpha,  beta,  gamma, m, n, k, history,
					   amount_history, &history_len, debug);


			free_dictionary(new_log_subperms);
			free_dictionary(old_log_subperms);

			PyErr_Format(PyExc_RuntimeError,
                 "Failed to compute permanent of sample s=%d\n", s
                 );

			return NULL;
		}

		if(debug){
			fprintf(stdout,"history len = %d\n", history_len);

			fprintf(stdout,"REDUCED SUBPERMS\n");
		}
		sparse_get_reduced_log_subperms( new_log_subperms,  alpha, beta, gamma,
						log_factorials, n,  m, k);

		dictionary * tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;



		if(debug){
			fprintf(stdout,"==========\nReverse reduction:\n==========\n");
		}
		dictionary * the_log_subperms = sparse_reverse_reduction(old_log_subperms, new_log_subperms, alpha,
						   beta,  gamma, m,  n, k,  history,
				           amount_history, &history_len, log_factorials);

		


		
		double logperm =  Csparse_log_sum_exp(the_log_subperms);
		log_perms[s] = logperm;
		if(debug){
			fprintf(stdout,"logperm = %f\n", logperm);

		}



	}
	free_dictionary(new_log_subperms);
	free_dictionary(old_log_subperms);


	npy_intp dims[1];
	dims[0] = S;
	
	PyMem_Free(a_union_b);
	PyMem_Free(alpha);
	PyMem_Free(beta);
	PyMem_Free(gamma);
	PyMem_Free(log_factorials);
	PyMem_Free(m);
	PyMem_Free(k);
	PyMem_Free(history);
	PyMem_Free(amount_history);
	PyMem_Free(a);
	PyMem_Free(b);
	return(PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, log_perms));

}

static PyObject *C_get_log_ML_bioassay(PyObject *self, PyObject *args) {

	PyArrayObject* successeso; // X (python object)
	PyArrayObject* trialso; // t (python object)
	PyArrayObject* log_permso; // y (python object)
	int n;
	int num_trials;
	int S;
	int debug;

	if (!PyArg_ParseTuple(args, "O!O!O!ii",&PyArray_Type, &log_permso, &PyArray_Type, &successeso,
		&PyArray_Type, &trialso, &n,&debug)){
        return NULL;
  	}
  	if(PyArray_NDIM(trialso)!= 1){
		PyErr_SetString(PyExc_ValueError, "trials must be 1-dimensional");
		return NULL;
  	}

  	if(PyArray_TYPE(trialso)!= NPY_INT32){
  		PyErr_SetString(PyExc_ValueError, "trials must have dtype numpy.int32");
  		return NULL;
  	}

  	if(PyArray_NDIM(log_permso)!= 1){
		PyErr_SetString(PyExc_ValueError, "log_perms must be 1-dimensional");
		return NULL;
  	}
  	if(!PyArray_ISFLOAT(log_permso)){
  		PyErr_SetString(PyExc_ValueError, "log_perms must be of type float");
  		return NULL;
  	}
  	if(PyArray_TYPE(log_permso)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "log_perms must have dtype numpy.float64");
  		return NULL;
  	}
  	if(PyArray_NDIM(successeso)!= 1){
		PyErr_SetString(PyExc_ValueError, "successes must be 1-dimensional");
		return NULL;
  	}
  	if(PyArray_TYPE(successeso)!= NPY_INT32){
  		PyErr_SetString(PyExc_ValueError, "successes must be np.int32");
  		return NULL;
  	}


  	

  	npy_intp * shapetrials = PyArray_SHAPE(trialso);
  	npy_intp * shapesuccesses = PyArray_SHAPE(successeso);
  	npy_intp * shapelog_perms = PyArray_SHAPE(log_permso);

  	num_trials = (int)(shapetrials[0]);


  	if( (int)shapesuccesses[0] != num_trials){
  		PyErr_SetString(PyExc_ValueError, "successes must have same length as num_trials");
  		return NULL;
  	}

  	S = (int)shapelog_perms[0];

  	

  	int *trials = PyArray_DATA(trialso);
    double *log_perms = PyArray_DATA(log_permso);
    int *successes = PyArray_DATA(successeso);


  	double maxval = -1;

  	for (int i = 0; i < S; ++i)
  	{
  		if(!npy_isnan(log_perms[i])){
  			if(log_perms[i]> maxval){
  				maxval = log_perms[i];
  			}
  		}
  		
  	}
  	double result = 0;
 	

  	if(maxval<=-1){
  		result = NPY_NAN;
  		//PyErr_Format(PyExc_RuntimeError,
        //         "Error! No non-zero perms in log_perms\n"
        //         );
  		return PyFloat_FromDouble(result);
  	}
  	result = Clog_sum_exp(log_perms, S, maxval) - log((double)S);

  	// compute log factorials
  	double * log_factorials =(double*) PyMem_Malloc(sizeof(double) * (n+1));
	memset(log_factorials, 0, sizeof(double)*(n+1));

	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	result = result -log_factorials[n];
  	
  	for (int j = 0; j < num_trials; ++j)
  	{
  		result = result + log_factorials[trials[j]] - log_factorials[successes[j]] - log_factorials[trials[j] - successes[j]];
  	}

  	PyMem_Free(log_factorials);

  	return PyFloat_FromDouble(result);

}

static PyObject *C_get_log_ML(PyObject *self, PyObject *args) {

	PyArrayObject* log_permso; // y (python object)
	int n;
	int S;
	int debug;

	if (!PyArg_ParseTuple(args, "O!ii", &PyArray_Type, &log_permso, &n,&debug)){
        return NULL;
  	}


  	if(PyArray_NDIM(log_permso)!= 1){
		PyErr_SetString(PyExc_ValueError, "log_perms must be 1-dimensional");
		return NULL;
  	}
  	if(!PyArray_ISFLOAT(log_permso)){
  		PyErr_SetString(PyExc_ValueError, "log_perms must be of type float");
  		return NULL;
  	}
  	if(PyArray_TYPE(log_permso)!= NPY_DOUBLE){
  		PyErr_SetString(PyExc_ValueError, "log_perms must have dtype numpy.float64");
  		return NULL;
  	}
  	
  	if( PyArray_NDIM(log_permso) != 1 ){
  		
		PyErr_SetString(PyExc_ValueError, "log_perms must be one-dimensional");
		return NULL;
  		
  	}
  	

  	npy_intp * shapelog_perms = PyArray_SHAPE(log_permso);

  	S = (int)shapelog_perms[0];

    double *log_perms = PyArray_DATA(log_permso);


  	double maxval = -1;

  	for (int i = 0; i < S; ++i)
  	{
  		if(!npy_isnan(log_perms[i])){
  			if(log_perms[i]> maxval){
  				maxval = log_perms[i];
  			}
  		}
  		
  	}
  	double result = 0;
 	

  	if(maxval<=-1){
  		result = NPY_NAN;
  		//PyErr_Format(PyExc_RuntimeError,
        //         "Error! No non-zero perms in log_perms\n"
        //         );
  		return PyFloat_FromDouble(result);
  	}
  	result = Clog_sum_exp(log_perms, S, maxval) - log((double)S);

  	// compute log factorials
  	double * log_factorials =(double*) PyMem_Malloc(sizeof(double) * (n+1));
	memset(log_factorials, 0, sizeof(double)*(n+1));

	log_factorials[0]=0.0;
	for (int i = 1; i <= n; ++i)
	{
		log_factorials[i] = log_factorials[i-1] +log((double)(i));
	}

	result = result -log_factorials[n];
  	



	PyMem_Free(log_factorials);
  	return PyFloat_FromDouble(result);

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
    
    double maxval = -NPY_INFINITY;

    for (int i = 0; i < totsize; ++i)
    {
    	if(npy_isnan(array[i])){
			continue;
		}
    	if(array[i]>maxval){
    		maxval = array[i];
    	}
    }

    if(maxval == -NPY_INFINITY){
    	return PyFloat_FromDouble(NPY_NAN);
    }
    double exp_result = 0;



	for (int i = 0; i < totsize; ++i)
	{
		if(npy_isnan(array[i])){
			continue;
		}

		exp_result += exp(array[i] - maxval);
	}

	////printf("res = %f\n", (maxval + log(exp_result)));
	return PyFloat_FromDouble(maxval + log(exp_result));

}

static PyMethodDef permsMethods[] = {
  {"get_log_perms", C_get_log_perms, METH_VARARGS, "get_log_perms(X, t, y, debug)\n\
\n\
Computes log permanents \n\
associated with simulated latent variables.\n\
\n\
Each row of the S x n matrix X contains a random sample of size n from\n\
the data model. If there is only a single covariate, then the\n\
observed data are represented as (t,y), where t is the observed\n\
values of the covariate and y is the vector of indicator variables.\n\
If there are more covariates or the problem is phrased as binary\n\
classification (see Section 5 in [1]), then t is an S x n numpy array,\n\
since the threshold values change in each iteration. The function returns\n\
a vector of log permanents corresponding to each sample in X.\n\
\n\
Parameters \n\
---------- \n\
X : ndarray \n\
    A numpy array of dimension S x n, in \n\
    which each row contains a sample from \n\
    the data model. \n\
t : ndarray\n\
    Either: A flat numpy array of length n\n\
    containing the observed values of\n\
    the covariate, \n\
    Or: A numpy array of dimension S x n (if\n\
    there are several covariates).\n\
y : ndarray\n\
    A flat binary numpy array of length n\n\
    indicating whether x_i<=t_i\n\
    for each i in the observed data.\n\
debug : Boolean\n\
    If true, debug information\n\
    is printed to stdout.\n\
\n\
Returns \n\
------- \n\
ndarray \n\
    Numpy array of log permanents,\n\
    each element associated to \n\
    the corresponding row in X.\n\
    A zero valued permanent is indicated\n\
    by a NaN.\n\
\n\
References\n\
----------\n\
[1] Christensen, D (2023). Inference for Bayesian nonparametric\n\
models with binary response data via permutation counting. \n\
Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.\n"},
  {"log_sum_exp", log_sum_exp, METH_VARARGS, "log_sum_exp(array)\n\
\n\
Computes the log sum exp of an array. \n\
\n\
Given input array = [x_1, ..., x_n], returns \n\
x_* + log(exp(x_1 - x_*) + ... + exp(x_n - x_*)), \n\
where x_* = max(x_1, ... x_n). Ignores entries\n\
with value NaN, as these correspond to vanishing\n\
permanents.\n\
\n\
Parameters \n\
---------- \n\
array : ndarray \n\
    input array \n\
\n\
Returns \n\
------- \n\
float \n"},
{"get_log_ML", C_get_log_ML, METH_VARARGS, "get_log_ML(log_perms, n, debug)\n\
\n\
Computes the log marginal likelihood of the data from the log permanents.\n\
\n\
Given the computed log permanents log_perms, this function\n\
computes the log marginal likelihood using the formula (2.3)\n\
in [1]. It is assumed that there are no repeated trials.\n\
If the data contain repeated trials, then the appropriate log\n\
binomial factor must be added to the output of this function.\n\
\n\
Parameters \n\
---------- \n\
log_perms : ndarray\n\
    A flat numpy array of length n\n\
    containing the computed log permanents,\n\
    where a zero permanent is indicated by \n\
    a NaN.\n\
n : int\n\
    Sample size.\n\
debug : Boolean\n\
    If true, debug information\n\
    is printed to stdout.\n\
\n\
Returns \n\
------- \n\
float \n\
    The estimated log marginal likelihood.\n\
References\n\
----------\n\
[1] Christensen, D (2023). Inference for Bayesian nonparametric\n\
models with binary response data via permutation counting. \n\
Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.\n"},
{"get_log_perms_bioassay", C_get_log_perms_bioassay, METH_VARARGS, "get_log_perms_bioassay(X, levels, successes, trials, debug)\n\
\n\
Computes log permanents associated with simulated latent variables X with\n\
bioassay data.\n\
\n\
Each row of the matrix X contains a random sample of size n from\n\
the data model. The observed data are represented as (levels, \n\
successes, trials), where levels are the different levels at which\n\
trials were conducted, successes is the vector of the number of\n\
successes per level, and trials is the vector of the total number of\n\
trials per level. The function returns a vector of log permanents\n\
corresponding to each sample.\n\
\n\
Parameters \n\
---------- \n\
X : ndarray \n\
    A numpy array of dimension S x n, in \n\
    which each row contains a sample from \n\
    the data model. \n\
levels : ndarray \n\
    A flat numpy array of length n containing\n\
    the levels at which trials were conducted.\n\
successes : ndarray \n\
    A flat numpy array of length n and dtype int32\n\
    containing the number of successful trials at\n\
    each level.\n\
trials : ndarray \n\
    A flat numpy array of length n and dtype int32\n\
    containing the number of trials at each level.\n\
debug : Boolean\n\
    If true, debug information\n\
    is printed to stdout.\n\
\n\
Returns \n\
------- \n\
ndarray \n\
    Numpy array of log permanents,\n\
    each element associated to \n\
    the corresponding row in X.\n\
    A zero valued permanent is indicated\n\
    by a NaN.\n"},
{"get_log_ML_bioassay", C_get_log_ML_bioassay, METH_VARARGS, "get_log_ML_bioassay(log_perms, successes, trials, n, debug)\n\
\n\
Computes the log marginal likelihood of bioassay data from the log permanents. \n\
\n\
Given the computed log permanents log_perms, this function\n\
computes the log marginal likelihood using the formula (2.3)\n\
in [1]. It takes care of repeated trials by adding the appropriate\n\
log binomial factor.\n\
\n\
Parameters \n\
---------- \n\
log_perms : ndarray\n\
    A flat numpy array of length n\n\
    containing the computed log permanents,\n\
    where a zero permanent is indicated by \n\
    a NaN.\n\
successes : ndarray \n\
    A flat numpy array of length n and dtype int32\n\
    contatining the number of successful trials at\n\
    each level.\n\
trials : ndarray \n\
    A flat numpy array of length n and dtype int32\n\
    containing the number of trials at each level.\n\
n : int\n\
    Sample size.\n\
debug : Boolean\n\
    If true, debug information\n\
    is printed to stdout.\n\
\n\
Returns \n\
------- \n\
float \n\
    The estimated log marginal likelihood.\n\
References\n\
----------\n\
[1] Christensen, D (2023). Inference for Bayesian nonparametric\n\
models with binary response data via permutation counting. \n\
Bayesian Analysis, Advance online publication, DOI: 10.1214/22-BA1353.\n"},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef perms = {
  PyModuleDef_HEAD_INIT,
  "perms",
  "Module for computing permanents of block rectangular matrices",
  -1,
  permsMethods
};

PyMODINIT_FUNC PyInit_perms(void)
{
    import_array();
    return PyModule_Create(&perms);
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