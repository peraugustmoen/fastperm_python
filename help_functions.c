#include "header.h"


size_t hash(pair p, int table_size) {
    
    int t[2];
    t[0] = p.x;
    t[1] = p.y;

    // apply xxhash:
    XXH64_hash_t hash = XXH64(t, sizeof(int)*2, 1337);
 
    return( hash % table_size);

}



double lookup_logperm(dictionary * dict, pair p){
	size_t h = hash(p, (*dict).table_size);
	while ((*dict).table[h] >= 0) {
        if ((*dict).array[(*dict).table[h]].x == p.x && (*dict).array[(*dict).table[h]].y == p.y) {
            return (*dict).value_array[(*dict).table[h]];
        }
        h = (h + 1) % (*dict).table_size;
    }

    return -1;
}


void add_to_dictionary(dictionary * dict, pair p, double val){
	// this function requires that the key (the pair)
	// is not already in the dictionary dict


	if ((RATIO*(*dict).used_len +1)> (*dict).table_size)
	{	
		
		expand_dictionary(dict);
	}

	size_t used_len = (*dict).used_len;

	size_t h = hash(p, (*dict).table_size);
    while ((*dict).table[h] >= 0) {
            h = (h + 1) % (*dict).table_size;  // Linear probing to find the next empty slot
    }
    (*dict).table[h] = used_len;
	(*dict).value_array[used_len] = val;
	(*dict).array[used_len] = p;
	(*dict).used_len++;

}

void update_dict(pair p, double value, dictionary * dict){
	// if the key (the pair) is not in the dictionary, 
	// then this key is added with value set to value.
	// if the key exists, then we log-sum-exp the new
	// and the old value



	// check if (r,s) is already key in dict:
	
	size_t h = hash(p, (*dict).table_size);
	int index = -1;
	while ((*dict).table[h] >= 0) {

        if ((*dict).array[(*dict).table[h]].x == p.x && (*dict).array[(*dict).table[h]].y == p.y) {
            index = (*dict).table[h];
            break;
        }
        h = (h + 1) % (*dict).table_size;
    }

   
    if(index ==-1 || ( (size_t) index) >=(*dict).used_len){
    	// the key (pair) is not contained in the dictionary, 
    	// adding here:

    	// check if dict needs to be expanded:
    	if (RATIO*((*dict).used_len +1)> (*dict).table_size)
		{	
			expand_dictionary(dict);
		}

		size_t used_len = (*dict).used_len;
	    (*dict).table[h] = used_len;
		(*dict).value_array[used_len] = value;
		(*dict).array[used_len] = p;
		(*dict).used_len++;
		return;
    }

    double existing_val = (*dict).value_array[index];

    // log sum exp:
    
    if(existing_val>value){
    	(*dict).value_array[index] += log(1.0 + exp(value - existing_val));
    }else{
    	(*dict).value_array[index] = value + log(1.0+ exp(existing_val - value));
    }



}



dictionary * init_dictionary(size_t init_size){
	// create dictionary and alloc memory etc

	dictionary * dict = (dictionary *) PyMem_Calloc(1, sizeof(dictionary));
	(*dict).array = (pair *) PyMem_Calloc(init_size,  sizeof(pair));
	(*dict).used_len = 0;
	(*dict).table_size = init_size;
	(*dict).value_array = (double*) PyMem_Calloc(init_size, sizeof(double));
	(*dict).table = (int*) PyMem_Calloc(init_size, sizeof(int));

	for (size_t i = 0; i < (*dict).table_size; i++) {
        (*dict).table[i] = -1;  // Initialize all entries to -1 (empty)
    }
	return dict;
}

void free_dictionary(dictionary * dict){
	PyMem_Free((*dict).array);
	PyMem_Free((*dict).value_array);
	PyMem_Free(dict);
}

void expand_dictionary(dictionary * dict){
	// expand the size of a dictionary. the function first tries to increase the memory size by
	// PyMem_Realloc(). this this fails then new memory is allocated. 

	pair * newarray = PyMem_Realloc((*dict).array, (size_t) 2*(*dict).table_size*sizeof(pair));
	if(newarray==NULL){
		newarray= (pair *) PyMem_Calloc( 2 * ((*dict).table_size),  sizeof(pair));
		memcpy(newarray, (*dict).array, sizeof(pair)*(*dict).used_len);
		PyMem_Free((*dict).array);
	}
	(*dict).array = newarray;
	
	
	double * newvalue_array = (double *) PyMem_Realloc((*dict).value_array, (size_t) 2 * ((*dict).table_size)*sizeof( double));
	if(newvalue_array==NULL){
		newvalue_array= (double *) PyMem_Calloc( 2 * ((*dict).table_size),  sizeof(double));
		memcpy(newvalue_array, (*dict).value_array, sizeof(double)*(*dict).used_len);
		PyMem_Free((*dict).value_array);
	}
		
	(*dict).value_array = newvalue_array;
	

	int * newtable = (int *) PyMem_Realloc((*dict).table, (size_t) 2 * ((*dict).table_size) *sizeof(int));
	if(newtable==NULL){
		newtable= (int *) PyMem_Calloc( 2 * ((*dict).table_size),  sizeof(int));
		memcpy(newtable, (*dict).table, sizeof(int)*(*dict).table_size);
		PyMem_Free((*dict).table);
		
	}
	(*dict).table = newtable;
	
	(*dict).table_size = 2*(*dict).table_size;


	for (size_t  i = 0; i < (*dict).table_size; i++) {
        (*dict).table[i] = -1;  // Initialize all entries to -1 (empty)
    }
    for (size_t i = 0; i < (*dict).used_len; i++) {
        size_t h = hash((*dict).array[i], (*dict).table_size);
        while ((*dict).table[h] >= 0) {
            h = (h + 1) % (*dict).table_size;  // Linear probing to find the next empty slot
        }
        (*dict).table[h] = i;
    }

}

void nullset_dictionary(dictionary * dict){
	(*dict).used_len = 0;
	for (size_t i = 0; i < (*dict).table_size; i++) {
        (*dict).table[i] = -1;  
    }
}




void print_sparse_matrix(dictionary * matrix) {

    // Print out the matrix with aligned columns
    for (size_t i = 0; i < (*matrix).used_len; i++) {
    	printf("(%d, %d) = %f\n", (*matrix).array[i].x, (*matrix).array[i].y, (*matrix).value_array[i]);

    }
}


void print_matrix(int rows, int cols, double * matrix) {
	return;

    // Print out the matrix with aligned columns
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%5f ", matrix[cord_spec(i,j,rows)]);
        }
    }
}

void print_int_vector(int len,  int * vec) {

    // Print out the matrix with aligned columns
    printf("(");
    for (int i = 0; i < (len); i++) {
        printf("%d ", vec[i]);
    }
    printf(")\n");
}
void print_float_vector(int len,  double * vec) {
	/*return;*/

    // Print out the matrix with aligned columns
    printf("(");
    for (int i = 0; i < (len); i++) {
        printf("%f ", vec[i]);
    }
    printf(")\n");
}



double Clog_sum_exp(double * array, int len, double maxval){

	// ignore NaN and values

	if(maxval<0){
		return NPY_NAN;
	}

	double exp_result = 0;



	for (int i = 0; i < len; ++i)
	{
		if(npy_isnan(array[i])){
			continue;
		}

		if(array[i]<0){
			continue;
		}

		exp_result += exp(array[i] - maxval);
	}

	return (maxval + log(exp_result));
}

double Csparse_log_sum_exp(dictionary * dict){

	// ignore neg values

	double maxval = -1;

	for (size_t z = 0; z < (*dict).used_len; ++z)
	{
		if((*dict).value_array[z]> maxval){
			 maxval = (*dict).value_array[z];
		}
	}

	if(maxval<0){
		return NPY_NAN;
	}

	double exp_result = 0;



	for (size_t z = 0; z < (*dict).used_len; ++z)
	{

		exp_result += exp((*dict).value_array[z] - maxval);
	}

	////printf("res = %f\n", (maxval + log(exp_result)));
	return (maxval + log(exp_result));
}

