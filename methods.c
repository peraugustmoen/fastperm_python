#include "header.h"




// reverse top trim
void reverse_tt(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
				int * k){


	// CHECK IF LOG IS NEG!!
	
	for (int i = 0; i < (n+1)*(n+1); ++i)
	{
		new_log_subperms[i]=NPY_NAN;
	}

	for (int r = amount; r <= alpha[0]; ++r)
	{
		double log_lower = log_factorials[alpha[0] - r];
		double log_upper = log_factorials[amount + alpha[0] - r];

		for (int s = 0; s <= alpha[*k-1]; ++s)
		{
			if(npy_isnan(old_log_subperms[cord_spec(r-amount,s,(n+1) )])){
				new_log_subperms[cord_spec(r,s, (n+1))] = NPY_NAN;
				continue;
			}
			if(old_log_subperms[cord_spec(r-amount,s,(n+1) )]<0){
				new_log_subperms[cord_spec(r,s, (n+1))] = NPY_NAN;
				continue;
			}

			new_log_subperms[cord_spec(r,s, (n+1))] = log_upper - log_lower +
			 	old_log_subperms[cord_spec(r-amount,s, (n+1))];
		}

	}


	// change alpha beta 
	
	beta[0] = amount;

	*m = (*m) + amount;

	return;
	
}



// Reverse bottom trim
void reverse_bt(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
				int * k){

	for (int i = 0; i < (n+1)*(n+1); ++i)
	{
		new_log_subperms[i]=NPY_NAN;
	}

	for (int r = 0; r <= alpha[0]; ++r)
	{

		for (int s = amount; s <= alpha[*k-1]; ++s)
		{
			if(npy_isnan(old_log_subperms[cord_spec(r,s-amount, (n+1) )])){
				new_log_subperms[cord_spec(r,s, (n+1))] = NPY_NAN;
				continue;
			}

			if(old_log_subperms[cord_spec(r,s-amount, (n+1))]<0){
				new_log_subperms[cord_spec(r,s, (n+1))] = NPY_NAN;
				continue;
			}
			double log_upper = log_factorials[(alpha[(*k)-1] - s + amount)];
			double log_lower = log_factorials[alpha[(*k)-1] - s];



			new_log_subperms[cord_spec(r,s, (n+1))] = log_upper - log_lower + old_log_subperms[cord_spec(r,s-amount, (n+1))];
		}

	}

	// change alpha beta osv osv
	
	gamma[(*k) - 2] = amount;
	*m = (*m) + amount;

	return;
	
}

// Reverse bottom split
void reverse_bs(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
				double * temp_vec, int * k){

	for (int i = 0; i < (n+1)*(n+1); ++i)
	{
		new_log_subperms[i]=NPY_NAN;
	}

	////printf("HEIIII\n");
	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3a=0;
	double log_lower3b=0;


	double maxval = -1;
	int lstart = 0;

	
	for (int r = 0; r <= alpha[0]; ++r)
	{
		for (int s = MAX(0, (*m + amount) +alpha[0] + alpha[(*k)-1]-n-r); s<= alpha[(*k)-1]; ++s)
		{
			maxval=-1;
			lstart = MAX(s-amount,0);
			for (int l = lstart; l <= s; ++l)
			{


				if(r == 1 && s == 1){

				}
				if(npy_isnan(old_log_subperms[cord_spec(r,l,(n+1))])){
					temp_vec[l-lstart] = NPY_NAN;
					continue;
				}

				if(old_log_subperms[cord_spec(r,l,(n+1))] < 0 ){
					temp_vec[l-lstart] = NPY_NAN;
					continue;
				}
				




				log_upper1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - ((*m) - r - l)];
				log_lower1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - (*m + amount) + r +s];
				log_upper2 = log_factorials[alpha[(*k)-1] -l  ];
				log_lower2 = log_factorials[alpha[(*k) -1] -s];
				log_upper3 = log_factorials[ amount ];
				log_lower3a = log_factorials[ amount - s + l ];
				log_lower3b = log_factorials[ s - l ];

				temp_vec[l-lstart] = log_upper1 + log_upper2 + log_upper3 - log_lower1 - log_lower2 -
							   log_lower3a - log_lower3b + old_log_subperms[cord_spec(r,l,(n+1))];

				if(r==1 && s==1){
					////printf("log_upper1 = %f\n", log_upper1);
					////printf("log_upper2 = %f\n", log_upper2);
					////printf("log_upper3 = %f\n", log_upper3);
					////printf("log_lower1 = %f\n", log_lower1);
					////printf("log_lower2 = %f\n", log_lower2);
					////printf("log_lower3a = %f\n", log_lower3a);
					////printf("log_lower3b = %f\n", log_lower3b);

				}

				if(temp_vec[l-lstart] > maxval){
					maxval = temp_vec[l-lstart];
				}



			}


			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, s-lstart+1, maxval);
			
		}
	}

	// change alpha beta 
	
	gamma[0] = amount;
	*m = (*m) + amount;

	return;
	
}

// Reverse top split
void reverse_ts(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
				double * temp_vec, int * k){

	for (int i = 0; i < (n+1)*(n+1); ++i)
	{
		new_log_subperms[i]=NPY_NAN;
	}

	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3a=0;
	double log_lower3b=0;


	double maxval = -1;
	int lstart = 0;

	for (int r = 0; r <= alpha[0]; ++r)
	{
		for (int s = MAX((*m+amount+alpha[0]+alpha[(*k)-1] -n -r),0); s<= alpha[(*k)-1]; ++s)
		{
			maxval=-1;
		  lstart = MAX(0, r - amount);
			for (int l = lstart; l <= r; ++l)
			{
				if(npy_isnan(old_log_subperms[cord_spec(l,s,(n+1))])){
					temp_vec[l-lstart] = NPY_NAN;
					continue;
				}

				if(old_log_subperms[cord_spec(l,s,(n+1))] < 0 ){
					temp_vec[l-lstart] = NPY_NAN;
					continue;
				}
				


				log_upper1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - ((*m)   - l - s )];
				log_lower1 = log_factorials[n - alpha[0] - alpha[(*k) -1] - (*m + amount) + r +s];
				log_upper2 = log_factorials[alpha[0] - l];
				log_lower2 = log_factorials[alpha[0] - r];
				log_upper3 = log_factorials[ amount ];
				log_lower3a = log_factorials[ amount - r + l ];
				log_lower3b = log_factorials[ r-l];

				temp_vec[l-lstart] = log_upper1 + log_upper2 + log_upper3 - log_lower1 - log_lower2 -
							   log_lower3a - log_lower3b + old_log_subperms[cord_spec(l,s,(n+1))];


				if(temp_vec[l-lstart] > maxval){
					maxval = temp_vec[l-lstart];
				}



			}


			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, r - lstart+1, maxval);
			
		}
	}

	// change alpha beta
	beta[(*k) - 2] = amount;
	*m = (*m) + amount;

	return;
	
}


// Reverse left merge
void reverse_lm(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
				double * temp_vec, int * k){



	for (int i = 0; i < (n+1)*(n+1); ++i)
	{
		new_log_subperms[i]=NPY_NAN;
	}

	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3=0;
	double log_upper4=0;
	double log_lower4a=0;
	double log_lower4b=0;


	double maxval = -1;

	log_upper1 = log_factorials[amount2];
	log_upper2 = log_factorials[amount];
	log_upper3 = log_factorials[amount + amount2];

	for (int r = 0; r <= (amount); ++r)
	{
		log_lower2 = log_factorials[amount - r];
		for (int s = 0; s<= alpha[(*k)-1]; ++s)
		{
			maxval=-1;
			for (int l = r; l <= (r+amount2); ++l)
			{
				if(npy_isnan(old_log_subperms[cord_spec(l,s,(n+1))])){
					temp_vec[l-r] = NPY_NAN;
					continue;
				}

				if(old_log_subperms[cord_spec(l,s,(n+1))] < 0 ){
					temp_vec[l-r] = NPY_NAN;
					continue;
				}


				log_lower1 = log_factorials[amount2 - l +r];
				log_lower3 = log_factorials[ amount +amount2 -l];
				log_upper4 = log_factorials[l];
				log_lower4a = log_factorials[l-r];
				log_lower4b = log_factorials[r];

			

				temp_vec[l-r] = log_upper1 + log_upper2 - log_upper3 + log_upper4
								- log_lower1 - log_lower2 + log_lower3
								- log_lower4a -log_lower4b + old_log_subperms[cord_spec(l,s,(n+1))];


				if(temp_vec[l-r] > maxval){
					maxval = temp_vec[l-r];
				}


			}


			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, (amount2+1), maxval);
			
		}
	}

	// change alpha beta
	
	for (int i = (*k); i>=2; --i)
	{
		alpha[i] = alpha[i-1];
	}
	alpha[0] = amount;
	alpha[1] = amount2;

	for (int i = ((*k)-1); i>=1; --i)
	{
		beta[i] = beta[i-1];
		gamma[i] = gamma[i-1];
	}

	beta[0] = 0;
	gamma[0]= 0;


	
	(*k)++;

	return;
	
}



// Reverse right merge
void reverse_rm(double * old_log_subperms, double * new_log_subperms, int amount, int amount2,
				int * alpha, int * beta, int * gamma, double * log_factorials, int n, int * m, 
				double * temp_vec, int * k){

	for (int i = 0; i < (n+1)*(n+1); ++i)
	{
		new_log_subperms[i]=NPY_NAN;
	}



	double log_upper1=0;
	double log_upper2=0;
	double log_upper3=0;
	double log_lower1=0;
	double log_lower2=0;
	double log_lower3=0;
	double log_upper4=0;
	double log_lower4a=0;
	double log_lower4b=0;


	double maxval = -1;

	log_upper1 = log_factorials[amount];
	log_upper2 = log_factorials[amount2];
	log_upper3 = log_factorials[amount + amount2];

	for (int r = 0; r <= alpha[0]; ++r)
	{
		log_lower2 = log_factorials[amount - r];
		for (int s = 0; s<= (amount2); ++s)
		{
			maxval=-1;
			for (int l = s; l <= (amount+s); ++l)
			{
				if(npy_isnan(old_log_subperms[cord_spec(r,l,(n+1))])){
					temp_vec[l-s] = NPY_NAN;
					continue;
				}

				if(old_log_subperms[cord_spec(r,l,(n+1))] < 0 ){
					temp_vec[l-s] = NPY_NAN;
					continue;
				}


				log_lower1 = log_factorials[amount - l +s];
				log_lower2 = log_factorials[amount2 - s];
				log_lower3 = log_factorials[amount +amount2 -l];
				log_upper4 = log_factorials[l];
				log_lower4a = log_factorials[l-s];
				log_lower4b = log_factorials[s];

			

				temp_vec[l-s] = log_upper1 + log_upper2 - log_upper3 + log_upper4
								- log_lower1 - log_lower2 + log_lower3
								- log_lower4a -log_lower4b + old_log_subperms[cord_spec(r,l,(n+1))];


				if(temp_vec[l-s] > maxval){
					maxval = temp_vec[l-s];
				}


			}

			// log sum exp 

			new_log_subperms[cord_spec(r,s,(n+1))] = Clog_sum_exp(temp_vec, (amount+1), maxval);
			
		}
	}

	// change alpha beta
	
	alpha[(*k)-1] = amount;
	alpha[(*k)] = amount2;

	beta[(*k)-1] = 0;
	gamma[(*k)-1] = 0;

	
	(*k)++;

	return;
}


void get_reduced_log_subperms(double * new_log_subperms, int * alpha, int * beta, int * gamma, 
				double * log_factorials, int n, int * m, int * k){
	// case 1': 
	if((*k)==2 && beta[0]==0 && gamma[0]==0){

		for (int r = 0; r < (n+1); ++r)
		{
			for (int s = 0; s < (n+1); ++s)
			{
				if(r+s==(*m)){
					new_log_subperms[cord_spec(r,s,(n+1))] = log_factorials[(*m)]- log_factorials[r] 
					- log_factorials[(*m)-r] + log_factorials[alpha[0]] - log_factorials[alpha[0]-r] 
					+ log_factorials[alpha[1]] - log_factorials[alpha[1]-s];
				}
				else{
					new_log_subperms[cord_spec(r,s,(n+1))] = NPY_NAN;
				}
			}
		}
	}

	// case 2: 

	else if((*k)==2 && beta[0]==0 && gamma[0]==(*m)){
		for (int r = 0; r < (n+1); ++r)
		{
			for (int s = 0; s < (n+1); ++s)
			{
				
				new_log_subperms[cord_spec(r,s,(n+1))] = NPY_NAN;
				
			}
		}
		new_log_subperms[cord_spec(0,(*m), (n+1))] = log_factorials[alpha[1]] - log_factorials[alpha[1]-(*m)];

	}

	// case 3:
	else if((*k)==2 && beta[0]==(*m) && gamma[0]==0){
		for (int r = 0; r < (n+1); ++r)
		{
			for (int s = 0; s < (n+1); ++s)
			{
				
				new_log_subperms[cord_spec(r,s,(n+1))] = NPY_NAN;
				
			}
		}
		new_log_subperms[cord_spec((*m),0,(n+1))] = log_factorials[alpha[0]] - log_factorials[alpha[0]-(*m)];
	}

	// case 4:
	else if((*k)==3 && beta[0]==0 && beta[1]==(*m) && gamma[0] == (*m) && gamma[1] == 0){
		for (int r = 0; r < (n+1); ++r)
		{
			for (int s = 0; s < (n+1); ++s)
			{
				
				new_log_subperms[cord_spec(r,s,(n+1))] = NPY_NAN;
				
			}
		}
		new_log_subperms[cord_spec(0,0,(n+1))] = log_factorials[alpha[1]] - log_factorials[alpha[1]-(*m)];
	}
	else{
	}
}

char check_if_reduced(int * alpha, int * beta, int * gamma, 
					  int * m, int * k){

	char reduced = 0;

	// case 1: 
	if((*k)==2 && beta[0]==0 && gamma[0]==0){
		reduced =1;
	}

	// case 2: 

	else if((*k)==2 && beta[0]==0 && gamma[0]==(*m)){
		reduced = 1;

	}

	// case 3:
	else if((*k)==2 && beta[0]==(*m) && gamma[0]==0){
		reduced = 1;
	}


	// case 4:
	else if((*k)==3 && beta[0]==0 && beta[1]==(*m) && gamma[0] == (*m) && gamma[1] == 0){
		reduced = 1;
	}

	return reduced;
}


int reduction(int * alpha, int * beta, int * gamma, int *m, int n, int *k, int * history,
			   int * amount_history, int * history_len, int debug){

	// returns 0 if success, otherwise -1

	
	if (history==NULL || amount_history==NULL)
	{
		return -1;
	}
	
	// assume history and amount_history are of size 3n and 6n, respectively
	// amount_history assumed 3n \times 2 matrix, column wise stored

	// history codes:
	// 0 = top trim
	// 1 = bottom split
	// 2 = left merge
	// 3 = bottom trim
	// 4 = top split
	// 5 = right merge

	memset(history, 0, sizeof(int)*3*n );
	memset(amount_history, 0, sizeof(int)*6*n );
	*history_len =0;


	char reduced = check_if_reduced(alpha, beta, gamma, 
					  m, k);

	int historycounter = 0;

	int cc = 0;

	if(!reduced){
		while(1){
			if(alpha[0]>n){
				printf("ERROR! alpha[0]>n\n");
				printf("=====\n");
				printf("alpha:\n");
				print_int_vector(n, alpha);
				printf("beta:\n");
				print_int_vector(n, beta);
				printf("gamma:\n");
				print_int_vector(n, gamma);
				printf("m:\n");
				print_int_vector(1, m);
				printf("k:\n");
				print_int_vector(1, k);
				printf("=====\n");
				return -1;
			}
			if((*k)<0){
				printf("ERROR! k<0\n");
				printf("=====\n");
				printf("alpha:\n");
				print_int_vector(n, alpha);
				printf("beta:\n");
				print_int_vector(n, beta);
				printf("gamma:\n");
				print_int_vector(n, gamma);
				printf("m:\n");
				print_int_vector(1, m);
				printf("k:\n");
				print_int_vector(1, k);
				printf("=====\n");
				return -1;
			}
			// if top trim possible
			if((*m)>beta[0] && beta[0] > 0){
				history[historycounter] = 0;
				amount_history[2*historycounter] = beta[0];
				(*m) = (*m)-beta[0];
				beta[0]=0;
				historycounter++;
				if(debug){
					printf("TOP TRIM ->\n");
					//print_int_vector(n, beta);
					printf("=====\n");
					printf("alpha:\n");
					print_int_vector(n, alpha);
					printf("beta:\n");
					print_int_vector(n, beta);
					printf("gamma:\n");
					print_int_vector(n, gamma);
					printf("m:\n");
					print_int_vector(1, m);
					printf("k:\n");
					print_int_vector(1, k);
					printf("=====\n");
				}
				
				if(check_if_reduced(alpha, beta, gamma, 
						  m, k)){break;}

			}
			// if bottom split possible
			else if((*m)>gamma[0] && gamma[0] > 0){
				history[historycounter] = 1;
				amount_history[2*historycounter] = gamma[0];
				(*m) = (*m) - gamma[0];
				gamma[0] = 0;
				historycounter++;
				if(debug){
					printf("BOTTOM SPLIT ->\n");
					//print_int_vector(n, beta);
					printf("=====\n");
					printf("alpha:\n");
					print_int_vector(n, alpha);
					printf("beta:\n");
					print_int_vector(n, beta);
					printf("gamma:\n");
					print_int_vector(n, gamma);
					printf("m:\n");
					print_int_vector(1, m);
					printf("k:\n");
					print_int_vector(1, k);
					printf("=====\n");
				}
				
				if(check_if_reduced(alpha, beta, gamma, 
						  m, k)){break;}
			}

			// if bottom left merge
			else if(beta[0]==0 && gamma[0]==0){
				history[historycounter] = 2;
				amount_history[2*historycounter] = alpha[0];
				amount_history[2*historycounter+1] = alpha[1];

				//shift alpha:
				alpha[0] = alpha[0]+ alpha[1];
				for (int i = 1; i < (*k-1); ++i)
				{
					alpha[i] = alpha[i+1];
				}

				for (int i = 0; i < (*k-2); ++i)
				{
					beta[i] = beta[i+1];
					gamma[i] = gamma[i+1];
				}
				(*k)--;
				historycounter++;
				if(debug){
					printf("LEFT MERGE ->\n");
					//print_int_vector(n, beta);
					printf("=====\n");
					printf("alpha:\n");
					print_int_vector(n, alpha);
					printf("beta:\n");
					print_int_vector(n, beta);
					printf("gamma:\n");
					print_int_vector(n, gamma);
					printf("m:\n");
					print_int_vector(1, m);
					printf("k:\n");
					print_int_vector(1, k);
					printf("=====\n");	
				}
				if(check_if_reduced(alpha, beta, gamma, 
						  m, k)){break;}
			}		
			
			// if bottom bottom trim possible
			else if((*m) > gamma[(*k)-2] && gamma[(*k)-2]>0){
				history[historycounter] = 3;
				amount_history[2*historycounter] = gamma[(*k)-2];

				(*m) = (*m) - gamma[(*k)-2];
				gamma[(*k)-2]=0;
				historycounter++;
				if(debug){
					printf("BOTTOM TRIM ->\n");
					//print_int_vector(n, beta);
					printf("=====\n");
					printf("alpha:\n");
					print_int_vector(n, alpha);
					printf("beta:\n");
					print_int_vector(n, beta);
					printf("gamma:\n");
					print_int_vector(n, gamma);
					printf("m:\n");
					print_int_vector(1, m);
					printf("k:\n");
					print_int_vector(1, k);
					printf("=====\n");

				}
				if(check_if_reduced(alpha, beta, gamma, 
						  m, k)){break;}
			}

			// if top split possible
			else if( (*m)> beta[(*k)-2] && beta[(*k)-2]>0 ){
				history[historycounter] = 4;
				amount_history[2*historycounter] = beta[(*k)-2];
				(*m) = (*m) - beta[(*k)-2];
				beta[(*k)-2]=0;
				historycounter++;

				if(debug){
					printf("TOP SPLIT ->\n");
					printf("=====\n");
					printf("alpha:\n");
					print_int_vector(n, alpha);
					printf("beta:\n");
					print_int_vector(n, beta);
					printf("gamma:\n");
					print_int_vector(n, gamma);
					printf("m:\n");
					print_int_vector(1, m);
					printf("k:\n");
					print_int_vector(1, k);
					printf("=====\n");

				}
				if(check_if_reduced(alpha, beta, gamma, 
						  m, k)){break;}
			}

			// if right merge possible
			else if(beta[(*k)-2]==0 && gamma[(*k)-2]==0 ){
				history[historycounter] = 5;
				amount_history[2*historycounter] = alpha[(*k)-2];
				amount_history[2*historycounter+1] = alpha[(*k)-1];


				alpha[(*k)-2] = alpha[(*k)-2] + alpha[(*k)-1];
				(*k)--;

				historycounter++;

				if(debug){
					printf("RIGHT MERGE ->\n");
					//print_int_vector(n, beta);
					printf("=====\n");
					printf("alpha:\n");
					print_int_vector(n, alpha);
					printf("beta:\n");
					print_int_vector(n, beta);
					printf("gamma:\n");
					print_int_vector(n, gamma);
					printf("m:\n");
					print_int_vector(1, m);
					printf("k:\n");
					print_int_vector(1, k);
					printf("=====\n");
				}
				if(check_if_reduced(alpha, beta, gamma, 
						  m, k)){break;}
			}
			else{
				printf("ERROR!");
				return -1;
				break;
			}		

			cc++;


		}
	}

	
	*history_len = historycounter;
	return 0;
}

double * reverse_reduction(double * old_log_subperms, double * new_log_subperms, int * alpha, 
					   int * beta, int * gamma, int *m, int n, int *k, int * history,
			           int * amount_history, int * history_len, double * log_factorials, 
			           double * temp_vec){
	// INPUT OLD LOG SUBPERMS = the result of the reduced case
	
	// returns pointer to the new log subperms

	double * tmp =NULL;

	int amount = 0;
	int amount2 = 0;

	int num_nonzeros = 0;
	for (int i =((*history_len)-1); i>=0; --i)
	{
		num_nonzeros = 0;

		for (int i = 0; i < (n+1)*(n+1); ++i)
		{
			if(old_log_subperms[i] >0){
				num_nonzeros++;
			}
		}
		if(history[i]==0){
			// reverse top trim
			amount = amount_history[2*i];
			reverse_tt(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				k);
		}
		else if(history[i] ==1){
			// reverse bottom split
			amount = amount_history[2*i];
			reverse_bs(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				temp_vec, k);
		}

		else if(history[i]==2){
			// reverse left merge
			amount = amount_history[2*i];
			amount2 = amount_history[2*i+1];
			reverse_lm(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				temp_vec, k);
		}

		else if(history[i]==3){
			// reverse bottom trim
			amount = amount_history[2*i];
			reverse_bt(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				k);
		}
		else if(history[i]==4){
			// reverse top split
			amount = amount_history[2*i];
			reverse_ts(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				temp_vec,k);
		}
		else if(history[i]==5){
			// reverse right merge
			amount = amount_history[2*i];
			amount2 = amount_history[2*i+1];
			reverse_rm(old_log_subperms, new_log_subperms, amount, amount2,
				alpha,  beta, gamma,log_factorials, n, m, 
				temp_vec,k);
		}
		

		tmp  = old_log_subperms;
		old_log_subperms = new_log_subperms;
		new_log_subperms = tmp;
	}
	

	return old_log_subperms;
}



