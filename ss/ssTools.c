#include "ss.h"
#include "rand.h"
#include "string.h"
#include <math.h>

/*
	Compute the distance of a member from all the members of set
 */
void distance_to_set(SSType *ssParams, Set *set, int set_size, individual *ind){

	double dist = 0;
	for (int i = 0; i < set_size; ++i)
	{
		dist += euclidean_distance(ssParams, &(set->members[i]), ind);
	}
	ind->distance = dist;
}

/*
	Compute the Eculidean Distance between two individuals
 */
double euclidean_distance(SSType *ssParams, individual *ind1, individual *ind2){

	double distance = 0;
	for (int i = 0; i < ssParams->nreal; ++i)
	{
		distance += (ind1->params[i] - ind2->params[i]) * (ind1->params[i] - ind2->params[i]);
	}
	return sqrt(distance);

}

int closest_member(SSType *ssParams, Set *set, int set_size, individual *ind, int ind_index){
	double dist;
	double min;
	int min_index;

	min = euclidean_distance(ssParams, ind, &(set->members[0]));
	min_index = 0;
	
	for (int i = 1; i < set_size; ++i)
	{
		if ( i == ind_index )
			continue;
		
		dist = euclidean_distance(ssParams, ind, &(set->members[i]));
		if (dist < min ){
			min = dist;
			min_index = i;
		}
	}
	return min_index;
}

/*
	Update the best set, containing all the best solution during the optimizations.
 */
void update_bestSet(SSType *ssParams, individual *best){


}

/*
	Return minimum value of an array with its `index`
 */
double min(const double *arr, int length, int *index) {

    double minimum = arr[0];
    for (int i = 1; i < length; ++i) {
        if (minimum > arr[i]) {
            minimum = arr[i];
            *index = i;
        }
    }
    return minimum;
}

/*
	Return the maximum value of an array with its `index`
 */
double max(const double *arr, int length, int *index) {

    double maximum = arr[0];
    for (int i = 1; i < length; ++i) {
        if (maximum < arr[i]) {
            maximum = arr[i];
            *index = i;
        }
    }
    return maximum;
}

void delete_and_shift(SSType *ssParams, Set *set, int set_size, int index_to_delete){

	for (int i = index_to_delete; i < set_size - 1; ++i)
	{
		copy_ind(ssParams, &(set->members[i]), &(set->members[i + 1]));
	}

}


bool is_equal(SSType *ssParams, individual *ind1, individual *ind2){

	bool isEqual = false;
	if ( euclidean_distance(ssParams, ind1, ind2) < ssParams->dist_epsilon )
		isEqual |= 1;	

	return isEqual;
}

/*
	Check if the individual exist in the set
 */
int is_exist(SSType *ssParams, Set *set, int set_size, individual *ind){

	int index = -1;
	for (int i = set_size - 1; i >= 0; --i)
	{
		if ( is_equal(ssParams, &(set->members[i]), ind) ){
			index = i;
			break;
		}

	}

	return index;
}

/* 
	Check if an subset is exist in an subsets_list
*/
bool is_subset_exist(SSType *ssParams, Set *subsets_list, int subsets_list_size, Set *subset, int subset_size, int member_length){
	// TODO: Implement it

	return false;
}

bool is_exist_in_subsets_list(SSType *ssParams, Set *subset){


	return false;
}

/*
	Copy src into the dest.
	Assume that the memory for dest is already allocated.
*/
void copy_ind(SSType *ssParams, individual *dest, individual *src){
	memcpy(dest->params, src->params, ssParams->nreal*sizeof(double));
	dest->cost =  src->cost;
	dest->distance = src->distance;
}

/*
	Check if the selected vector is in the flatzone 
 */
bool is_in_flatzone(SSType *ssParams, Set *set, int set_size, individual *ind){

	bool isInFlatzone = false;

	/* 
		The loop doesn't check the last item since it is the best sol and every good solution in
		comparison to that is in flatzone coverd by that!
	 */
	for (int i = set_size - 1; i > 0; --i)	// set: ref_set
	{
		/* code */
		// printf("%f, %f, %f\n", ind->cost , set->members[i].cost, set->members[i].cost * ( 1 - ssParams->fitness_epsilon) );
		if ( ind->cost > set->members[i].cost * ( 1 - ssParams->fitness_epsilon) )
		{
			isInFlatzone |= 1;
			break;
		}

	}
	return isInFlatzone;
}


void parse_double_row(SSType *ssParams, char *line, double *row){

    int i = 0;
    const char* tok;
    for (tok = strtok(line, "\t"); tok && *tok; i++, tok = strtok(NULL, "\t\n"))
    {
        row[i] = atof(tok);
    }
}


void parse_int_row(SSType *ssParams, char *line, int *row){

    int i = 0;
    const char* tok;
    for (tok = strtok(line, "\t"); tok && *tok; i++, tok = strtok(NULL, "\t\n"))
    {
        row[i] = atoi(tok);
    }
}

void warm_start(SSType *ssParams){
	int i;
    char line[4098];
    printf("Loading the data to perform warm start...\n");
	
	// Read refSet
    FILE* refSetStream = fopen("ref_set_final.csv", "r");
    ssParams->ref_set = (Set *)malloc(sizeof(Set));
	allocate_set_memory(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);
 
    i = 0;
    while (fgets(line, 4098, refSetStream) && (i < ssParams->ref_set_size))
    {
    	double row[ssParams->nreal + 1];
        char* tmp = strdup(line);
        parse_double_row(ssParams, tmp, row);
        for (int j = 0; j < ssParams->nreal; ++j)
        {
        	/* code */
        	ssParams->ref_set->members[i].params[j] = row[j];
        }
        // memcpy(ssParams->ref_set->members[i].params, row, (ssParams->nreal - 1) * sizeof(double));
        ssParams->ref_set->members[i].cost = row[ssParams->nreal];
        // ssParams->ref_set->members[i].distance = row[ssParams->nreal - 1];
        free(tmp);
        i++;
    }
    ssParams->best = (individual *)malloc(sizeof(individual));
	ssParams->best = &(ssParams->ref_set->members[0]);				// The first members of ref_set is always the best

    print_set(ssParams, ssParams->ref_set, ssParams->ref_set_size, ssParams->nreal);

	// Read freqMat
    FILE* freqMatStream = fopen("freq_mat_final.csv", "r");
    i = 0;
    while (fgets(line, 4098, freqMatStream))
    {
    	int row[ssParams->p];
        char* tmp = strdup(line);
        parse_int_row(ssParams, tmp, row);
        memcpy(ssParams->freqs_matrix[i], (int*)row, ssParams->p * sizeof(int));
        free(tmp);
        i++;
    }
    print_int_matrix(ssParams, ssParams->freqs_matrix, ssParams->nreal, ssParams->p);


	// Read probMat
    FILE* probMatStream = fopen("prob_mat_final.csv", "r");
    i = 0;
    while (fgets(line, 4098, probMatStream))
    {
    	double row[ssParams->p];
        char* tmp = strdup(line);
        parse_double_row(ssParams, tmp, row);
        memcpy(ssParams->probs_matrix[i], row, ssParams->p * sizeof(double));
        free(tmp);
        i++;
    }

    print_double_matrix(ssParams, ssParams->probs_matrix, ssParams->nreal, ssParams->p);

}
