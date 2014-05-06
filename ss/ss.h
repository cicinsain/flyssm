#ifndef SS_INCLUDED
#define SS_INCLUDED

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>

#include "maternal.h"

#define STATS

#define eul  2.71828182845905
#define pi 3.14159265358979

#ifndef MAX
	#define MAX(x, y) (((x) > (y)) ? (x) : (y))
	#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

typedef enum { false, true } bool;


typedef struct individual
{
	// int index;			// The idea is to use this index for sorting purposes
	double* params;
	double cost;
	double distance; 		// The distance of an individual to a set. Only uses in select_ref_set function.

} individual;

typedef struct Set
{
	individual *members;
	// int size;

} Set;

// typedef struct subSetItem
// {
// 	individual *pairs;
// 	struct subSetItem *next;
// } subSetItem;

typedef struct SSType
{
	int seed;
	int max_iter;
	double step_size;					// Only to perfom take_step local search
	int max_no_improve;
	int max_elite;						// = ref_set_size / 2

	int nreal;
	double sol;							// Possible guess or information about the value of minima,
										// for instace, if it is minimization, the 0 would be a guess.
										// It will be used to avoid unneccessary local searches

	double* min_real_var;
	double* max_real_var;
	
	int p;								// The number of sub regions for boundary intervals
	double **min_boundary_matrix;
	double **max_boundary_matrix;

	individual *best;

	int ref_set_size;					
	Set *ref_set;

	int diverse_set_size;				// Banga's version always has: diverse_set_size = 10 * ref_set_size
	Set *diverse_set;
	
	int pair_size;						// The size of the pairs in subSetItem; normally `2`
	int subsets_list_size;				// The size of the linked list containing the candidates
	Set *subsets_list;

	// Set *recombined_list;			// Basically the subsets_list, first will be recombined
	Set	*candidates_set;				// and then refined to produce candidates_list
	int candidates_set_size;
	
	double dist_epsilon;
	double fitness_epsilon;	 


	/* stats */
	int n_refinement;					// Number of local search performed
	int n_ref_set_update;				// Number of substitution in refSet
	
	int **freqs_matrix;
	double **probs_matrix;

	double local_search_criteria;

	double stop_criteria;

	int perform_warm_start;
	int perform_local_search;
	int perform_flatzone_detection;

	char *ref_set_final_filename;
	char *freq_mat_final_filename;
	char *prob_mat_final_filename;


} SSType;


/*
				Output Files [Global]
 */
extern FILE *ref_set_history_file;
extern FILE *best_sols_history_file;
extern FILE *freqs_matrix_file;
extern FILE *freq_mat_final_file;
extern FILE *prob_mat_final_file;
extern FILE *ref_set_final_file;

/*
				Functions Prototypes
 */

// input.c
void read_input(SSType *ssParams, int argc, char const *argv[]);

// ss.c
void InitSS(Input *inp, SSType *ssParams, char *inname);
void RunSS(Input *inp, SSType *ssParams, char *inname);

// init.c
void init_ssParams(SSType *ssParams);
void init_diverse_set(SSType *ssParams, Set *set);
void diversify(SSType *ssParams, Set *set, int set_size);

// recombine.c
void recombine_subsets_list(SSType *ssParams);
void recombine_subset(SSType *ssParams, Set *subset);
void recombine_ind(SSType *ssParams, individual *ind, double step);

// select.c
void select_subsets_list(SSType *ssParams, Set *set, int set_size/*, subSetItem *subsets_list*/);
void select_ref_set(SSType *ssParams, Set *set, int set_size);
bool select_best(SSType *ssParams);
void init_scatter_set(SSType *ssParams, Set *set);
void init_ref_set(SSType *ssParams);
void generate_candiates(SSType *ssParams);
void generate_ind_candidate(SSType *ssParams, individual *base, individual *candidate,  double *steps, char type);

// refine.c
bool is_in_flatzone(SSType *ssParams, Set *set, int set_size, individual *ind);
void update_ref_set(SSType *ssParams);

// allocate.c
void allocate_ind_memory(SSType *ssParams, individual *ind, int member_length);
void allocate_set_memory(SSType *ssParams, Set *set, int set_size, int member_length);
void deallocate_ind_memory(SSType *ssParams, individual *ind);
void deallocate_set_memory(SSType *ssParams, Set *set, int set_size);
void deallocate_subsets_list_memory(SSType *ssParams);
void deallocate_ssParam(SSType *ssParams);

// sort.c
void quick_sort_set(SSType *ssParams, Set *set, int set_size, char kay);
void quick_sort(SSType *ssParams, Set *set, int set_size, double *numbers, int left, int right);
void insertion_sort(SSType *ssParam, Set *set, int set_size, char key);
// Maybe more function for quick_sort


// local_search.c
void nelder_mead(SSType *ssParams, individual *ind);
void refine_subsets_list(SSType *ssParams, char method, Input *inp, ScoreOutput *out);
void refine_set(SSType *ssParams, Set *set, int set_size, char method, Input *inp, ScoreOutput *out);
void refine_individual(SSType *ssParams, Set *set, int set_size, individual *ind, char method, Input *inp, ScoreOutput *out);
void take_step(SSType *ssParams, double *params, double *new_params);


// ssTools.c
void distance_to_set(SSType *ssParams, Set *set, int set_size, individual *ind);
double euclidean_distance(SSType *ssParams, individual *ind1, individual *ind2);
double dist(double param1, double param2);
void random_ind(SSType *ssParams, individual *ind);
void generate_random_set(SSType *ssParams, Set *set, int set_size);
double rndreal(double low, double high);
void copy_params(SSType *ssParams, individual *ind1, individual *ind2);
void update_bestSet(SSType *ssParams, individual *best);
bool is_equal(SSType *ssParams, individual *ind1, individual *ind2);
int is_exist(SSType *ssParams, Set *set, int set_size, individual *ind);
bool is_exist_in_subsets_list(SSType *ssParams, Set *subset);
bool is_subset_exist(SSType *ssParams, Set *subsets_list, int subsets_list_size, Set *subset, int subset_size, int member_length);
double min(const double *arr, int length, int *index);
double max(const double *arr, int length, int *index);
void delete_and_shift(SSType *ssParams, Set *set, int set_size, int index_to_delete);
int closest_member(SSType *ssParams, Set *set, int set_size, individual *ind, int ind_index);
void copy_ind(SSType *ssParams, individual *src, individual *dest);
void parse_double_row(SSType *ssParams, char *line, double *row);
void parse_int_row(SSType *ssParams, char *line, int *row);
void warm_start(SSType *ssParams);

// report.c
void init_report_files(SSType *ssParams);
void write_set(SSType *ssParams, Set *set, int set_size, int member_length, FILE *fpt, int iter, char mode);
void write_ind(SSType *ssParams, individual *ind,  int member_length, FILE *fpt, int iter, char mode);
void print_set(SSType *ssParams, Set *set, int set_size, int member_length);
void print_ind(SSType *ssParams, individual *ind, int member_length);
void print_subsets_list(SSType *ssParams);
void print_double_matrix(SSType *ssParams, double **matrix, int row, int col);
void print_int_matrix(SSType *ssParams, int **matrix, int row, int col);
void write_int_matrix(SSType *ssParams, int **matrix, int row, int col, FILE *ftp, int iter, char mode);
void write_double_matrix(SSType *ssParams, double **matrix, int row, int col, FILE *ftp, int iter, char mode);
void loadBar(int x, int n, int r, int w);
void write_params_to_fly_output_standard(SSType *ssParams, Input *inp, char *inname);


// stats.c
void update_frequency_matrix(SSType *ssParams, individual *ind);

// problemdef.c
void init_sample_params(SSType *ssParams);
void sample_obj_func(SSType *ssParams, individual *ind);
double objective_function(double *s, SSType *ssParams, Input *inp, ScoreOutput *out);

// evaluate.c
void evaluate_ind(SSType *ssParams, individual *ind, Input *inp, ScoreOutput *out);
void evaluate_set(SSType *ssParams, Set *set, int set_size, Input *inp, ScoreOutput *out);

/*
// linkedlist.c
void allocate_list_item(SSType *ssParams, subSetItem *item);
void deallocate_list_item(SSType *ssParams, subSetItem *item);
void print_list(SSType *ssParams, subSetItem *head);
int pop_first_item(SSType *ssParams, subSetItem **head);
void push_at_head(SSType *ssParams, subSetItem **head, individual *newItem);

*/

#endif