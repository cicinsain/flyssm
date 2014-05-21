#ifndef SS_INCLUDED
#define SS_INCLUDED

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>

#include "maternal.h"

#define STATS
#define DEBUG

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

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
	int max_elite;						// = ref_set_size / 2

	int nreal;
	double sol;							// Possible guess or information about the value of minima,
										// for instace, if it is minimization, the 0 would be a guess.
										// It will be used to avoid unneccessary local searches

	double* min_real_var;
	double* max_real_var;
	
	int p;								// The number of sub-regions for boundary intervals
	double **min_boundary_matrix;
	double **max_boundary_matrix;

	individual *best;

	int ref_set_size;					
	Set *ref_set;

	int scatter_set_size;				// Banga's version always has: scatter_set_size = 10 * ref_set_size
	Set *scatter_set;
	
	int pair_size;						// The size of the pairs in subSetItem; normally `2`
	int subsets_list_size;				// The size of the linked list containing the candidates
	Set *subsets_list;

	// Set *recombined_list;			// Basically the subsets_list, first will be recombined
	Set	*candidates_set;				// and then refined to produce candidates_list
	int candidates_set_size;
	
	double dist_epsilon;
	double fitness_epsilon;	 

	int perform_ref_set_regen;
	int ref_set_regen_freq;


	/* Stats */
	int n_refinement;					// Number of local search performed
	int n_ref_set_update;				// Number of substitution in refSet
	int n_duplicates;					// Number of duplicates deteceted in candidateSet
	int n_flatzone_detected;			// Number of flatzone detected during the updating procedure
	int n_function_evals;
	int n_regen;
	
	int **freqs_matrix;					// Frequencies of variables being in sub-regions
	double **probs_matrix;

	int perform_stop_criteria;			// Usually it's not a good idea to activate it!
	double stop_criteria;

	int perform_warm_start;


	int perform_flatzone_detection;


	/* Output */
	char *ref_set_final_filename;
	char *freq_mat_final_filename;
	char *prob_mat_final_filename;


	/* Local Search Parameters */
	int perform_local_search;
	char local_search_method;

	int local_search_1_filter;
	double local_search_f1_criteria;

	int local_search_2_filter;
	double local_search_f2_criteria;

	double min_distance_for_local_search;

	int max_no_improve;
	double step_size;					// Only to perfom take_step local search


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
extern FILE *can_set_history_file;

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
void init_scatter_set(SSType *ssParams, Set *set);
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
void generate_ind_candidate(SSType *ssParams, individual *base, individual *candidate,  double *dists, char type);

// refine.c
bool is_in_flatzone(SSType *ssParams, Set *set, int set_size, individual *ind);
void update_ref_set(SSType *ssParams);
void replace(SSType *ssParams, individual *dest, individual *src, char sort_to_perform);
void compute_Mt(SSType *ssParams, Set *set, int set_size, double **M, int m_row, int m_col);
void re_gen_ref_set(SSType *ssParams, Set *set, int set_size, char type, Input *inp, ScoreOutput *out);

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
bool is_exist_in_subsets_list(SSType *ssParams, individual *ind1, individual *ind2);
bool is_subset_exist(SSType *ssParams, Set *subsets_list, int subsets_list_size, Set *subset, int subset_size, int member_length);
double min(const double *arr, int length, int *index);
double max(const double *arr, int length, int *index);
void delete_and_shift(SSType *ssParams, Set *set, int set_size, int index_to_delete);
int closest_member(SSType *ssParams, Set *set, int set_size, individual *ind, int ind_index);
void copy_ind(SSType *ssParams, individual *src, individual *dest);
void parse_double_row(SSType *ssParams, char *line, double *row);
void parse_int_row(SSType *ssParams, char *line, int *row);
void warm_start(SSType *ssParams);
void matrix_product(SSType *ssParams, double **A, int a_row, int a_col, double **B, int b_row, int b_col, double **P, int p_row, int p_col);

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