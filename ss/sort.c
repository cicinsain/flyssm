#include "ss.h"

/*
	Sort a set based with respect to 'c': cost or 'd':distance
 */
void quick_sort_set(SSType *ssParams, Set *set, int set_size, char key){

	double *numbers;
	numbers = (double *)malloc(set_size * sizeof(double));
	// allocate_ind_memory(ssParams, &(pivot_ind), ssParams->nreal);
			// TODO: It could be as a static variable

	for (int i = 0; i < set_size; ++i)
	{
		if (key == 'c')
			numbers[i] = set->members[i].cost;
		// else if (key == 'd')
		// 	numbers[i] = set->members[i].distance;
	}

	quick_sort(ssParams, set, set_size, numbers, 0, set_size - 1);

	// deallocate_ind_memory(ssParams, &(pivot_ind));
	free(numbers);

}

    individual pivot_ind;								// Might cause problem
/*
	Perform QuickSort in-place, based on the numbers array which could be cost or distance.
 */
void quick_sort(SSType *ssParams, Set *set, int set_size, double* numbers, int left, int right){

	double pivot;
	int l_hold, r_hold;
														// TODO: Memory leak!

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];

	pivot_ind = set->members[left];

	while (left < right)
	{
		while ((numbers[right] >= pivot) && (left < right))
			right--;

		if (left != right)
		{
			numbers[left] = numbers[right];
			set->members[left] = set->members[right];
			// TODO: Check the assigment
			left++;
		}

		while ((numbers[left] <= pivot) && (left < right))
			left++;

		if (left != right)
		{
			numbers[right] = numbers[left];
			set->members[right] = set->members[left];
			// TODO: Check the assignment
			right--;
		}
	}

	numbers[left] = pivot;

	set->members[left] = pivot_ind;

	pivot = left;
	left = l_hold;
	right = r_hold;

	// free(&pivot_ind);

	if (left < pivot)
		quick_sort(ssParams, set, set_size, numbers, left, pivot-1);

	if (right > pivot)
		quick_sort(ssParams, set, set_size, numbers, pivot+1, right);



}

/*
	Remove the last item of the set and put the ind at it's correct place, assuming that 
	the list is already a sorted list.
	So, it's basically a reverse Insertion Sort!
*/
void insertion_sort(SSType *ssParams, Set *set, int set_size, char key){
        
	individual temp;
    // allocate_ind_memory(ssParams, &(temp), ssParams->nreal);
    int j =  set_size - 1;			// The iteration starts from `set_size - 1` since the last
    								// item should be replaced anyway.
 
    while (j > 0 && set->members[j].cost < set->members[j - 1].cost) {

    	temp 				= set->members[j - 1];
    	set->members[j - 1] = set->members[j];
    	set->members[j] 	= temp;

     	j--;

    }

    // free(&temp);
    // deallocate_ind_memory(ssParams, &(temp));

       
}

