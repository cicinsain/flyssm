// #include "ss.h"

// void print_list(SSType *ssParams, subSetItem *head){

// 	subSetItem *current = head;

// 	while( current != NULL ){
// 		printf("%lf, %lf \n", current->pairs[0].cost, current->pairs[1].cost);
// 		current = current->next;
// 		// TODO: Need to be generalized!
// 	}
// }


// void push_at_head(SSType *ssParams, subSetItem **head, individual *newItem){

// 	subSetItem *new_subSet;
// 	new_subSet = malloc(sizeof(subSetItem));
// 	allocate_list_item(ssParams, new_subSet);


// 	// TODO: Make it in a smarter way
// 	for (int i = 0; i < ssParams->pair_size; ++i)
// 	{
// 		/* Copy newPairs into the new_subSet */
// 		new_subSet->pairs[i] = newItem[i];
// 	}

// 	new_subSet->next = *head;
// 	*head = new_subSet;

// 	ssParams->subsets_list_size++;
// }

// int pop_first_item(SSType *ssParams, subSetItem **head){

// 	subSetItem *next_node;

// 	if (*head == NULL)
// 		return -1;		// Unsuccessful 

// 	next_node = (*head)->next;
// 	// TODO: maybe freeing head
// 	*head = next_node;
// 	// TODO: maybe returning the values of the popped item

// 	ssParams->subsets_list_size--;

// 	return 1;	// TODO: Return the value of node
// }

// /*
// 	Allocate memory for a pair in a subSetItem
//  */
// void allocate_list_item(SSType *ssParams, subSetItem *item){

// 	item->pairs = (individual *) malloc( ssParams->pair_size * sizeof(individual));
// 	for (int i = 0; i < ssParams->pair_size; ++i)
// 	{
// 		allocate_ind_memory(ssParams, &(item->pairs[i]), ssParams->pair_size);
// 	}
// 	item->next = NULL;
// }

// void deallocate_list_item(SSType *ssParams, subSetItem *item){


// }