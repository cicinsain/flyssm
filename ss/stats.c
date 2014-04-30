#include "ss.h"

void update_frequency_matrix(SSType *ssParams, individual *ind){

	for (int i = 0; i < ssParams->nreal; ++i)
	{
		for (int j = 0; j < ssParams->p; ++j)
		{
			if (ind->params[i] > ssParams->min_boundary_matrix[i][j] && ind->params[i] < ssParams->max_boundary_matrix[i][j])
			{
				ssParams->freqs_matrix[i][j]++;
				break;
			}
		}
	}
}