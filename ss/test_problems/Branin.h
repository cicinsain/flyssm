/*Provide dimension of the problem.  Integer >= 1.*/
#define TEST_PROBLEM
#define N 2
/*Provide the diameter of the feasible region.*/
/*Provide maximum number of iterations.*/
#define MAXITER 100*N*N/*500*N */

#define PI 3.14159265359
#define SOL 0.3979

void bounds(double lb[], double ub[]);
int feasible(double x[]);
double objfn(double x[]);

void bounds(double lb[], double ub[])
/*Provide lower and upper bounds for each of N variables.
 Number of bounds is equal to N.*/
{
  lb[0] = -5;
  ub[0] = 10;
  lb[1] = 0;
  ub[1] = 15;
 
}

/*Test feasibility of x.  Return 1 if feasible, 0 if not.*/
int feasible(double x[])

{
	return 1;
}

/*Calculate objective function value of x[].*/
double objfn(double x[])
{
	double sum=0.;
	double a=1., b=5.1/(4*PI*PI), c=5/PI;
	double d=6., e=10., f=1/(8*PI);

	sum = a*pow((x[1]-b*x[0]*x[0]+c*x[0]-d),2) + e*(1-f)*cos(x[0]) + e;

	return (sum);

}
