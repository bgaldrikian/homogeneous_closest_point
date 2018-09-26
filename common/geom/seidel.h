/*
	Implementation of Seidel's LP algorithm adapted from

	http://www.cs.sunysb.edu/~algorith/implement/linprog/distrib/source.a

	See copyright notice below.

	Prefix SLP stands for Seidel Linear Program
*/

/* 
 *
 * Copyright (c) 1990 Michael E. Hohmeyer,
 *       hohmeyer@icemcfd.com
 * Permission is granted to modify and re-distribute this code in any manner
 * as long as this notice is preserved.  All standard disclaimers apply.
 *       
 */


/* Define to test validity */
//#define SLP_CHECK


/* SLP_EPS() defined as a templated function */
template<typename F>	inline F		SLP_EPS()	{ return (F)0; }
template<>				inline float	SLP_EPS()	{ return 1.0e-7f; }
template<>				inline double	SLP_EPS()	{ return 2.0e-16; }


/* SLP_ABS() defined as a templated function */
template<typename F>	inline F		SLP_ABS(F x)	{ return x > (F)0 ? x : -x; }


/* SLP_NOT_ZERO() defined as a templated function */
template<typename F>	inline bool		SLP_NOT_ZERO(F x)	{ return x >= 2*SLP_EPS<F>() || x <= -2*SLP_EPS<F>(); }


/* Math operators */
#define dot2(a,b)		((a)[0]*(b)[0] + (a)[1]*(b)[1])
#define cross2(a,b) 	((a)[0]*(b)[1] - (a)[1]*(b)[0])


/* status from lp_intersect or linprog */ 
#define SLP_INFEASIBLE 0
#define SLP_MINIMUM 1
#define SLP_UNBOUNDED 2
#define SLP_AMBIGUOUS 3

/* status from plane_down */
#define SLP_REDUNDANT 0
#define SLP_PROPER 1


template<typename F>
inline int
unit2(F a[],F b[],F eps)
{
	F size;
	size = sqrt(a[0]*a[0] + a[1]*a[1]);
	if(size < 2*eps) return(1);
	b[0] = a[0]/size;
	b[1] = a[1]/size;
	return(0);
}


/* returns the index of the plane that is in i's place */
inline int
move_to_front(int i, int next[], int prev[], int m)
{
	int previ;
	if (i == 0 || i == next[0]) return i;
	previ = prev[i];
	/* remove i from its current position */
	next[prev[i]] = next[i];
	prev[next[i]] = prev[i];
	/* put i at the front */
	next[i] = next[0];
	prev[i] = 0;
	prev[next[i]] = i;
	next[0] = i;
	return(previ);
}


template<typename F>
int
wedge(
	F halves[][2],
	int m,
	int next[],
	int prev[],
	F cw_vec[],
	F ccw_vec[],
	int *degen
)
{
	int i;
	F d_cw, d_ccw;
	int offensive;

	*degen = 0;
	for(i=0;i!=m;i = next[i]) {
		if(!unit2(halves[i],ccw_vec,SLP_EPS<F>())) {
/* clock-wise */
			cw_vec[0] = ccw_vec[1];
			cw_vec[1] = -ccw_vec[0];
/* counter-clockwise */
			ccw_vec[0] = -cw_vec[0];
			ccw_vec[1] = -cw_vec[1];
			break;
		}
	}
	if(i==m) return(SLP_UNBOUNDED);
	i = 0;
	while(i!=m) {
		offensive = 0;
		d_cw = dot2(cw_vec,halves[i]);
		d_ccw = dot2(ccw_vec,halves[i]);
		if(d_ccw >= 2*SLP_EPS<F>()) {
			if(d_cw <= -2*SLP_EPS<F>()) {
				cw_vec[0] = halves[i][1];
				cw_vec[1] = -halves[i][0];
				(void)unit2(cw_vec,cw_vec,SLP_EPS<F>());
				offensive = 1;
			}
		} else if(d_cw >= 2*SLP_EPS<F>()) {
			if(d_ccw <= -2*SLP_EPS<F>()) {
				ccw_vec[0] = -halves[i][1];
				ccw_vec[1] = halves[i][0];
				(void)unit2(ccw_vec,ccw_vec,SLP_EPS<F>());
				offensive = 1;
			}
		} else if(d_ccw <= -2*SLP_EPS<F>() && d_cw <= -2*SLP_EPS<F>()) {
			return(SLP_INFEASIBLE);
		} else if((d_cw <= -2*SLP_EPS<F>()) ||
			(d_ccw <= -2*SLP_EPS<F>()) ||
			(cross2(cw_vec,halves[i]) < (F)0)) {
/* degenerate */
			if(d_cw <= -2*SLP_EPS<F>()) {
				(void)unit2(ccw_vec,cw_vec,SLP_EPS<F>());
			} else if(d_ccw <= -2*SLP_EPS<F>()) { 
				(void)unit2(cw_vec,ccw_vec,SLP_EPS<F>());
			}
			*degen = 1;
			offensive = 1;
		}
/* place this offensive plane in second place */
		if(offensive) i = move_to_front(i,next,prev,m);
		i = next[i];
		if(*degen) break;
	}
	if(*degen) {
		while(i!=m) {
			d_cw = dot2(cw_vec,halves[i]);
			d_ccw = dot2(ccw_vec,halves[i]);
			if(d_cw < -2*SLP_EPS<F>()) {
				if(d_ccw < -2*SLP_EPS<F>()) {
					return(SLP_INFEASIBLE);
				} else {
					cw_vec[0] = ccw_vec[0];
					cw_vec[1] = ccw_vec[1];
				}
			} else if(d_ccw < -2*SLP_EPS<F>()) {
				ccw_vec[0] = cw_vec[0];
				ccw_vec[1] = cw_vec[1];
			}
			i = next[i];
		}
	}
	return(SLP_MINIMUM);
}


template<typename F>
void
lp_min_lin_rat
(
	int degen,
	F cw_vec[2],
	F ccw_vec[2],
	F n_vec[2],
	F d_vec[2],
	F opt[2]
)
{
	F d_cw, d_ccw, n_cw, n_ccw;

/* linear rational function case */
	d_cw = dot2(cw_vec,d_vec);
	d_ccw = dot2(ccw_vec,d_vec);
	n_cw = dot2(cw_vec,n_vec);
	n_ccw = dot2(ccw_vec,n_vec);
	if(degen) {
/* if degenerate simply compare values */
		if(n_cw/d_cw < n_ccw/d_ccw) {
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
		} else {
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		}
/* check that the clock-wise and counter clockwise bounds are not near a poles */
	} else if(SLP_NOT_ZERO(d_cw) && SLP_NOT_ZERO(d_ccw)) {
/* the valid region does not contain a poles */
		if(d_cw*d_ccw > (F)0) {
/* find which end has the minimum value */
			if(n_cw/d_cw < n_ccw/d_ccw) {
				opt[0] = cw_vec[0];
				opt[1] = cw_vec[1];
			} else {
				opt[0] = ccw_vec[0];
				opt[1] = ccw_vec[1];
			}
		} else {
/* the valid region does contain a poles */
			if(d_cw > (F)0) {
				opt[0] = -d_vec[1];
				opt[1] = d_vec[0];
			} else {
				opt[0] = d_vec[1];
				opt[1] = -d_vec[0];
			}
		}
	} else if(SLP_NOT_ZERO(d_cw)) {
/* the counter clockwise bound is near a pole */
		if(n_ccw*d_cw > (F)0) {
/* counter clockwise bound is a positive pole */
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
		} else {
/* counter clockwise bound is a negative pole */
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		}
	} else if(SLP_NOT_ZERO(d_ccw)) {
/* the clockwise bound is near a pole */
		if(n_cw*d_ccw > 2*SLP_EPS<F>()) {
/* clockwise bound is at a positive pole */
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		} else {
/* clockwise bound is at a negative pole */
			 opt[0] = cw_vec[0];
			 opt[1] = cw_vec[1];
		} 
	} else {
/* both bounds are near poles */
		if(cross2(d_vec,n_vec) > (F)0) {
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
		} else {
			opt[0] = ccw_vec[0];
			opt[1] = ccw_vec[1];
		}
	}
}


/* unitize a d+1 dimensional point */
template<typename F>
inline int
lp_d_unit(int d, F a[], F b[])
{
	int i;
	F size;

	size = (F)0;
	for (i = 0; i <= d; i++)
		size += a[i] * a[i];
	if (size < (d + 1)*SLP_EPS<F>()*SLP_EPS<F>()) return(1);
	size = (F)1 / sqrt(size);
	for (i = 0; i <= d; i++)
		b[i] = a[i] * size;
	return(0);
}


/* optimize the objective function when there are no contraints */
template<typename F>
int
lp_no_constraints(int d, F n_vec[], F d_vec[], F opt[])
{
	int i;
	F n_dot_d, d_dot_d;

	n_dot_d = (F)0;
	d_dot_d = (F)0;
	for (i = 0; i <= d; i++) {
		n_dot_d += n_vec[i] * d_vec[i];
		d_dot_d += d_vec[i] * d_vec[i];
	}
	if (d_dot_d < SLP_EPS<F>()*SLP_EPS<F>()) {
		d_dot_d = (F)1;
		n_dot_d = (F)0;
	}
	for (i = 0; i <= d; i++) {
		opt[i] = -n_vec[i] + d_vec[i] * n_dot_d / d_dot_d;
	}
	/* normalize the optimal point */
	if (lp_d_unit(d, opt, opt)) {
		opt[d] = (F)1;
		return(SLP_AMBIGUOUS);
	}
	else {
		return(SLP_MINIMUM);
	}
}


/*
 * return the minimum on the projective line
 *
 */
template<typename F>
int
lp_base_case
(
	F halves[][2], 	/* halves --- half lines */
	int m, 				/* m      --- terminal marker */
	F n_vec[2],			/* n_vec  --- numerator funciton */
	F d_vec[2],			/* d_vec  --- denominator function */
	F opt[2],			/* opt    --- optimum  */
	int next[],
	int prev[]			/* next, prev  --- 
					double linked list of indices */
)
{
	F cw_vec[2], ccw_vec[2];
	int degen;
	int status;
	F ab;

/* find the feasible region of the line */
	status = wedge(halves,m,next,prev,cw_vec,ccw_vec,&degen);

	if(status==SLP_INFEASIBLE) return(status);
/* no non-trivial constraints one the plane: return the unconstrained
** optimum */
	if(status==SLP_UNBOUNDED) {
		return(lp_no_constraints(1,n_vec,d_vec,opt));
	}
        ab = SLP_ABS(cross2(n_vec,d_vec));
	if(ab < 2*SLP_EPS<F>()*SLP_EPS<F>()) {
		if(dot2(n_vec,n_vec) < 2*SLP_EPS<F>()*SLP_EPS<F>() ||
			 dot2(d_vec,d_vec) > 2*SLP_EPS<F>()*SLP_EPS<F>()) {
/* numerator is zero or numerator and denominator are linearly dependent */
			opt[0] = cw_vec[0];
			opt[1] = cw_vec[1];
			status = SLP_AMBIGUOUS;
		} else {
/* numerator is non-zero and denominator is zero 
** minimize linear functional on circle */
			if(!degen && cross2(cw_vec,n_vec) <= (F)0 &&
			cross2(n_vec,ccw_vec) <= (F)0 ) {
/* optimum is in interior of feasible region */
				opt[0] = -n_vec[0];
				opt[1] = -n_vec[1];
			} else if(dot2(n_vec,cw_vec) > dot2(n_vec,ccw_vec) ) {
/* optimum is at counter-clockwise boundary */
				opt[0] = ccw_vec[0];
				opt[1] = ccw_vec[1];
			} else {
/* optimum is at clockwise boundary */
				opt[0] = cw_vec[0];
				opt[1] = cw_vec[1];
			}
			status = SLP_MINIMUM;
		}
	} else {
/* niether numerator nor denominator is zero */
		lp_min_lin_rat(degen,cw_vec,ccw_vec,n_vec,d_vec,opt);
		status = SLP_MINIMUM;
	}
#ifdef CHECK
	for(int i=0; i!=m; i=next[i]) {
		F d_cw = dot2(opt,halves[i]);
		if(d_cw < -2*SLP_EPS<F>()) {
			printf("error at base level\n");
			exit(1);
		}
	}
#endif
	return(status);
}


/* find the largest coefficient in a plane */
template<typename F>
void
findimax(F pln[],int idim,int *imax)
{
	F rmax;
	int i;

	*imax = 0;
	rmax = SLP_ABS(pln[0]);
	for(i=1; i<=idim; i++) {
		F ab = SLP_ABS(pln[i]);
		if(ab>rmax) {
			*imax = i;
			rmax = ab;
		}
	}
}


template<typename F>
inline void
vector_up(F equation[],int ivar,int idim,F low_vector[],F vector[])
{
	int i;

	vector[ivar] = (F)0;
	for(i=0; i<ivar; i++) {
		vector[i] = low_vector[i];
		vector[ivar] -= equation[i]*low_vector[i];
	}
	for(i=ivar+1; i<=idim; i++) {
		vector[i] = low_vector[i-1];
		vector[ivar] -= equation[i]*low_vector[i-1];
	}
	vector[ivar] /= equation[ivar];
}


template<typename F>
inline void
vector_down(F elim_eqn[], int ivar, int idim,F old_vec[], F new_vec[])
{
	int i;
	F fac, ve, ee;
	ve = (F)0;
	ee = (F)0;
	for(i=0; i<=idim; i++) {
		ve += old_vec[i]*elim_eqn[i];
		ee += elim_eqn[i]*elim_eqn[i];
	}
	fac = ve/ee;
	for(i=0; i<ivar; i++) {
		new_vec[i] = old_vec[i] - elim_eqn[i]*fac;
	}
	for(i=ivar+1; i<=idim; i++) {
		new_vec[i-1] = old_vec[i] - elim_eqn[i]*fac;
	}
}


template<typename F>
inline void
plane_down(F elim_eqn[], int ivar, int idim,F old_plane[], F new_plane[])
{
	register F crit;
 	register int i;

	crit = old_plane[ivar]/elim_eqn[ivar];
	for(i=0; i<ivar; i++)  {
		new_plane[i] = old_plane[i] - elim_eqn[i]*crit;
	}
	for(i=ivar+1; i<=idim; i++)  {
		new_plane[i-1] = old_plane[i] - elim_eqn[i]*crit;
	}
}


template<typename F>
int
linprog
(
	F halves[], /* halves  --- half spaces */
	int istart,     /* istart  --- should be zero
				 unless doing incremental algorithm */
	int m,  	/* m       --- terminal marker */
	F n_vec[], 	/* n_vec   --- numerator vector */
	F d_vec[], 	/* d_vec   --- denominator vector */
	int d, 		/* d       --- projective dimension */
	F opt[],	/* opt     --- optimum */
	F work[], 	/* work    --- work space (see below) */
	int next[], 	/* next    --- array of indices into halves */
	int prev[], 	/* prev    --- array of indices into halves */
	int max_size 	/* max_size --- size of halves array */
)
/*
**
** half-spaces are in the form
** halves[i][0]*x[0] + halves[i][1]*x[1] + 
** ... + halves[i][d-1]*x[d-1] + halves[i][d]*x[d] >= 0
**
** coefficients should be normalized
** half-spaces should be in random order
** the order of the half spaces is 0, next[0] next[next[0]] ...
** and prev[next[i]] = i
**
** halves[max_size][d+1]
**
** the optimum has been computed for the half spaces
** 0 , next[0], next[next[0]] , ... , prev[istart]
** the next plane that needs to be tested is istart
**
** m is the index of the first plane that is NOT on the list
** i.e. m is the terminal marker for the linked list.
**
** the objective function is dot(x,nvec)/dot(x,dvec)
** if you want the program to solve standard d dimensional linear programming
** problems then n_vec = ( x0, x1, x2, ..., xd-1, 0)
** and           d_vec = (  0,  0,  0, ...,    0, 1)
** and halves[0] = (0, 0, ... , 1)
**
** work points to (max_size+3)*(d+2)*(d-1)/2 F space
*/
{
	int status;
	int i, j, imax;
	F *new_opt, *new_n_vec, *new_d_vec,  *new_halves, *new_work;
	F *plane_i;
	F val;

	if(d==1 && m!=0) {
		return(lp_base_case((F (*)[2])halves,m,n_vec,d_vec,opt,
			next,prev));
	} else {
		int d_vec_zero;
		val = (F)0;
		for(j=0; j<=d; j++) val += d_vec[j]*d_vec[j];
		d_vec_zero = (val < (d+1)*SLP_EPS<F>()*SLP_EPS<F>());

/* find the unconstrained minimum */
		if(!istart) {
			status = lp_no_constraints(d,n_vec,d_vec,opt); 
		} else {
			status = SLP_MINIMUM;
		}
		if(m==0) return(status);
/* allocate memory for next level of recursion */
		new_opt = work;
		new_n_vec = new_opt + d;
		new_d_vec = new_n_vec + d;
		new_halves = new_d_vec + d;
		new_work = new_halves + max_size*d;
		for(i = istart; i!=m; i=next[i]) {
#ifdef SLP_CHECK
			if(i<0 || i>=max_size) {
				printf("index error\n");
				exit(1);
			}
#endif
/* if the optimum is not in half space i then project the problem
** onto that plane */
			plane_i = halves + i*(d+1);
/* determine if the optimum is on the correct side of plane_i */
			val = (F)0;
			for(j=0; j<=d; j++) val += opt[j]*plane_i[j];
			if(val<-(d+1)*SLP_EPS<F>()) {
/* find the largest of the coefficients to eliminate */
			    findimax(plane_i,d,&imax);
/* eliminate that variable */
			    if(i!=0) {
				F fac;
				fac = (F)1/plane_i[imax];
				for(j=0; j!=i; j=next[j]) {
					F *old_plane, *new_plane;
					int k;
					F crit;

					old_plane = halves + j*(d+1);
					new_plane = new_halves + j*d;
					crit = old_plane[imax]*fac;
					for(k=0; k<imax; k++)  {
						new_plane[k] = old_plane[k] - plane_i[k]*crit;
					}
					for(k=imax+1; k<=d; k++)  {
						new_plane[k-1] = old_plane[k] - plane_i[k]*crit;
					}
				}
			    }
/* project the objective function to lower dimension */
			    if(d_vec_zero) {
				vector_down(plane_i,imax,d,n_vec,new_n_vec);
				for(j=0; j<d; j++) new_d_vec[j] = (F)0;
			    } else {
			        plane_down(plane_i,imax,d,n_vec,new_n_vec);
			        plane_down(plane_i,imax,d,d_vec,new_d_vec);
			    }
/* solve sub problem */
			    status = linprog(new_halves,0,i,new_n_vec,
			    new_d_vec,d-1,new_opt,new_work,next,prev,max_size);
/* back substitution */
			    if(status!=SLP_INFEASIBLE) {
				    vector_up(plane_i,imax,d,new_opt,opt);
				{
/* in line code for unit */
				F size;
				size = (F)0;
				for(j=0; j<=d; j++) 
				    size += opt[j]*opt[j];
				size = (F)1/sqrt(size);
				for(j=0; j<=d; j++)
				    opt[j] *= size;
				}
			    } else {
				    return(status);
			    }
/* place this offensive plane in second place */
			    i = move_to_front(i,next,prev,m);
#ifdef SLP_CHECK
			    j=0;
			    while(1) {
/* check the validity of the result */
				val = (F)0;
				for(int k=0; k<=d; k++) 
					val += opt[k]*halves[j*(d+1)+k];
				if(val <-(d+1)*SLP_EPS<F>()) {
				    printf("error\n");
				    exit(1);
				}
				if(j==i)break;
				j=next[j];
			    }
#endif
			} 
		}
 		return(status);
	}
}


template<typename F>
class SeidelLP
{
public:
	// Pass in half-space vectors in a flattened array, with outward-pointing normals.  No normalization is assumed.
	SeidelLP() : m_D(0), m_N(0), m_halves(NULL), m_work(NULL), m_next(NULL), m_prev(NULL), m_n_vec(NULL), m_d_vec(NULL), m_opt(NULL) {}
	SeidelLP(int D, int N) : m_D(0), m_N(0), m_halves(NULL), m_work(NULL), m_next(NULL), m_prev(NULL), m_n_vec(NULL), m_d_vec(NULL), m_opt(NULL) { init(D, N); }

	~SeidelLP() { term(); }

	void
	init(int D, int N)
	{
		if (D < 0) D = 0;
		if (N < 0) N = 0;
		if (D == m_D && N == m_N) return;
		term();
		if (D == 0 || N == 0) return;
		m_D = D;
		m_N = N;
		m_halves = new F[(m_N+1)*(m_D+1)];	// To accomodate halves[0] = (0,0,0,...,1)
		m_work = m_D > 1 ? new F[(m_N+3)*(m_D+2)*(m_D-1)/2] : NULL;
		m_next = new int[m_N];
		m_prev = new int[m_N+1];
		m_n_vec = new F[m_D+1];
		m_d_vec = new F[m_D+1];
		m_opt = new F[m_D+1];
		for (int d = 0; d <= m_D; ++d) m_n_vec[d] = m_d_vec[d] = (F)0;
		m_d_vec[m_D] = (F)1;
	}

	void
	set_half_spaces(const F* half_spaces, int half_space_count, int D)
	{
		init(D, half_space_count);
		if (half_space_count > 0 && half_spaces != NULL)
		{
			m_N = half_space_count;
			set_halves(half_spaces, half_space_count);
		}
	}

	void
	set_objective_fn_coefficients(const F* c)
	{
		for (int d = 0; d < m_D; ++d) m_n_vec[d] = c[d];
	}

	int
	solve(real* solution = NULL)
	{
		const int status = linprog(m_halves, 0, m_N, m_n_vec, m_d_vec, m_D, m_opt, m_work, m_next, m_prev, m_N);
		if (solution != NULL) for (int d = 0; d <= m_D; ++d) solution[d] = m_opt[d];
		return status;
	}
	
	int
	is_feasible()
	{
		switch (solve())
		{
		case SLP_INFEASIBLE:
			return 0;
		case SLP_MINIMUM:
		case SLP_UNBOUNDED:
		case SLP_AMBIGUOUS:
			return 1;
		default:
			return -1;
		}
	}

private:
	void
	set_halves(const F* half_spaces, int half_space_count)
	{
		F* half_space = m_halves;
		for (int i = 0; i < m_D; ++i) *half_space++ = (F)0;
		*half_space++ = (F)1;	// Initial half_space = (0,0,0,...,1)
		for (int i = 0; i < half_space_count; ++i, half_space += m_D+1)
		{
			F norm = (F)0;
			for (int d = 0; d <= m_D; ++d)
			{
				half_space[d] = -(*half_spaces++);	// Negate; this implementation of Seidel's algorithm assumes inward-pointing normals
				norm += half_space[d]*half_space[d];
			}
			// This implementation of Seidel's algorithm assumes normalized half-space vectors, so normalize
			if (norm > (F)0)
			{
				norm = (F)1/sqrt(norm);
				for (int d = 0; d <= m_D; ++d)
				{
					half_space[d] *= norm;
				}
			}
			// Set next, prev arrays as required by implementation
			m_next[i] = i+1;
			m_prev[i] = i-1;
		}
		m_prev[half_space_count] = half_space_count-1;
	}

	void
	term()
	{
		delete [] m_opt;
		m_opt = NULL;
		delete [] m_d_vec;
		m_d_vec = NULL;
		delete [] m_n_vec;
		m_n_vec = NULL;
		delete [] m_prev;
		m_prev = NULL;
		delete [] m_next;
		m_next = NULL;
		delete [] m_work;
		m_work = NULL;
		delete [] m_halves;
		m_halves = NULL;
		m_N = 0;
	}

	int		m_D;
	int		m_N;
	F*		m_halves;
	F*		m_work;
	int*	m_next;
	int*	m_prev;
	F*		m_n_vec;
	F*		m_d_vec;
	F*		m_opt;
};
