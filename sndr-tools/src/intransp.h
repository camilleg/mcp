// intransp.h -- In-place matrix transposition
// language: C++
// author  : Paris Smaragdis

/* Original author:
	http://www.netlib.org/toms/513
	http://calgo.acm.org/
	ACM Collected Algorithms
	file	513
	keywords	transposition in place, matrix transposition, permutation
	gams	D1b3
	title	TRANS
	for	in-situ matrix transposition
	alg	makes use of the cyclic structure of the transposition mapping
	by	E.G. Cate and D.W. Twigg
	ref	ACM TOMS 3 (1977) 104-110
	#	revision of algorithm 380
*/

#ifndef __INTRANSP_H__
#define __INTRANSP_H__

#include <algorithm>

// In-place transposition
template <class T>
void intp( T *a, const int m, const int n)
{
	// Allocate temp memory
	const int iwrk = (m+n)/2;
	int *move = new int[iwrk];
	std::fill_n(move, iwrk, 0);

	// Pretend that arrays are 1-based, not 0-based.
	--a;
	--move;

	int ncount = 2;
	if( m>=3 && n>=3)
	{
		// Calculate the number of fixed points, Euclid's algorithm for GCD(m-1,n-1).
		int ir2 = m - 1;
		int ir1 = n - 1;
		int ir0;
		do {
			ir0 = ir2 % ir1;
			ir2 = ir1;
			ir1 = ir0;
		} while (ir0 != 0);
		ncount += ir2 - 1;
	}

	const int mn = m*n;
	const int k = mn - 1;
	int i1, i2, i1c, i2c;

	// Initial values for search
	int i = 1;
	int im = m;
	// Rearrange at least one loop, by jumping into the middle of the big iteration that follows.
	goto L80;

	// Search for loops to rearrange
	do {
		{
			const int mx = k - i;
			if( ++i > mx)
				return;
			im += m;
			if( im > k)
				im -= k;
			i2 = im;
			if( i == i2)
				continue;
			if( i > iwrk)
				goto L60;
			if( move[i] == 0)
				goto L80;
			continue;

			// Reachable only by "goto L60" from above.
			for (;;)
			{
				i2 = m * i1 - k * (i1 / n);
L60:
				if( i2 <= i || i2 >= mx)
					break;
				i1 = i2;
			}
		}

		if( i2 != i)
			continue;

		// Rearrange the elements of a loop and its companion loop
L80:
		i1 = i;
		const int kmi = k - i;
		i1c = kmi;
		T b = a[i1  + 1];
		T c = a[i1c + 1];
		for (;;) {
			i2 = m * i1 - k * (i1 / n);
			i2c = k - i2;
			if( i1 <= iwrk)
				move[i1] = 2;
			if( i1c <= iwrk)
				move[i1c] = 2;
			ncount += 2;
			if( i2 == i)
				break;
			if( i2 == kmi) {
				std::swap(b, c);
				break;
			}
			a[i1  + 1] = a[i2  + 1];
			a[i1c + 1] = a[i2c + 1];
			i1  = i2;
			i1c = i2c;
		}

		// Store
		a[i1  + 1] = b;
		a[i1c + 1] = c;
	} while (ncount < mn);
	// The original ACM algorithm had more code here to specially handle square matrices.
}

#endif
