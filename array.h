// array.h -- A multidimensional array class
// language: C++
// author  : Paris Smaragdis, 2004
// reviser : Camille Goudeseune, 2013

#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream> // for to_str()
#ifdef __CHECK
#include <limits> // for NaN
#endif

// C++-11 deprecates this with std::to_string().
// Or we could use boost::lexical_cast<std::string>().
// todo: use to_str in other throw()'s.
template <typename T> std::string to_str(const T& t) { std::ostringstream os; os << t; return os.str(); }

// Multidimensional array

template <class T>
class array{
public:
	size_t m, n, k;		// Dimensions.
	T* v;			// Linear array storing the values.
	const bool alias;	// False means *this owns v.  True means caller owns v, and resize() is forbidden.

	// Construct/destruct
	array() : m(0), n(0), k(0), v(NULL), alias(false) {}

	array( const size_t _m, const size_t _n = 1, const size_t _k = 1) 
	 : m(_m), n(_n), k(_k), alias(false)
	{
	  if ((long)m < 0 || (long)n < 0 || (long)k < 0)
	    throw std::runtime_error( "array::array(): negative dimension");
	  v = empty() ? NULL : new T[size()];
	}

	array( T *p, const size_t M, const size_t N = 1, const size_t K = 1) : m(M), n(N), k(K), v(p), alias(true) {}

	array( const array<T>& b) : m(b.m), n(b.n), k(b.k), v(NULL), alias(false)
	{
		try{
			v = new T[size()];
		}
		catch( std::bad_alloc& ){
			throw std::runtime_error( "array::array(): error allocating memory");
		}
		std::copy(b.v, b.v+b.size(), v);
		// memcpy fails when T is a std::string.  So use std::copy even if it might be slower.
	}

	~array() { if( !alias) delete [] v;}

	// Assign
	array<T>& operator=( const array<T>& b)
	{
		resize( b.m, b.n, b.k);
		std::copy(b.v, b.v+b.size(), v);
		return *this;
	}

#ifdef __CHECK
	// foo(1,1,42), while in spirit a vector, may not be single-indexed like a(41).
	// It may only be triple-indexed: a(0,0,41).
	// Conversely, bar(42) *may* be triple-indexed; but only as bar(41,0,0).
#endif

	// Triple indexing
	T& operator()( const size_t i1, const size_t i2, const size_t i3)
	{
#ifdef __CHECK
		if (i1>=m || i2>=n || i3>=k)
			throw std::runtime_error( "array::operator()(i1,i2,i3): index out of bounds");
#endif
		return v[i1+i2*m+i3*(m*n)];
	}

	T operator()( const size_t i1, const size_t i2, const size_t i3) const
	{
#ifdef __CHECK
		if (i1>=m || i2>=n || i3>=k)
			throw std::runtime_error( "array::operator()(i1,i2,i3): index out of bounds");
#endif
		return v[i1+i2*m+i3*(m*n)];
	}

	// Double indexing
	T& operator()( const size_t i1, const size_t i2)
	{
#ifdef __CHECK
		if (k > 1)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + ", " + to_str(i2) + "): missing third index");
		if (i1>=m || i2>=n)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + ", " + to_str(i2) + "): index overflows dimensions " + to_str(m) + ", " + to_str(n));
#endif
		return v[i1+i2*m];
	}

	T operator()( const size_t i1, const size_t i2) const
	{
#ifdef __CHECK
		if (k > 1)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + ", " + to_str(i2) + "): missing third index");
		if (i1>=m || i2>=n)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + ", " + to_str(i2) + "): index overflows dimensions " + to_str(m) + ", " + to_str(n));
#endif
		return v[i1+i2*m];
	}

	// Single indexing
	T& operator()( const size_t i1)
	{
#ifdef __CHECK
		if (empty())
			throw std::runtime_error( "array::operator()(" + to_str(i1) + "): empty array");
		if (n > 1 || k > 1)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + "): missing second or third index");
		if (i1>=m)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + ") overflows dimension " + to_str(m));
#endif
		return v[i1];
	}

	T operator()( const size_t i1) const
	{
#ifdef __CHECK
		if (n > 1 || k > 1)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + "): missing second or third index");
		if (i1>=m)
			throw std::runtime_error( "array::operator()(" + to_str(i1) + ") overflows dimension " + to_str(m));
#endif
		return v[i1];
	}

	// Indexing linearly, ignoring dimensionality.
	T& operator[]( const size_t i1)
	{
#ifdef __CHECK
		if (empty())
			throw std::runtime_error( "array::operator[" + to_str(i1) + "]: empty array");
		if( i1 >= size())
			throw std::runtime_error( "array::operator[" + to_str(i1) + "]: index overflows size " + to_str(size()));
#endif
		return v[i1];
	}

	T operator[]( const size_t i1) const
	{
#ifdef __CHECK
		if (empty())
			throw std::runtime_error( "array::operator[" + to_str(i1) + "]: empty array");
		if( i1 >= size())
			throw std::runtime_error( "array::operator[" + to_str(i1) + "]: index overflows size " + to_str(size()));
#endif
		return v[i1];
	}

	// Convert to a pointer, e.g. for single-value indexing the flattened elements
	// to e.g. multiply two matrices elementwise, or pass an "iterator" to std::copy().
	// (Don't extend this to const pointers, because ambiguity results.)
	operator T*() { return v; }

	// Resize
	void resize(const size_t i1, const size_t i2 = 1, const size_t i3 = 1)
	{
		if( m != i1 || n != i2 || k != i3){
			if( alias)
				throw std::runtime_error( "array::resize(): forbidden for aliased array");
			m = i1; n = i2; k = i3;
			delete [] v; // ok if NULL
			try{
				v = new T[size()];
			}
			catch( std::bad_alloc& ){
				throw std::runtime_error( "array::resize() failed to allocate " + to_str(m) + " x " + to_str(n) + " x " + to_str(k) + " elements\n" );
			}
		}
#if 0
#ifdef __CHECK
		for (size_t i=0; i<size(); ++i)
			v[i] = static_cast<T>(std::numeric_limits<T>::signaling_NaN());
			// This fails to compile if float isn't compatible with T,
			// e.g. T is std::string or has a pile of member variables.
			// ( __builtin_types_compatible_p() works only with C, not C++.)
#endif
#endif
	}

	// Todo: remove special treatment of vectors (1,n,1) and (1,1,k), leaving only (m,1,1).
	// Much cleaner.  And operator<< always did this.

	// If this vector is short, make it longer.
	void grow(const size_t i)
	{
		if (!is_vector())
			throw std::runtime_error( "array::grow(): not a vector");			
		if (i < size())
			return;
		try{
			T* v2 = new T[i];
			std::copy(v, v+size(), v2);
#ifdef __CHECK
			for (size_t j=size(); j<i; ++j)
				v2[j] = std::numeric_limits<T>::signaling_NaN();
#endif
			delete [] v;
			v = v2;
			if (empty()) m=i, n=k=1;
			else if(m == size()) m = i;
			else if(n == size()) n = i;
			else if(k == size()) k = i;
		}
		catch( std::bad_alloc& ){
			throw std::runtime_error( "array::grow(): error allocating memory");
		}
		assert(size() >= i);
	}

	// Append to a vector.  Expensive: linear in the vector's length!
	void push_back(const T a)
	{
		if (!is_vector())
			throw std::runtime_error( "array::push_back(): not a vector");			
		try{
			T *v2 = new T[size()+1];
			std::copy(v, v+size(), v2);
			delete [] v;
			v = v2;
			if (empty()) m=n=k=1;
			else if( m == size()) ++m;
			else if( n == size()) ++n;
			else if( k == size()) ++k;
		}
		catch( std::bad_alloc& ){
			throw std::runtime_error( "array::push_back(): error allocating memory");
		}
		v[size()-1] = a;
	}
	const T& back() const
	{
		if (!is_vector())
			throw std::runtime_error( "array::back(): not a vector");			
		return v[size()-1];
	}

	// Scale elements so they sum to unity.
	void normalize()
	{
	  const T a = sum();
	  if (a == 0)
	    return;
	  for (size_t i=0; i<size(); ++i) v[i] /= a;
	}

	// Scale log-elements so they sum to unity.
	void normalize_log()
	{
	  T sum = 0;
	  for (size_t i=0; i<size(); ++i) sum += exp(v[i]);
	  if (sum == 0)
	    return;
	  for (size_t i=0; i<size(); ++i) v[i] = log(exp(v[i]) / sum);
	}

	// Normalize log-elements in another way.
	void normalize_log_also()
	{
	  const T a = log(sum());
	  for (size_t i=0; i<size(); ++i) v[i] -= a;
	}

	// Scale elements so largest element is unity.
	// Assumes that all elements are nonnegative.
	void normalize_max1()
	{
	  T a = 0;
	  for (size_t i=0; i<size(); ++i) a = std::max(a, v[i]);
	  if (a == 0)
	    return;
	  for (size_t i=0; i<size(); ++i) v[i] /= a;
	}

	// Scale elements in j'th row so they sum to unity.
	void normalize(const size_t j)
	{
	  const T t = sum(j);
	  if (t == 0)
	    return;
	  for (size_t i=0; i<m; ++i) v[i+j*m] /= t;
	}

	// Sum of elements.
	T sum() const
	{
	  T a = 0;
	  for (size_t i=0; i<size(); ++i) a += v[i];
	  return a;
	}

	// Sum of elements in j'th row.
	T sum(const size_t j) const
	{
	  if (j >= n)
	    throw std::runtime_error( "array::sum(" + to_str(j) + ") arg exceeds second dimension " + to_str(n));			
	  T a = 0;
	  for (size_t i=0; i<m; ++i) a += v[i+j*m]; // v(i,j);
	  return a;
	}

	// Number of elements
	size_t size() const { return m*n*k; }
	bool empty() const { return size() == 0; }
	bool is_vector() const { return m == size() || n == size() || k == size() || empty(); }
	void clear() { resize(0,0,0); }

	// Return min and max element
	T min() const { return *std::min_element(v, v+size()); }
	T max() const { return *std::max_element(v, v+size()); }

	// todo: normalize values of m,n,k so if one of them is zero, all of them are;
	// i.e., forbid multiple representations with the same meaning (empty()).
	
	void write(std::ofstream& f) const { f.write((const char*)v,  size()*sizeof(T)); }

	void read(std::ifstream& f, const size_t i1, const size_t i2=1, const size_t i3=1)
	{
	  resize(i1, i2, i3);
	  f.read((char*)v, size()*sizeof(T));
	}
};

// ostream interface
template <class T>
std::ostream& operator<<( std::ostream& o, const array<T>& a)
{
	// Vector
	if( a.n == 1 && a.k == 1){
		for( size_t i = 0 ; i < a.m ; i++)
			o << a(i) << ' ';
		o << std::endl;

	// Matrix
	}else if( a.k == 1){
		for( size_t i = 0 ; i < a.m ; i++){
			for( size_t j = 0 ; j < a.n ; j++)
				o << a(i,j) << ' ';
			o << std::endl;
		}

	// And beyond ...
	}else{
		for( size_t i = 0 ; i < a.m ; i++){
			for( size_t j = 0 ; j < a.n ; j++)
				for( size_t k = 0 ; k < a.k ; k++){
					o << a(i,j,k) << ' ';
				o << std::endl;
			}
		}
	}
	return o;
}

#endif
