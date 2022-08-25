/****************************************************************************
* Indirect predicates for geometric constructions					        *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2019: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU Lesser General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or (at  *
* your option) any later version.                                           *
*                                                                           *
* This program is distributed in the hope that it will be useful, but       *
* WITHOUT ANY WARRANTY; without even the implied warranty of                *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  *
* General Public License for more details.                                  *
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  *
* along with this program.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/

#ifndef NUMERICS_H
#define NUMERICS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <fenv.h>
#include <iostream>

	inline void ip_error(const char* msg)
	{
		fprintf(stderr, msg);
		exit(0);
	}

#if INTPTR_MAX == INT64_MAX
#define	IS64BITPLATFORM
#endif

#ifdef _MSC_VER
#define	ISVISUALSTUDIO
#endif

	// 64-bit
#ifdef IS64BITPLATFORM
#define USE_SIMD_INSTRUCTIONS
#endif

#ifdef ISVISUALSTUDIO

#pragma fenv_access (on)

	inline void setFPUModeToRoundUP() { _controlfp(_RC_UP, _MCW_RC); }
	inline void setFPUModeToRoundNEAR() { _controlfp(_RC_NEAR, _MCW_RC); }

#else

#pragma STDC FENV_ACCESS ON

	inline void setFPUModeToRoundUP() { fesetround(FE_UPWARD); }
	inline void setFPUModeToRoundNEAR() { fesetround(FE_TONEAREST); }
#endif

#ifdef USE_SIMD_INSTRUCTIONS
#include <emmintrin.h>

#pragma intrinsic(fma)

	//inline void setFPUModeToRoundUP() { _MM_SET_ROUNDING_MODE(_MM_ROUND_UP); }
	//inline void setFPUModeToRoundNEAR() { _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST); }
	

	class interval_number
	{
		__m128d interval; // interval[1] = min_low, interval[0] = high

		static __m128d zero;
		static __m128i sign_low_mask, sign_high_mask;

	public:
		static void init();

		const double *getInterval() const { return (const double*)&interval; }

		inline interval_number() { }
		inline interval_number(const double a) : interval(_mm_set_pd(-a, a)) {}
		inline interval_number(const double minf, const double sup) : interval(_mm_set_pd(minf, sup)) {}
		inline interval_number(const __m128d& i) : interval(i) {}
		inline interval_number(const interval_number& b) : interval(b.interval) {}

		inline double inf() const { return -((double*)(&interval))[1]; }
		inline double sup() const { return ((double*)(&interval))[0]; }
		inline double width() const { return sup() - inf(); }

		inline void invert() { interval = _mm_shuffle_pd(interval, interval, 1); }

		inline bool isNegative() const { return _mm_comilt_sd(interval, zero); }
		inline bool isPositive() const { return _mm_comilt_sd(_mm_shuffle_pd(interval, interval, 1), zero); }

		inline bool signIsReliable() const { return (isNegative() || isPositive()); } // Zero is not accounted for
		inline int sign() const { return (isNegative()) ? (-1) : (1); } // Zero is not accounted for

		inline bool isNAN() const {	return sup() != sup(); }

		inline bool operator<(const double b) const { return (_mm_comilt_sd(interval, _mm_set_sd(b))); }

		inline interval_number& operator=(const interval_number& b) { interval = b.interval; return *this; }

		inline interval_number operator+(const interval_number& b) const { return interval_number(_mm_add_pd(interval, b.interval)); }

		inline interval_number operator-(const interval_number& b) const { return interval_number(_mm_add_pd(interval, _mm_shuffle_pd(b.interval, b.interval, 1))); }

		inline interval_number operator*(const interval_number& b) const
		{
			__m128i ssg;
			__m128d llhh, lhhl, ip;

			switch ((_mm_movemask_pd(interval) << 2) | _mm_movemask_pd(b.interval))
			{
			case 0:
				llhh = _mm_mul_pd(interval, b.interval);
				lhhl = _mm_mul_pd(interval, _mm_shuffle_pd(b.interval, b.interval, 1));
				return interval_number(_mm_max_pd(_mm_unpacklo_pd(llhh, lhhl), _mm_unpackhi_pd(llhh, lhhl)));
			case 1:
				return interval_number(_mm_mul_pd(_mm_shuffle_pd(b.interval, b.interval, 3), _mm_shuffle_pd(interval, interval, 1)));
			case 2:
				return interval_number(_mm_mul_pd(_mm_shuffle_pd(b.interval, b.interval, 0), interval));
			case 4:
				return interval_number(_mm_mul_pd(_mm_shuffle_pd(interval, interval, 3), _mm_shuffle_pd(b.interval, b.interval, 1)));
			case 5:
				ip = _mm_mul_pd(_mm_castsi128_pd(_mm_xor_si128(_mm_castpd_si128(interval), sign_high_mask)), b.interval);
				return interval_number(_mm_shuffle_pd(ip, ip, 1));
			case 6:
				ssg = _mm_xor_si128(_mm_castpd_si128(b.interval), sign_low_mask);
				return interval_number(_mm_mul_pd(interval, _mm_shuffle_pd(_mm_castsi128_pd(ssg), _mm_castsi128_pd(ssg), 1)));
			case 8:
				return interval_number(_mm_mul_pd(_mm_shuffle_pd(interval, interval, 0), b.interval));
			case 9:
				ssg = _mm_xor_si128(_mm_castpd_si128(interval), sign_low_mask);
				return interval_number(_mm_mul_pd(b.interval, _mm_shuffle_pd(_mm_castsi128_pd(ssg), _mm_castsi128_pd(ssg), 1)));
			case 10:
				return interval_number(_mm_mul_pd(interval, _mm_castsi128_pd(_mm_xor_si128(_mm_castpd_si128(b.interval), sign_low_mask))));
			}

			return interval_number(NAN);
		}

	};

#else // USE_SIMD_INSTRUCTIONS

	// Interval_number
	class interval_number
	{
		typedef union error_approx_type_t
		{
			double d;
			uint64_t u;

			inline error_approx_type_t() {}
			inline error_approx_type_t(double a) : d(a) {}
			inline uint64_t is_negative() const { return u >> 63; }
		} casted_double;

	public:
		double min_low, high;

		static void init() {}

		const double* getInterval() const { return (const double*)&min_low; }

		inline interval_number() { }
		inline interval_number(double a) : min_low(-a), high(a) {}
		inline interval_number(double minf, double sup) : min_low(minf), high(sup) {}
		inline interval_number(const interval_number& b) : min_low(b.min_low), high(b.high) {}

		inline double inf() const { return -min_low; }
		inline double sup() const { return high; }
		inline double width() const { return sup() - inf(); }

		inline bool isNegative() const { return (high < 0); }
		inline void invert() { double tmp = min_low; min_low = high; high = tmp; }

		inline bool signIsReliable() const { return (min_low < 0 || high < 0); }
		inline int sign() const { return (min_low < 0) ? (1) : ((high < 0) ? (-1) : (0)); }

		inline bool isNAN() const { return (high != high); }

		inline double center() const { return 0.5 * (high - min_low); }

		inline bool operator<(const double b) const { return (high < b); }

		inline interval_number& operator=(const interval_number& b) { min_low = b.min_low; high = b.high; return *this; }

		inline interval_number operator+(const interval_number& b) const { return interval_number(min_low + b.min_low, high + b.high); }

		inline interval_number operator-(const interval_number& b) const { return interval_number(b.high + min_low, high + b.min_low); }

		inline interval_number operator*(const interval_number& b) const
		{
			casted_double l1(min_low), h1(high), l2(b.min_low), h2(b.high);
			uint64_t cfg = (l1.is_negative() << 3) + (h1.is_negative() << 2) + (l2.is_negative() << 1) + (h2.is_negative());

			switch (cfg)
			{
			case 10: return interval_number(min_low * (-b.min_low), high * b.high);
			case 8: return interval_number(high * b.min_low, high * b.high);
			case 9: return interval_number(high * b.min_low, (-min_low) * b.high);
			case 2: return interval_number(min_low * b.high, high * b.high);
			case 0:
				double ll, lh, hl, hh;
				ll = min_low * b.min_low; lh = (min_low * b.high); hl = (high * b.min_low); hh = high * b.high;
				if (hl > lh) lh = hl;
				if (ll > hh) hh = ll;
				return interval_number(lh, hh);
			case 1: return interval_number(high * b.min_low, min_low * b.min_low);
			case 6: return interval_number(min_low * b.high, high * (-b.min_low));
			case 4: return interval_number(min_low * b.high, min_low * b.min_low);
			case 5: return interval_number((-high) * b.high, min_low * b.min_low);
			};

			return interval_number(NAN);
		}
	};
#endif // USE_SIMD_INSTRUCTIONS

	inline std::ostream& operator<<(std::ostream& os, const interval_number& p)
	{
		os << "[ " << p.inf() << ", " << p.sup() << " ]";
		return os;
	}


	// The following macros are fast implementations of basic expansion arithmetic due
	// to Dekker, Knuth, Priest, Shewchuk, and others.

	// See Y. Hida, X. S. Li,  D. H. Bailey "Algorithms for Quad-Double Precision Floating Point Arithmetic"

	// Sums
#define Quick_Two_Sum(a, b, x, y) x = a + b; y = b - (x - a)
#define Two_Sum(a, b, x, y) x = a + b; _bv = x - a; y = (a - (x - _bv)) + (b - _bv)
#define Two_One_Sum(a1, a0, b, x2, x1, x0) Two_Sum(a0, b , _i, x0); Two_Sum(a1, _i, x2, x1)

// Differences
#define Two_Diff(a, b, x, y) x = a - b; _bv = a - x; y = (a - (x + _bv)) + (_bv - b)
#define Two_One_Diff(a1, a0, b, x2, x1, x0) Two_Diff(a0, b , _i, x0); Two_Sum( a1, _i, x2, x1)

// Products
#define Split(a, _ah, _al) _c = 1.3421772800000003e+008 * a; _ah = _c - (_c - a); _al = a - _ah
#define Two_Prod_PreSplit(a, b, _bh, _bl, x, y) x = a * b; Split(a, _ah, _al); y = (_al * _bl) - (((x - (_ah * _bh)) - (_al * _bh)) - (_ah * _bl))
#define Two_Product_2Presplit(a, _ah, _al, b, _bh, _bl, x, y) x = a * b; y = (_al * _bl) - (((x - _ah * _bh) - (_al * _bh)) - (_ah * _bl))


// An instance of the following must be created to access functions for expansion arithmetic
	class expansionObject
	{
		// Temporary vars used in low-level arithmetic
		double _bv, _c, _ah, _al, _bh, _bl, _i, _j, _k, _l, _0, _1, _2, _u3;

	public:
		expansionObject() {}

		inline void two_Sum(const double a, const double b, double* xy) { Two_Sum(a, b, xy[1], xy[0]); }

		inline void two_Diff(const double a, const double b, double* xy) { Two_Diff(a, b, xy[1], xy[0]); }

		// [x,y] = [a]*[b]		 Multiplies two expansions [a] and [b] of length one
		inline void Two_Prod(const double a, const double b, double& x, double& y)
		{
			x = a * b;
//			y = fma(a, b, -x);
			Split(a, _ah, _al); Split(b, _bh, _bl);
			y = ((_ah * _bh - x) + _ah * _bl + _al * _bh) + _al * _bl;
		}
		inline void Two_Prod(const double a, const double b, double* xy) { Two_Prod(a, b, xy[1], xy[0]); }


		// [x,y] = [a]^2		Squares an expansion of length one
		inline void Square(const double a, double& x, double& y)
		{
			x = a * a;
			Split(a, _ah, _al);
			y = (_al * _al) - ((x - (_ah * _ah)) - ((_ah + _ah) * _al));
		}
		inline void Square(const double a, double* xy) { Square(a, xy[1], xy[0]); }

		// [x2,x1,x0] = [a1,a0]-[b]		Subtracts an expansion [b] of length one from an expansion [a1,a0] of length two
		inline void two_One_Diff(const double a1, const double a0, const double b, double& x2, double& x1, double& x0)
		{ Two_One_Diff(a1, a0, b, x2, x1, x0); }
		inline void two_One_Diff(const double* a, const double b, double* x) { two_One_Diff(a[1], a[0], b, x[2], x[1], x[0]); }

		// [x3,x2,x1,x0] = [a1,a0]*[b]		Multiplies an expansion [a1,a0] of length two by an expansion [b] of length one
		inline void Two_One_Prod(const double a1, const double a0, const double b, double& x3, double& x2, double& x1, double& x0)
		{
			Split(b, _bh, _bl);
			Two_Prod_PreSplit(a0, b, _bh, _bl, _i, x0); Two_Prod_PreSplit(a1, b, _bh, _bl, _j, _0);
			Two_Sum(_i, _0, _k, x1); Quick_Two_Sum(_j, _k, x3, x2);
		}
		inline void Two_One_Prod(const double* a, const double b, double* x) { Two_One_Prod(a[1], a[0], b, x[3], x[2], x[1], x[0]); }

		// [x3,x2,x1,x0] = [a1,a0]+[b1,b0]		Calculates the sum of two expansions of length two
		inline void Two_Two_Sum(const double a1, const double a0, const double b1, const double b0, double& x3, double& x2, double& x1, double& x0)
		{
			Two_One_Sum(a1, a0, b0, _j, _0, x0); Two_One_Sum(_j, _0, b1, x3, x2, x1);
		}
		inline void Two_Two_Sum(const double* a, const double* b, double* xy) { Two_Two_Sum(a[1], a[0], b[1], b[0], xy[3], xy[2], xy[1], xy[0]); }

		// [x3,x2,x1,x0] = [a1,a0]-[b1,b0]		Calculates the difference between two expansions of length two
		inline void Two_Two_Diff(const double a1, const double a0, const double b1, const double b0, double& x3, double& x2, double& x1, double& x0)
		{
			Two_One_Diff(a1, a0, b0, _j, _0, x0); Two_One_Diff(_j, _0, b1, _u3, x2, x1); x3 = _u3;
		}
		inline void Two_Two_Diff(const double* a, const double* b, double* x) { Two_Two_Diff(a[1], a[0], b[1], b[0], x[3], x[2], x[1], x[0]); }

		// Calculates the second component 'y' of the expansion [x,y] = [a]-[b] when 'x' is known
		inline void Two_Diff_Back(const double a, const double b, double& x, double& y) { _bv = a - x; y = (a - (x + _bv)) + (_bv - b); }
		inline void Two_Diff_Back(const double a, const double b, double* xy) { Two_Diff_Back(a, b, xy[1], xy[0]); }

		// [h] = [a1,a0]^2		Squares an expansion of length 2
		// 'h' must be allocated by the caller with 6 components.
		void Two_Square(const double& a1, const double& a0, double* x);

		// [h7,h6,...,h0] = [a1,a0]*[b1,b0]		Calculates the product of two expansions of length two.
		// 'h' must be allocated by the caller with eight components.
		void Two_Two_Prod(const double a1, const double a0, const double b1, const double b0, double* h);
		inline void Two_Two_Prod(const double* a, const double* b, double* xy) { Two_Two_Prod(a[1], a[0], b[1], b[0], xy); }

		// [h7,h6,...,h0] = [a1,a0]*[b1,b0]		Calculates the product of two expansions of length two.
		// 'h' must be allocated by the caller with eight components.
		//void Two_Two_Prod(const double a1, const double a0, const double b1, const double b0, double *h);
		//inline void Two_Two_Prod(const double *a, const double *b, double *xy) { Two_Two_Prod(a[1], a[0], b[1], b[0], xy); }

		// [e] = -[e]		Inplace inversion
		static void Gen_Invert(const int elen, double* e) { for (int i = 0; i < elen; i++) e[i] = -e[i]; }

		// [h] = [e] + [f]		Sums two expansions and returns number of components of result
		// 'h' must be allocated by the caller with at least elen+flen components.
		int Gen_Sum(const int elen, const double* e, const int flen, const double* f, double* h);

		// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
		inline int Gen_Sum_With_Alloc(const int elen, const double* e, const int flen, const double* f, double** h)
		{
			*h = (double*)malloc((elen + flen) * sizeof(double));
			return Gen_Sum(elen, e, flen, f, *h);
		}

		// [h] = [e] + [f]		Subtracts two expansions and returns number of components of result
		// 'h' must be allocated by the caller with at least elen+flen components.
		int Gen_Diff(const int elen, const double* e, const int flen, const double* f, double* h);

		// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
		inline int Gen_Diff_With_Alloc(const int elen, const double* e, const int flen, const double* f, double** h)
		{
			*h = (double*)malloc((elen + flen) * sizeof(double));
			return Gen_Diff(elen, e, flen, f, *h);
		}

		// [h] = [e] * b		Multiplies an expansion by a scalar
		// 'h' must be allocated by the caller with at least elen*2 components.
		int Gen_Scale(const int elen, const double* e, const double& b, double* h);

		// [h] = [e] * 2		Multiplies an expansion by 2
		// 'h' must be allocated by the caller with at least elen components. This is exact up to overflows.
		inline void Double(const int elen, const double* e, double* h) const { for (int i = 0; i < elen; i++) h[i] = 2 * e[i]; }
		
		// [h] = [e] * 2		Multiplies an expansion by n
		// If 'n' is a power of two, the multiplication is exact
		inline static void ExactScale(const int elen, double* e, const double n) { for (int i = 0; i < elen; i++) e[i] *= n; }

		// [h] = [a] * [b]
		// 'h' must be allocated by the caller with at least 2*alen*blen components.
		int Sub_product(const int alen, const double* a, const int blen, const double* b, double* h);

		// [h] = [a] * [b]
		// 'h' must be allocated by the caller with at least MAX(2*alen*blen, 8) components.
		int Gen_Product(const int alen, const double* a, const int blen, const double* b, double* h);

		// Same as above, but 'h' is allocated internally. The caller must still call 'free' to release the memory.
		inline int Gen_Product_With_Alloc(const int alen, const double* a, const int blen, const double* b, double** h)
		{
			int h_len = alen * blen * 2;
			if (h_len < 8) h_len = 8;
			*h = (double*)malloc(h_len * sizeof(double));
			return Gen_Product(alen, a, blen, b, *h);
		}


		// Assume that *h is pre-allocated with hlen doubles.
		// If more elements are required, *h is re-allocated internally.
		// In any case, the function returns the size of the resulting expansion.
		// The caller must verify whether reallocation took place, and possibly call 'free' to release the memory.
		// When reallocation takes place, *h becomes different from its original value.

		inline int Double_With_PreAlloc(const int elen, const double* e, double** h, const int hlen)
		{
			int newlen = elen;
			if (hlen < newlen) *h = (double*)malloc(newlen * sizeof(double));
			//if (hlen < newlen) printf("REALLOC %d bytes\n", newlen);
			Double(elen, e, *h);
			return newlen;
		}

		inline int Gen_Scale_With_PreAlloc(const int elen, const double* e, const double& b, double** h, const int hlen)
		{
			int newlen = elen * 2;
			if (hlen < newlen) *h = (double*)malloc(newlen * sizeof(double));
			return Gen_Scale(elen, e, b, *h);
		}

		inline int Gen_Sum_With_PreAlloc(const int elen, const double* e, const int flen, const double* f, double** h, const int hlen)
		{
			int newlen = elen + flen;
			if (hlen < newlen) *h = (double*)malloc(newlen * sizeof(double));
			return Gen_Sum(elen, e, flen, f, *h);
		}

		inline int Gen_Diff_With_PreAlloc(const int elen, const double* e, const int flen, const double* f, double** h, const int hlen)
		{
			int newlen = elen + flen;
			if (hlen < newlen) *h = (double*)malloc(newlen * sizeof(double));
			return Gen_Diff(elen, e, flen, f, *h);
		}

		inline int Gen_Product_With_PreAlloc(const int alen, const double* a, const int blen, const double* b, double** h, const int hlen)
		{
			int newlen = alen * blen * 2;
			if (hlen < newlen || hlen < 8)
			{
				if (newlen < 8) newlen = 8;
				*h = (double*)malloc(newlen * sizeof(double));
			}
			return Gen_Product(alen, a, blen, b, *h);
		}

		// Approximates the expansion to a double
		static double To_Double(const int elen, const double* e);

		static void print(const int elen, const double* e) { for (int i = 0; i < elen; i++) printf("%e ", e[i]); printf("\n");}
	};


	/////////////////////////////////////////////////////////////////////
	// 	   
	// 	   B I G   N A T U R A L
	// 
	/////////////////////////////////////////////////////////////////////

	// A bignatural is an arbitrarily large non-negative integer.
	// It is made of a sequence of digits in base 2^32.
	// Leading zero-digits are not allowed.
	// The value 'zero' is represented by an empty digit sequence.

	const uint32_t static_digit_count = 8;

	// Uses the stack for up to 'static_digit_count' elements, so that
	// 'small' numbers can be processed efficiently.
	// When a number grows above the threshold digits, the heap is used.
	class bignatural {
		uint32_t m_static[static_digit_count];	// Local vector of digits
		uint32_t* digits;	// Ptr to the digits. This may point to m_static
		uint32_t m_size;		// Actual number of digits
		uint32_t m_capacity;	// Current vector capacity

		inline void checkAndRelease() { if (digits != m_static) free(digits); }

		inline void init(const bignatural& m) {
			m_size = m.m_size;
			m_capacity = m.m_capacity;
			if (m.digits != m.m_static) digits = m.dupDigits();
			else {
				for (uint32_t i = 0; i < static_digit_count; i++) m_static[i] = m.m_static[i];
				digits = m_static;
			}
		}

		inline void init(const uint64_t m) {
			m_capacity = static_digit_count;

			if (m == 0) {
				m_size = 0;
			}
			else if (m <= UINT32_MAX) {
				m_size = 1;
				m_static[0] = (uint32_t)m;
			}
			else {
				m_size = 2;
				m_static[0] = (uint32_t)(m >> 32);
				m_static[1] = (uint32_t)(m & UINT32_MAX);
			}
			digits = m_static;
		}

	public:
		// Creates a 'zero'
		inline bignatural() : m_size(0), m_capacity(static_digit_count) { 
			digits = m_static;
		}

		inline ~bignatural() { checkAndRelease(); }

		inline bignatural(const bignatural& m) { init(m); }

		// Creates from uint64_t
		inline bignatural(uint64_t m) { init(m); }

		// If the number fits a uint64_t convert and return true
		inline bool toUint64(uint64_t& n) const {
			if (m_size == 0) n = 0;
			else if (m_size == 1) n = digits[0];
			else if (m_size == 2) { n = (((uint64_t)digits[0]) << 32) + digits[1]; }
			else return false;

			return true;
		}

		// If the number fits a uint32_t convert and return true
		inline bool toUint32(uint32_t& n) const {
			if (m_size == 0) n = 0;
			else if (m_size == 1) n = digits[0];
			else return false;

			return true;
		}

		inline bignatural& operator=(const bignatural& m) {
			checkAndRelease();
			init(m);
			return *this;
		}

		inline bignatural& operator=(const uint64_t m) { 
			checkAndRelease();
			init(m);
			return *this;
		}

		inline const uint32_t& back() const { return digits[m_size - 1]; }

		inline const uint32_t& operator[](int i) const { return digits[i]; }

		inline uint32_t size() const { return m_size; }

		inline bool empty() const { return m_size == 0; }

		// Left-shift by n bits and possibly add limbs as necessary
		void operator<<=(uint32_t n) {
			uint32_t s = n & 0x0000001f;
			uint32_t lz = countLeadingZeroes();
			if (lz < s) { // Need a further limb
				push_back(0);
				s = 32 - s;
				for (int i = (int)m_size - 1; i > 0; i--) {
					digits[i] >>= s;
					digits[i] |= (digits[i - 1] << (32 - s));
				}
				digits[0] >>= s;
			}
			else if (s) { // Leading zeroes are enough
				for (int i = 0; i < (int)m_size - 1; i++) {
					digits[i] <<= s;
					digits[i] |= (digits[i + 1] >> (32 - s));
				}
				back() <<= s;
			}

			while (n >= 32) {
				push_back(0);
				n -= 32;
			}
		}

		// Right-shift by n bits
		void operator>>=(uint32_t n) {
			while (n >= 32) {
				pop_back();
				n -= 32;
			}
			for (uint32_t i = m_size; i > 1; i--) {
				digits[i - 1] >>= n;
				digits[i - 1] |= (digits[i - 2] << (32 - n));
			}
			digits[0] >>= n;
			if (digits[0] == 0) pop_front();
		}

		inline bool operator>=(const bignatural& b) const {
			const int s = (size() > b.size()) - (size() < b.size());
			if (s) return (s > 0);
			uint32_t i;
			for (i = 0; i < size() && digits[i] == b.digits[i]; i++);
			return (i==size() || digits[i] > b.digits[i]);
		}

		inline bool operator>(const bignatural& b) const {
			const int s = (size() > b.size()) - (size() < b.size());
			if (s) return (s > 0);
			uint32_t i;
			for (i = 0; i < size() && digits[i] == b.digits[i]; i++);
			return (i != size() && digits[i] > b.digits[i]);
		}

		inline bignatural operator+(const bignatural& b) const {
			bignatural result;
			result.toSum(*this, b);
			return result;
		}

		// Assume that b is smaller than or equal to this number!
		inline bignatural operator-(const bignatural& b) const {
			bignatural result;
			result.toDiff(*this, b);
			return result;
		}

		inline bignatural operator*(const bignatural& b) const {
			bignatural result;
			result.toProd(*this, b);
			return result;
		}

		// Short division algorithm
		inline bignatural divide_by(const uint32_t D, uint32_t& remainder) const {
			if (D == 0) ip_error("Division by zero\n");
			if (m_size == 0) return 0;

			// If both dividend fits into 64 bits, use hardware division
			uint64_t n;
			if (toUint64(n)) {
				remainder = n % D;
				return n / D;
			}

			bignatural Q;
			uint32_t next_digit = 0;
			uint64_t dividend = digits[next_digit++];
			for (;;) {
				uint64_t tmp_div = dividend / D;
				if (!Q.empty() || tmp_div) Q.push_back((uint32_t)tmp_div);
				dividend -= (tmp_div * D);
				if (next_digit < m_size) {
					dividend <<= 32;
					dividend += digits[next_digit++];
				}
				else break;
			}
			remainder = (uint32_t)dividend;

			return Q;
		}

		uint32_t getNumSignificantBits() const {
			if (!m_size) return 0;
			int nsb = 31;
			while (!(digits[0] & (1 << nsb))) nsb--;
			nsb++;
			return nsb + (m_size - 1) * 32;
		}

		inline bool getBit(uint32_t b) const {
			const uint32_t dig = (m_size - (b >> 5)) - 1;
			const uint32_t bit = b & 31;
			return (digits[dig] & (1 << bit));
		}

		// Long division
		inline bignatural divide_by(const bignatural& divisor, bignatural& remainder) const {
			if (divisor.empty()) ip_error("Division by zero\n");
			if (empty()) return 0;

			// If divisor fits into 32 bits, revert to short division
			uint32_t d32, rem;
			if (divisor.toUint32(d32)) {
				bignatural q = divide_by(d32, rem);
				remainder = rem;
				return q;
			}

			// If both dividend and divisor fit into 64 bits, use hardware division
			uint64_t n, d;
			if (toUint64(n) && divisor.toUint64(d)) {
				remainder = n % d;
				return n / d;
			}

			// If divisor is greater than dividend...
			if (divisor > *this) {
				remainder = *this;
				return 0;
			}

			// Use binary (per-bit) long division
			const bignatural& dividend = *this;

			bignatural quotient, loc_dividend;
			uint32_t next_dividend_bit = dividend.getNumSignificantBits();

			do {
				loc_dividend.push_bit_back(dividend.getBit(--next_dividend_bit));
				if (loc_dividend >= divisor) {
					loc_dividend = loc_dividend - divisor;
					quotient.push_bit_back(1);
				}
				else if (!quotient.empty()) quotient.push_bit_back(0);
			} while (next_dividend_bit);

			remainder = loc_dividend;

			return quotient;
		}

		// Greatest common divisor (Euclidean algorithm)
		inline bignatural GCD(const bignatural& D) const {
			bignatural A = *this;
			bignatural B = D;
			bignatural R;
			while (!A.empty() && !B.empty()) {
				A.divide_by(B, R);
				A = B;
				B = R;
			}
			if (A.empty()) return B;
			else return A;
		}

		// String representation in decimal form
		std::string get_dec_str() const {
			std::string st;
			bignatural N = *this;
			uint32_t R;
			if (N.empty()) return "0";
			while (!N.empty()) {
				N = N.divide_by(10, R);
				st += ('0' + R);
			}
			std::reverse(st.begin(), st.end());

			return st;
		}

		// String representation in binary form
		std::string get_str() const {
			std::string st;
			char s[33];
			s[32] = 0;
			for (uint32_t j = 0; j < m_size; j++) {
				for (int i = 0; i < 32; i++)
					s[i] = (digits[j] & (((uint32_t)1) << (31 - i))) ? '1' : '0';
				st += s;
			}
			return st;
		}

		// Count number of zeroes on the right (least significant binary digits)
		uint32_t countEndingZeroes() const {
			if (m_size == 0) return 0;
			uint32_t i = m_size - 1;
			uint32_t shft = 0;
			while (!digits[i]) {
				i--; shft += 32;
			}

			uint32_t s = UINT32_MAX;
			uint32_t m = digits[i];
			while ((s & m) == m) {
				s <<= 1;
				shft++;
			}
			return shft - 1;
		}

	protected:
		inline uint32_t& back() { return digits[m_size - 1]; }

		inline void pop_back() { m_size--; }

		inline uint32_t& operator[](int i) { return digits[i]; }

		inline void push_back(uint32_t b) {
			if (m_size == m_capacity) increaseCapacity(m_capacity * 2);
			digits[m_size++] = b;
		}

		inline void push_bit_back(uint32_t b) {
			if (m_size) {
				operator<<=(1);
				back() |= b;
			}
			else if (b) push_back(1);
		}

		inline void reserve(uint32_t n) {
			if (n > m_capacity) increaseCapacity(n);
		}

		inline void resize(uint32_t n) {
			reserve(n);
			m_size = n;
		}

		inline void fill(uint32_t v) {
			for (uint32_t i = 0; i < m_size; i++) digits[(int)i] = v;
		}

		inline void pop_front() {
			for (uint32_t i = 1; i < m_size; i++) digits[i - 1] = digits[i];
			pop_back();
		}

		// Count number of zeroes on the left (most significant digits)
		uint32_t countLeadingZeroes() const {
			uint32_t s = UINT32_MAX;
			const uint32_t m = digits[0];
			uint32_t shft = 0;
			while ((s & m) == m) {
				s >>= 1;
				shft++;
			}
			return shft - 1;
		}

		inline void pack() {
			while (m_size && digits[0] == 0) pop_front();
		}

		// a and b must NOT be this number!
		void toSum(const bignatural& a, const bignatural& b) {
			if (a.m_size == 0) operator=(b);
			else if (b.m_size == 0) operator=(a);
			else {
				const uint32_t a_s = a.m_size;
				const uint32_t b_s = b.m_size;
				uint32_t res_size = a_s;
				if (b_s > res_size) res_size = b_s;
				resize(res_size + 1);

				uint64_t carry = 0;
				for (uint32_t i = 1; i <= res_size; i++) {
					const uint64_t da = (i <= a_s) ? (a.digits[(int)(a_s - i)]) : (0);
					const uint64_t db = (i <= b_s) ? (b.digits[(int)(b_s - i)]) : (0);
					const uint64_t sm = da + db + carry;
					digits[(int)(res_size + 1 - i)] = sm & UINT32_MAX;
					carry = (sm >> 32);
				}
				digits[0] = (uint32_t)carry;
			}
			pack();
		}

		// a and b must NOT be this number!
		// Assume that b is smaller or equal than a!
		void toDiff(const bignatural& a, const bignatural& b) {
			if (b.m_size == 0) operator=(a);
			else {
				const uint32_t a_s = a.m_size;
				const uint32_t b_s = b.m_size;
				uint32_t res_size = a_s;
				if (b_s > res_size) res_size = b_s;
				resize(res_size);

				uint64_t debt = 0;
				for (uint32_t i = 1; i <= res_size; i++) {
					const uint64_t da = (i <= a_s) ? (a.digits[(int)(a_s - i)]) : (0);
					const uint64_t db = ((i <= b_s) ? (b.digits[(int)(b_s - i)]) : (0)) + debt;
					debt = !(da >= db);
					if (debt) digits[(int)(res_size - i)] = (uint32_t)((da + (((uint64_t)1) << 32)) - db);
					else digits[(int)(res_size - i)] = (uint32_t)(da - db);
				}
			}
			pack();
		}

		// a and b must NOT be this number!
		void toProd(const bignatural& a, const bignatural& b) {
			if (a.empty()) operator=(a);
			else if (b.empty()) operator=(b);
			else {
				resize(a.m_size + b.m_size);
				fill(0);

				uint32_t ls = 0;
				for (uint32_t i = b.m_size; i > 0; i--)
					a.addmul(b[(int)(i - 1)], ls++, *this);

				pack();
			}
		}

	private:

		// Multiplies by a single limb, left shift, and add to accumulator. Does not pack!
		void addmul(uint32_t b, uint32_t left_shifts, bignatural& result) const {
			uint64_t carry = 0;
			int d = (int)(result.m_size - m_size - left_shifts);
			for (uint32_t i = m_size; i > 0; i--) {
				uint64_t pm = ((uint64_t)digits[(int)(i - 1)]) * b + carry + result[(int)i + d - 1];
				result[(int)i + d - 1] = (uint32_t)pm;
				carry = pm >> 32;
			}
			result[d - 1] = (uint32_t)carry;
		}

		// Creates a copy of the array
		inline uint32_t* dupDigits() const {
			uint32_t *m = (uint32_t*)malloc(sizeof(uint32_t) * m_capacity);
			memcpy(m, digits, m_size * sizeof(uint32_t));
			return m;
		}

		void increaseCapacity(uint32_t new_capacity) {
			m_capacity = new_capacity;
			if (digits == m_static) {
				digits = (uint32_t*)malloc(sizeof(uint32_t) * m_capacity);
				for (uint32_t i = 0; i < static_digit_count; i++) digits[i] = m_static[i];
			}
			else 
				digits = (uint32_t*)realloc(digits, sizeof(uint32_t) * m_capacity);
		}

		friend class bigfloat;
	};


	/////////////////////////////////////////////////////////////////////
	// 	   
	// 	   B I G   F L O A T
	// 
	/////////////////////////////////////////////////////////////////////

	// A bigfloat is a floting point number with arbitrarily large mantissa.
	// In principle, we could have made the exponent arbitrarily large too,
	// but in practice this appears to be useless.
	// Exponents are in the range [-INT32_MAX, INT32_MAX]
	//
	// A bigfloat f evaluates to f = sign * mantissa * 2^exponent
	//
	// mantissa is a bignatural whose least significant bit is 1.
	// Number is zero if mantissa is empty.

	class bigfloat {
		bignatural mantissa; // .back() is less significant. Use 32-bit limbs to avoid overflows using 64-bits
		int32_t exponent; // In principle we might still have under/overflows, but not in practice
		int32_t sign;	// Redundant but keeps alignment

	public:
		// Default constructor creates a zero-valued bigfloat
		inline bigfloat() : sign(0), exponent(0) {}

		// Lossless conversion from double
		bigfloat(const double d) {
			sign = (d > 0) - (d < 0);

			if (sign) {
				uint64_t dn = *((uint64_t*)(&d));
				const uint64_t m = (dn & 0x000fffffffffffff) + 0x0010000000000000;
				mantissa.push_back(m >> 32);
				mantissa.push_back(m & 0x00000000ffffffff);
				dn <<= 1;
				dn >>= 53;
				exponent = ((int32_t)dn) - 1075; // Exp

				pack();
			}
			else exponent = 0;
		}

		// Truncated approximation
		double get_d() const {
			uint64_t dn = 0;
			if (mantissa.empty()) return 0.0;

			uint64_t m;
			int32_t e;
			uint32_t shft;

			if (mantissa.size() == 1) {
				m = ((uint64_t)mantissa[0]);
				shft = mantissa.countLeadingZeroes() + 21;
				m <<= shft;
				e = exponent - shft;
			}
			else {
				m = (((uint64_t)mantissa[0]) << 32) | ((uint64_t)mantissa[1]);
				e = exponent + 32 * ((uint32_t)mantissa.size() - 2);
				shft = mantissa.countLeadingZeroes();

				if (shft < 11) {
					m >>= (11 - shft);
					e += (11 - shft);
				}
				if (shft > 11) {
					m <<= (shft - 11);
					e -= (shft - 11);
					if (mantissa.size() > 2) m |= (mantissa[2] >> (43 - shft));
				}
			}
			m &= (~0x0010000000000000); // Remove implicit digit
			e += 52;

			if (e < (-1022)) return 0.0;
			if (e > 1023) return sign * INFINITY;

			if (sign < 0) dn |= 0x8000000000000000; // Set sign
			dn |= (((uint64_t)(e + 1023)) << 52); // Set exponent
			dn |= m; // Set mantissa

			return *((double*)(&dn));
		}
		
		bigfloat operator+(const bigfloat& b) const {
			if (mantissa.empty()) return b;
			if (b.mantissa.empty()) return *this;

			if (exponent == b.exponent) {
				bigfloat result;

				if (sign == b.sign) {
					result.mantissa.toSum(mantissa, b.mantissa);
					result.sign = sign;
				}
				else if (b.mantissa >= mantissa) {
					result.mantissa.toDiff(b.mantissa, mantissa);
					result.sign = b.sign;
				}
				else {
					result.mantissa.toDiff(mantissa, b.mantissa);
					result.sign = sign;
				}

				result.exponent = exponent;
				result.pack();
				return result;
			}
			else if (exponent > b.exponent) {
				bigfloat op(*this);
				op.leftShift(exponent - b.exponent);
				return op + b;
			}
			else { // exponent < b.exponent
				bigfloat op(b);
				op.leftShift(b.exponent - exponent);
				return op + *this;
			}
		}

		bigfloat operator-(const bigfloat& b) const {
			if (mantissa.empty()) return b.inverse();
			if (b.mantissa.empty()) return *this;

			if (exponent == b.exponent) {
				bigfloat result;

				if (sign != b.sign) {
					result.mantissa.toSum(mantissa, b.mantissa);
					result.sign = sign;
				}
				else if (b.mantissa >= mantissa) {
					result.mantissa.toDiff(b.mantissa, mantissa);
					result.sign = -sign;
				}
				else {
					result.mantissa.toDiff(mantissa, b.mantissa);
					result.sign = sign;
				}

				result.exponent = exponent;
				result.pack();
				return result;
			}
			else if (exponent > b.exponent) {
				bigfloat op(*this);
				op.leftShift(exponent - b.exponent);
				return op - b;
			}
			else { // exponent < b.exponent
				bigfloat op(b);
				op.leftShift(b.exponent - exponent);
				return *this - op;
			}
		}


		bigfloat operator*(const bigfloat& b) const {
			if (mantissa.empty() || b.mantissa.empty()) return 0;

			// Left-shift operator with highest exponent
			if (exponent == b.exponent) {
				bigfloat result;
				result.mantissa.toProd(mantissa, b.mantissa);
				result.exponent = exponent;
				result.sign = sign * b.sign;
				result.leftShift(result.exponent - exponent);
				result.exponent *= 2;

				result.pack();
				return result;
			}
			else if (exponent > b.exponent) {
				bigfloat op(*this);
				op.leftShift(exponent - b.exponent);
				return op * b;
			} // exponent < b.exponent
			else {
				bigfloat op(b);
				op.leftShift(b.exponent - exponent);
				return op * *this;
			}
		}

		inline void invert() { sign = -sign; }

		bigfloat inverse() const {
			bigfloat r = *this;
			r.invert();
			return r;
		}

		inline int sgn() const { return sign; }

		std::string get_str() const {
			std::string s;
			if (sign == 0) s += "0";
			if (sign < 0) s += "-";
			s += mantissa.get_str();
			s += " * 2^";

			char es[1024];
			s += _itoa(exponent, es, 10);
			return s;
		}

		const bignatural& getMantissa() const { return mantissa; }
		int32_t getExponent() const { return exponent; }

	private:

		// Right-shift as long as the least significant bit is zero
		void pack() {
			while (!mantissa.empty() && mantissa.back() == 0) {
				mantissa.pop_back();
				exponent += 32;
			}

			if (mantissa.empty()) {
				sign = exponent = 0;
				return;
			}

			const uint32_t s = mantissa.countEndingZeroes();
			if (s) {
				for (int i = (int)mantissa.size() - 1; i > 0; i--) {
					mantissa[i] >>= s;
					mantissa[i] |= (mantissa[i - 1] << (32 - s));
				}
				mantissa[0] >>= s;
				exponent += s;
			}

			mantissa.pack();
		}

		// Left-shift the mantissa by n bits and reduce the exponent accordingly
		void leftShift(uint32_t n) {
			mantissa <<= n;
			exponent -= n;
		}
	};

	inline int sgn(const bigfloat& f) {
		return f.sgn();
	}

	inline bigfloat operator-(const bigfloat& f) {
		return f.inverse();
	}

	inline bigfloat operator*(double d, const bigfloat& f) {
		return f * bigfloat(d);
	}


/////////////////////////////////////////////////////////////////////
// 	   
// 	   B I G   R A T I O N A L
// 
/////////////////////////////////////////////////////////////////////

// A bigrational is a fraction of two bignaturals with a sign.
// Number is zero if sign is zero

	class bigrational {
		bignatural numerator, denominator;
		int32_t sign;	// Redundant but keeps alignment

		bigrational(const bignatural& num, const bignatural& den, int32_t s) :
			numerator(num), denominator(den), sign(s) {
			canonicalize();
		}

		void invert() { std::swap(numerator, denominator); }

		bigrational inverse() const { bigrational r = *this; r.invert(); return r; }

	public:
		bigrational() : sign(0) {}

		bigrational(const bigfloat& f) {
			if (f.sgn() == 0) sign = 0;
			else {
				sign = f.sgn();
				numerator = f.getMantissa();
				denominator = 1;
				int32_t e = f.getExponent();
				if (e >= 0) numerator <<= e;
				else denominator <<= (-e);
			}
		}

		void compress() {
			const uint32_t nez = numerator.countEndingZeroes();
			const uint32_t dez = denominator.countEndingZeroes();
			const uint32_t s = std::min(nez, dez);
			numerator >>= s;
			denominator >>= s;
		}

		void canonicalize() {
			if (sign) {
				if (numerator.empty()) {
					numerator = denominator = 0;
					sign = 0;
				}
				else {
					compress();
					bignatural r;
					const bignatural gcd = numerator.GCD(denominator);
					numerator = numerator.divide_by(gcd, r);
					denominator = denominator.divide_by(gcd, r);
				}
			}
		}

		bigrational operator*(const bigrational& r) const {
			if (sign == 0 || r.sign == 0) return bigrational();
			else return bigrational(numerator * r.numerator, denominator * r.denominator, sign * r.sign);
		}

		bigrational operator/(const bigrational& r) const {
			if (!r.sign) ip_error("bigrational::operator/ : division bby zero!\n");
			return operator*(r.inverse()); 
		}

		bigrational operator+(const bigrational& r) const {
			if (sign == 0) return r;
			else if (r.sign == 0) return *this;
			else {
				const bignatural left_num = numerator * r.denominator;
				const bignatural right_num = r.numerator * denominator;
				if (sign > 0 && r.sign > 0)	return bigrational(left_num + right_num, denominator * r.denominator, 1);
				else if (sign < 0 && r.sign < 0) return bigrational(left_num + right_num, denominator * r.denominator, -1);
				else if (sign > 0 && r.sign < 0) {
					if (left_num >= right_num) return bigrational(left_num - right_num, denominator * r.denominator, 1);
					else return bigrational(right_num - left_num, denominator * r.denominator, -1);
				}
				else /*if (sign < 0 && r.sign > 0)*/ {
					if (left_num >= right_num) return bigrational(left_num - right_num, denominator * r.denominator, -1);
					else return bigrational(right_num - left_num, denominator * r.denominator, 1);
				}
			}
		}

		std::string get_dec_str() const {
			std::string st;
			if (sign < 0) st = "-";
			else if (sign == 0) return "0";
			st += numerator.get_dec_str();
			const std::string dens = denominator.get_dec_str();
			if (dens != "1") {
				st += "/";
				st += dens;
			}
			return st;
		}

		std::string get_str() const {
			std::string st;
			if (sign < 0) st = "-";
			else if (sign == 0) return "0";
			st += numerator.get_str();
			st += "/";
			st += denominator.get_str();
			return st;
		}
	};

#endif //NUMERICS_H