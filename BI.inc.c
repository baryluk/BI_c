#ifndef __BI_C_H__
#define __BI_C_H__ 1

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#ifndef KARATSUBA_MUL_P
#define KARATSUBA_MUL_P 64
#endif
#ifndef KARATSUBA_SQR_P
#define KARATSUBA_SQR_P 64
#endif

#define PRESHIFT 1

#ifndef DIV_LIMIT
#define DIV_LIMIT 4
#endif


#ifndef BI_RABIN_MILLER_TESTS
#define BI_RABIN_MILLER_TESTS 6
#endif


/* realloc array only if it is more than REALLC_TRESH bigger than optimal */
#define BI_REALLOC_TRESH 32


#define _ do { fprintf(stderr, "%s:%d in %s\n",__FILE__,__LINE__,__func__); } while (0);
#define __(x) do { fprintf(stderr, "%s:%d in %s  %s\n",__FILE__,__LINE__,__func__,##x); } while (0);
/*#define DD(...) do { fprintf(stderr, __VA_ARGS__); } while (0);*/
#define DD(...) do {  } while (0);


/*** Obsluga duzych liczba */


#define max(_a_, _b_) ((_a_) > (_b_) ? (_a_) : (_b_))
#define min(_a_, _b_) ((_a_) < (_b_) ? (_a_) : (_b_))


/*
#define BIBit (32)
#define BIByte (4)
#define BIType uint32_t
#define BIMask (~(0xffffffffu<<BIBit))
*/

#define BIBit (16)
#define BIByte (2)
#define BIType uint16_t
#define BIMask (~(0xffffffffu<<BIBit))

struct BI_ {
	BIType *b;
	int l;
	int sign;
	int ref_count;
	int radix_point;
};
typedef struct BI_ BI;

/* should work in ANSI mode and C99 */
#define BIsmallconst(x) ( (struct BI_) { .b=(BIType[]) { (x) } , .l=1 , .sign=1 } )

BI BInew(int l_) {
	BI z;
	z.l = l_;
	z.b = malloc(sizeof(BIType)*l_);
	if (z.b == NULL) {
		perror("malloc");
		exit(1);
	}
	z.sign = 1;
	z.ref_count = 1;
	z.radix_point = 0;
	return z;
}

void BIfree(BI x) {
	free(x.b);
}

int BImsb(BI x) {
	int msb = (x.l-1)*BIBit;
	BIType v = x.b[x.l-1];

	/* fix last bits */

	/* slow method */
	/*
	while (v) {
		v >>= 1;
		msb++;
	}
	return msb;
	*/


/* faster method 32 bit */
/*
	uint32_t r =     (v > 0xFFFF) << 4; v >>= r;
	uint32_t shift = (v > 0xFF  ) << 3; v >>= shift; r |= shift;
	         shift = (v > 0xF   ) << 2; v >>= shift; r |= shift;
	         shift = (v > 0x3   ) << 1; v >>= shift; r |= shift;
	                                                 r |= (v >> 1);
*/

/* faster method 16 bit */
	uint32_t r     = (v > 0xFF  ) << 3; v >>= r;
	uint32_t shift = (v > 0xF   ) << 2; v >>= shift; r |= shift;
	         shift = (v > 0x3   ) << 1; v >>= shift; r |= shift;
	                                                 r |= (v >> 1);

	/*return msb + r; */
	return msb | r;
}

BI BIclone(BI x) {
	BI z = BInew(x.l);
	z.sign = x.sign;
	memcpy(z.b, x.b, sizeof(BIType)*x.l);
	return z;
}

int BIeq(BI x, BI y) {
	if ((x.l != y.l) || (x.sign != y.sign)) {
		return 0;
	}
	if (memcmp(x.b, y.b, sizeof(BIType)*x.l) == 0) {
		return 1;
	} else {
		return 0;
	}
}


int is_zero(BI x) {
	return (x.l == 1) && (x.b[0] == 0);
}

int BIcmp(BI x, BI y) {
	int i, j;
	int t;
	/* check signs */
	if (x.sign > 0 && y.sign > 0) {
		/* check lengths of numbers */
		if (x.l > y.l) {
			return 1;
		}
		if (x.l < y.l) {
			return -1;
		}
		/* TODO: memcmp */
		/* check digits of numbers, starting from most significant, and terminate early */
		for (i = 0, j = x.l-1; i < x.l; i++, j--) {
			t = (int32_t)x.b[j] - (int32_t)y.b[j];
			if (t) {
				return t;
			}
		}
		return 0;
	}
	assert(0);
}

BI BIfix(BI x, int l_) {
	assert(l_ < 1000000);
	if (l_ > x.l) {
		exit(2);
	}

	/* reduce back leading zeros */
	while (l_ && !x.b[l_]) {
		l_--;
	}
	l_++;

	if (l_ < 1) {
		l_ = 1;
	}
	if (l_ > x.l) {
		exit(52);
	}
	if (x.l > l_ + BI_REALLOC_TRESH) {
		DD("reallok\n");
		x.b = realloc(x.b, sizeof(BIType)*l_);
		if (x.b == NULL) {
			perror("realloc");
			exit(3);
		}
	}
	x.l = l_;
	return x;
}


BI BI_0, BI_1, BI_2, BI_3, BI_5, BI_7, BI_11, BI_13, BI_17, BI_23, BI_31, BI_37, BI_61, BI_10000;


/* O(1) */
BI BIzero() {
	BI z = BInew(1);
	z.b[0] = 0;
	z.sign = 1;
	return z;
}

/* O(1) */
BI BIone() {
	BI z = BInew(1);
	z.b[0] = 1;
	z.sign = 1;
	return z;
}

/* skonweruj "mala" liczbe do BI */
/* O(1) */
BI BIsmall(uint64_t x) {
	BI z = BInew(4);
	uint64_t c = x;
	z.b[0] = c & BIMask;
	c >>= BIBit;
	z.b[1] = c & BIMask;
	c >>= BIBit;
	z.b[2] = c & BIMask;
	c >>= BIBit;
	z.b[3] = c & BIMask;
	c >>= BIBit;
	z.sign = 1;
	return BIfix(z, 4-1);
}

/* skonwertuj bardzo mala liczbe, mniejsza niz 65535 */
/* O(1) */
BI BIvsmall(BIType x) {
	BI z = BInew(1);
	uint64_t c = x;
	z.b[0] = c & BIMask;
	z.sign = 1;
	return z;
}



/* O(n) */
BI BIshift(BI x, int n) {
	int i, j;
	int offset;
	int bit_offset_a;
	int bit_offset_b;
	BI z;
	assert(x.sign == 1);
	if (n > x.l*BIBit) {
		return BIzero();
	}
	if (n == 0) {
		return BIclone(x);
	}
	if (n > 0) { /* right shift */
		z = BInew(x.l + ((-n+1)/BIBit));
		offset = n/BIBit;
		bit_offset_a = n-(offset*BIBit);
		bit_offset_b = BIBit-bit_offset_a;
		z.b[z.l-1] = 0;
		if (bit_offset_a == 0) {
			for (i = 0, j = offset; i < z.l; i++, j++) {
				z.b[i] = x.b[j];
			}
		} else {
			for (i = 0, j = offset; (i < z.l) && (j+1 < x.l); i++, j++) {
				z.b[i] = (x.b[j+1] << bit_offset_b) | (x.b[j] >> bit_offset_a);
			}
			z.b[z.l-1] |= (x.b[x.l-1] >> bit_offset_a);
		}
	} else { /* left shift */
		z = BInew(x.l + ((-n+BIBit-1)/BIBit));
		offset = (-n)/BIBit;
		bit_offset_a = (-n)-(offset*BIBit);
		bit_offset_b = BIBit-bit_offset_a;
		for (i = 0; i < offset; i++) {
			z.b[i] = 0;
		}
		z.b[z.l-1] = 0;
		if (bit_offset_a == 0) {
			for (i = offset, j = 0; j < x.l; i++, j++) {
				z.b[i] = x.b[j];
			}
		} else {
			z.b[offset] = (x.b[0] << bit_offset_a);
			for (i = offset+1, j = 1; (i < z.l) && (j < x.l); i++, j++) {
				z.b[i] =  (x.b[j-1] >> bit_offset_b) | (x.b[j] << bit_offset_a);
			}
			z.b[z.l-1] |= (x.b[x.l-1] >> bit_offset_b);
		}
	}
	return BIfix(z,z.l-1);
}



/* O(n) */
BI BIadd(BI x, BI y) {
	BI z;
	int i;
	int m;
	uint32_t temp;

	assert(x.sign == 1);
	assert(y.sign == 1);

	/* adding zero */
	if (y.l == 1 && y.b[0] == 0) {
		return BIclone(x);
	}
	if (x.l == 1 && x.b[0] == 0) {
		return BIclone(y);
	}

	/* swap to make x bigger than y */
	if (x.l < y.l) {
		z = x;
		x = y;
		y = z;
	}


	z = BInew(x.l+1); /* max(x.l, y.l) */
	m = y.l; /* min(x.l, y.l) */

	/* perform addition, from right to left, carring carry in temp */
	temp = 0;
	for (i = 0; i < m; i++) {
		temp += (uint32_t)x.b[i] + (uint32_t)y.b[i];
		z.b[i] = (temp & BIMask);
		temp >>= BIBit;
	}

	/* perform rest of addition in x (which is bigger than y) */
	for (; i < x.l; i++) {
		temp += x.b[i];
		z.b[i] = (temp & BIMask);
		temp >>= BIBit;
	}

	/* write most significant digit (possible if carry spiled over most significant digit in x) */
	z.b[i] = (temp & BIMask);

	return BIfix(z, i);
}

/* todo: add(x, y, offset)  */
/* todo: sub(x, y, offset)  */

/* O(n) */
BI BIsub(BI x, BI y) {
	BI z;
	int i;
	int m;
	uint32_t temp;
	uint32_t borrow;

	assert(x.sign == 1);
	assert(y.sign == 1);
	i = BIeq(x, y);
	if (i == 1) {
		return BIzero();
	}
/*
	i = BIcmp(x, y);
	assert(i >= 0);
*/

	/* substracting */

	if (y.l == 1 && y.b[0] == 0) {
		return BIclone(x);
	}

#if 0
	/* swap if result will be negative */
	if (i < 0) {
		z = x;
		x = y;
		y = z;
		z.sign = -1;
	}
#endif

	z = BInew(x.l); /* max(x.l, y.l) */
	m = y.l; /* min(x.l, y.l) */

	/* perform substraction, from right to left, carring borrow */
	borrow = 0;
	temp = 0;
	for (i = 0; i < m; i++) {
		if ((uint32_t)x.b[i] >= borrow+(uint32_t)y.b[i]) {
			temp = (uint32_t)x.b[i] - (uint32_t)y.b[i] - borrow;
			borrow = 0;
		} else {
			temp = (1<<BIBit) + (uint32_t)x.b[i] - (uint32_t)y.b[i] - borrow;
			borrow = 1;
		}
		z.b[i] = (temp & BIMask);
	}

	/* last borrowing from x */
	if (borrow) {
		/* we assume borrow is always equal 1 in this loop, and break ASAP it is 0 */
		for (; i < x.l; i++) {
			if ((uint32_t)x.b[i] >= 1) {
				temp = (uint32_t)x.b[i] - /*borrow*/ 1;
				z.b[i] = (temp & BIMask);
				/*borrow = 0;*/ /* we will not use borrow anymore */
				i++;
				break;
			} else {
				temp = (1<<BIBit) + (uint32_t)x.b[i] - /*borrow*/ 1;
				/*borrow = 1;*/ /* borrow is already 1 */
				z.b[i] = (temp & BIMask);
			}
		}
	}

	/* rest is without borrows */
	for (; i < x.l; i++) {
		z.b[i] = x.b[i];
	}

	return BIfix(z, i-1);
}


/* O(n) */
void toadd(BI *x, BI y) {
	BI z;
	z = BIadd(*x, y);
	BIfree(*x);
	*x = z;
}

/* O(n) */
void tosub(BI *x, BI y) {
	BI z;
	z = BIsub(*x, y);
	BIfree(*x);
	*x = z;
}

/* return integer composed of digits of x beetwen a and b
 *
 * Note: returned integer can still need fixing leading zeros to be correct
 *
 * Do not free returned BI, just free underlaying integer x.
 *
 * O(1)
 */
BI BIslice(BI x, uint32_t a, uint32_t b) {
	BI r;
	assert(0 <= a);
	assert(a < b);
	assert(b <= x.l);
	r.sign = x.sign;
	r.b = x.b + a;
	r.l = (b-a);
	r.ref_count = 1;
	r.radix_point = 0;
	return r;
}


/* base case, just like school book but column oriented */
/* O(n^2) */
BI BImul_basecase_columns(BI x, BI y) {
	BI z;
	int n;
	int i, j;
	uint64_t temp;
	int a, b, top;

	if (x.l == 1 && x.b[0] == 0) {
		return BIzero();
	}
	if (y.l == 1 && y.b[0] == 0) {
		return BIzero();
	}

	n = max(x.l, y.l)+1;
	z = BInew(x.l+y.l+1);

	z.sign = x.sign*y.sign; /* or us sign ^ sign, but currently we have -1, +1 */

	temp = 0;
	for (i = 0; i < z.l-1; i++) {
		a = min(i, n-1);
		b = max(0, i+1-n);
		top = max(a,b);
		for (j = min(a,b); j <= top; j++) {
			if (j >= 0 && j < x.l && (i >= j) && (i < j+y.l)) {
				temp += ((uint32_t)x.b[j])*y.b[i-j];
			}
		}
		z.b[i] = (temp & BIMask);
		temp >>= BIBit;
	}
	z.b[i] = (temp & BIMask);
	return BIfix(z, i);
}

/* forward declaration */
BI BImul(BI x, BI y);

int indent = 0;
void I() {
	int i;
	for (i = 0; i < indent; i++) printf("    ");
}


/* O(n^1.585)
 *
 * Based on simple equation (cf. Wikipedia):
 *
 * z = (a*x1 + x0)(a*y1 + y0)
 * z = z2*a*a + z1*a + z0
 *
 * where z0 = x0*y0
 *       z1 = x1*y1
 *       z2 = (x0*y1 + x1*y) === (x0+x1)(y0+y1) - z0 - z1
 *
 */
BI BImul_karatsuba(BI x, BI y) {
	BI x1, x0;
	BI y1, y0;
	BI z;
	BI z0;
	BI z1;
	BI z2;
	BI t, tx, ty;

	assert(x.sign == 1);
	assert(y.sign == 1);

	if (x.l <= BI_KARATSUBA_MUL_P || y.l <= BI_KARATSUBA_MUL_P) {
		z = BImul_basecase_columns(x, y);
		return z;
	}

	/* make x larger, by swaping */
	if (x.l < y.l) {
		t = y;
		y = x;
		x = t;
	}

	int nn = x.l/2;
	/* x is larger */
	if (nn >= y.l) {
		/* unbalanced. but perform division of operands, so in deeper recursion level, it can become balanced.  y1 == 0 */
		y0 = y;

		x0 = BIslice(x, 0, nn);
		while ((x0.l > 1) && !x0.b[x0.l-1]) { /* fix trailing zeros possible in x0 */
			x0.l--;
		}

		z0 = BImul(x0, y0);

		x1 = BIslice(x, nn, x.l);

		z1 = BImul(x1, y0);

		/* x0, x1, y0, could be freeded here, but they are just aliases */
		t = BIshift(z1, -BIBit*nn);
		BIfree(z1);

		z = t;
		toadd(&z, z0);
		BIfree(z0);
	} else {
		x0 = BIslice(x, 0, nn);
		while ((x0.l > 1) && !x0.b[x0.l-1]) { /* fix trailing zeros possible in x0 */
			x0.l--;
		}

		y0 = BIslice(y, 0, nn);
		while ((y0.l > 1) && !y0.b[y0.l-1]) { /* fix trailing zeros possible in y0 */
			y0.l--;
		}

		z0 = BImul(x0, y0);

		x1 = BIslice(x, nn, x.l);

		y1 = BIslice(y, nn, y.l);

		z2 = BImul(x1, y1);

		t = BIadd(z0, z2);

		tx = BIadd(x1, x0);
		ty = BIadd(y1, y0);

		/* x0, x1, y0, y1, could be freeded here, but they are just aliases */
		z1 = BImul(tx, ty);
		BIfree(tx);
		BIfree(ty);

		tosub(&z1, t);
		BIfree(t);

		/* this add, add, can be done much more efficient, remembering that we really do not need to perform shift, and only perform partial addition */

		/*    (z2 << 2n) + (z1 << n) + z0 */
		/* == ((z2 << n) + z1) << n) + z0 */

		z = BIshift(z2, -BIBit*nn);
		BIfree(z2);
		toadd(&z, z1);
		BIfree(z1);
		t = BIshift(z, -BIBit*nn);
		BIfree(z);
		z = t;

		toadd(&z, z0);
		BIfree(z0);
	}
	return BIfix(z, z.l-1);
}



/* O(n^1.585) */
BI BImul(BI x, BI y) {
	return BImul_karatsuba(x, y);
}


/* O(n^2) */
BI BIsqr_basecase_columns(BI x) {
	BI z;
	int n;
	int i, j;
	uint64_t temp;
	int a, b, top;

	n = x.l+1;
	z = BInew(2*x.l+1);
	z.sign = 1;

	temp = 0;
	for (i = 0; i < z.l-1; i++) {
		a = min(i, n-1);
		b = max(0, i+1-n);
		top = max(a,b);
		for (j = min(a,b); j <= top; j++) {
			if (j >= 0 && j < x.l && (i >= j) && (i < j+x.l)) {
				temp += ((uint32_t)x.b[j])*x.b[i-j];
			}
		}
		z.b[i] = (temp & BIMask);
		temp >>= BIBit;
	}
	z.b[i] = (temp & BIMask);
	return BIfix(z, i);
}

BI BIsqr(BI x);

BI BIsqr_karatsuba(BI x) {
	BI x1, x0;
	BI z;
	BI z0;
	BI z1;
	BI z2;
	BI t, tx;

	assert(x.sign == 1);

	if (x.l <= BI_KARATSUBA_SQR_P) {
		z = BIsqr_basecase_columns(x);
		return z;
	}

	int nn = x.l/2;

		x0 = BIslice(x, 0, nn);
		while ((x0.l > 1) && !x0.b[x0.l-1]) { /* fix trailing zeros possible in x0 */
			x0.l--;
		}

		z0 = BIsqr(x0);

		x1 = BIslice(x, nn, x.l);

		z2 = BIsqr(x1);

		t = BIadd(z0, z2);

		tx = BIadd(x1, x0);

		/* x0, x1, y0, y1, could be freeded here, but they are just aliases */
		z1 = BIsqr(tx);
		BIfree(tx);

		tosub(&z1, t);
		BIfree(t);

		/* this add, add, can be done much more efficient, remembering that we really do not need to perform shift, and only perform partial addition */

		/*    (z2 << 2n) + (z1 << n) + z0 */
		/* == ((z2 << n) + z1) << n) + z0 */

		z = BIshift(z2, -BIBit*nn);
		BIfree(z2);
		toadd(&z, z1);
		BIfree(z1);
		t = BIshift(z, -BIBit*nn);
		BIfree(z);
		z = t;

		toadd(&z, z0);
		BIfree(z0);

	return BIfix(z, z.l-1);
}

/* O(N^1.585) */
BI BIsqr(BI x) {
	return BIsqr_karatsuba(x);
}

void tomul(BI *x, BI y) {
	BI z;
	z = BImul(*x, y);
	BIfree(*x);
	*x = z;
}

/* division of positive x and y */
/* x = y*q + r */
/* start from:
   x = 0*q + r; where r==x
 *
 * O(n*m)
 *
 * School book division, bit by bit, own idea and optimisations.
 */
void BIdivmod_basecase_binary(BI x, BI y, BI *q_, BI *r_) {
	int shift;
	int test;
	BI temp, temp2;
	BI q, r;
#if PRESHIFT
	BI shifts[16];
	int i;
#else
	BI bi_shift;
#endif

	/* todo: eventually shifting of y, can be done in place, or in the same buffer,
	 * without allocating, and allocating
	 */

	assert(x.sign == 1);
	assert(y.sign == 1);

	if (is_zero(y)) {
		fprintf(stderr, "aa division by zero\n");
		*q_ = BIzero();
		*r_ = BIzero();
		exit(14);
		return;
	}
	assert(!is_zero(y));

	if (x.l < y.l) {
		*q_ = BIzero();
		*r_ = BIclone(x);
		return;
	}


	/*q = BIzero();*/
	q = BInew(x.l-y.l+2);
	memset(q.b, 0, sizeof(BIType)*q.l);

	r = BIclone(x);

#if PRESHIFT
	/*
		precalculate all shifts modulo 16 of y, 0, 1, 2, ...., 15.
		then take one of them with based on shift mod 16, and shift pointer to b
		more importantly this table can be shared beetween calls for the same y
		(so in example in BImod(., n), for the same n).
	*/
	shift = BIBit*(x.l - y.l); /* slightly bigger than biggest possible shift in main loop */
	for (i = 0; i < 16; i++) {
		shifts[i] = BIshift(y, -shift-i);
	}
#endif

	/*test = 1;*/ /* possible that test == 0, is actually correct, but it is safe to assume on entry a 1 */
	test = BIcmp(r, y);
	while (test >= 0) {
		if (test == 0) {
			BIfree(r);
			r = BIzero();   /* r==y => temp=y => tosub(&r, temp) == tosub(&r, r) => r=0 */
			q.b[0] |= 1;    /* r==y  => shift == 0 */
			break;
		}

		shift = BImsb(r) - BImsb(y);

#if PRESHIFT
		temp.sign = 1;
		temp.l = y.l + (shift/BIBit) + (shifts[shift%BIBit].l-shifts[0].l);
		temp.b = shifts[shift%BIBit].b + (shifts[shift%BIBit].l - y.l - (shift/BIBit)) - (shifts[shift%BIBit].l-shifts[0].l);
#else
		temp = BIshift(y, -shift);
#endif

		if (BIcmp(temp, r) > 0) {
			shift--;

#if PRESHIFT
			temp2.sign = 1;
			temp2.l = y.l + (shift/BIBit) + (shifts[shift%BIBit].l-shifts[0].l);
			temp2.b = shifts[shift%BIBit].b + (shifts[shift%BIBit].l - y.l - (shift/BIBit)) - (shifts[shift%BIBit].l-shifts[0].l);
			temp = temp2;
#else
			temp2 = BIshift(temp, 1);
			BIfree(temp);
			temp = temp2;
#endif
		}

		tosub(&r, temp);
#if !PRESHIFT
		BIfree(temp);
#endif

		/* this is equivalent in our case to doing toadd(q, BIshift(BI_1, -shift)), + cleaning, but much faster */
		q.b[shift/BIBit] |= (1 << (shift&(BIBit-1)));

		test = BIcmp(r, y);
	}
	*q_ = BIfix(q,q.l-1);
	*r_ = r;

#if PRESHIFT
	for (i = 0; i < 16; i++) {
		BIfree(shifts[i]);
	}
#endif

	return;
}


BI BI_beta;

/** O(n*m) */
void BIdivmod_basecase_c_sub(BI x, BI y, BI *q_, BI *r_) {
	BI temp;
	BI temp2, qq;
	BI t;
	uint32_t q;
	BI addition;
	int free_x = 0;
	int free_temp2 = 0;

	assert(!is_zero(y));

	assert(x.l <= y.l+1);
	assert(BImsb(y) == y.l*BIBit-1); /* leading bit set */

	assert(y.b[y.l-1] & (1 << (BIBit-1))); /* leading bit set */

	if (x.l < y.l) {
		*q_ = BIzero();
		*r_ = BIclone(x);
		return;
	}

	temp.sign = x.sign;
	temp.b = x.b + 1;
	temp.l = x.l - 1;

	addition = BIzero();
	while (BIcmp(temp, y) >= 0) { /* while (x >= BI_beta*y) { */
		temp2 = BIshift(y, -BIBit);
		free_temp2 = 1;
		temp = BIsub(x, temp2);
		if (free_x) {
			BIfree(x);
		}
		x = temp;
		free_x = 1;

		if (x.l < y.l) {
			*q_ = BIzero();
			*r_ = BIclone(x);
			return;
		}

		temp.sign = x.sign;
		temp.b = x.b + 1;
		temp.l = x.l - 1;

		toadd(&addition, BI_beta);
		return;
	}
	if (free_temp2) {
		BIfree(temp2);
	}

	/** asserts should still be correct */
	assert(x.l <= y.l+1);
	assert(BImsb(y) == y.l*BIBit-1); /* leading bit set */
	assert(y.b[y.l-1] & (1 << (BIBit-1))); /* leading bit set */

	/*assert(y.l == x.l-1);*/
	/*assert(x.l >= 2);*/
	/*q = (((uint32_t)(x.b[y.l])) << BIBit) + (uint32_t)(x.b[y.l-1]); */
	q = 0;
	if (x.l == y.l+1) {
		q += (((uint32_t)(x.b[y.l])) << BIBit);
	}
	if (x.l >= y.l) {
		q += (uint32_t)(x.b[y.l-1]);
	}
	q /= y.b[y.l-1];
	if (q >= (1<<BIBit)) { /* can be tested using mask */
		q = (1<<BIBit)-1;
	}
	qq = BIvsmall(q);
	t = BImul_basecase_columns(y, qq);
	BIfree(qq);
	if (BIcmp(t, x) > 0) { /* too large, fix */
		q--;
		tosub(&t, y);
		if (BIcmp(t, x) > 0) { /* still too large, fix again */
			q--;
			tosub(&t, y);
			assert(!(BIcmp(t, x) > 0)); /* there should be no fixes needed again */
		}
	}
	*q_ = BIvsmall(q);
	if (!is_zero(addition)) {
		toadd(q_, addition);
	}
	BIfree(addition);
	*r_ = BIsub(x, t);
	if (free_x) {
		BIfree(x);
	}
	BIfree(t);

	return;
}

/* O(n*m)
 */
void BIdivmod_basecase_c(BI x, BI y, BI *q_, BI *r_) {
	int shift;
	BI temp;
	int c;
	BI xprim, s;
	BI qprim, rprim;

	assert(x.sign == 1);
	assert(y.sign == 1);

	assert(!is_zero(y));

	if (x.l < y.l) {
		*q_ = BIzero();
		*r_ = BIclone(x);

		return;
	}

	shift = BImsb(y);
	/* shift y, so it have 1-bit set at the most significant bit at most signicifcant digit
	 * shift x, by the same amount, to make x/y still the same value
	 * if it already is correct, do nothing,
	 * if we need to shift, remember about freeing x, and y at the return
	 */
	if ((shift % 16) != 15) {
		x = BIshift(x, (shift % 16) - 15);
		y = BIshift(y, (shift % 16) - 15);
		assert(!(y.l == 1 && y.b[0] == 0));
	}

	if (x.l == y.l) {
		c = BIcmp(x, y);
		if (c < 0) {
			*q_ = BIzero();
			*r_ = BIclone(x);
		} else {
			*q_ = BIone();
			*r_ = BIsub(x, y);
		}

		goto ending_unshift;
	}
	if (x.l == y.l+1) {
		/* helper call */
		BIdivmod_basecase_c_sub(x, y, q_, r_);

		goto ending_unshift;
	}
	assert(x.l > y.l+1);

	/* split x at x.l-y.l-1 position */

	xprim.sign = x.sign;
	xprim.b = x.b + (x.l-y.l-1);
	xprim.l = y.l+1;  /* == x.l - (x.l-y.l-1); */
	s.sign = x.sign;
	s.b = x.b;
	s.l = (x.l-y.l-1);

	/* helper call */
	BIdivmod_basecase_c_sub(xprim, y, &qprim, &rprim);

	temp = BInew(s.l + rprim.l + 1); /* one more so we would reuse it later */
	temp.l--;
	memcpy(temp.b, s.b, sizeof(BIType)*s.l);
	memcpy(temp.b + s.l, rprim.b, sizeof(BIType)*rprim.l);
	BIfree(rprim);

	/* recursion */
	BIdivmod_basecase_c(temp, y, q_, r_);

	assert(q_->l <= s.l);

	temp.l++; /* we can do this because we allocated one sparse digit */

	memcpy(temp.b, q_->b, sizeof(BIType)*(q_->l));
	memset(temp.b + q_->l, 0, sizeof(BIType)*(s.l - q_->l));
	assert(qprim.l + s.l <= temp.l);
	memcpy(temp.b + s.l, qprim.b, sizeof(BIType)*qprim.l);
	temp.l = qprim.l + s.l;

	BIfree(*q_);
	BIfree(qprim);

	*q_ = temp;

ending_unshift:

	if ((shift % 16) != 15) {
		/* todo: assert that last 15 - (shift % 16) bits are equal zero */
		temp = BIshift(*r_, 15 - (shift % 16)); /* shift back reminder */
		BIfree(*r_);
		*r_ = temp;
	} else {
		*r_ = *r_;
	}


	if ((shift % 16) != 15) {
		BIfree(x);
		BIfree(y);
	}

	return;
}

void BIdivmod_basecase(BI x, BI y, BI *q_, BI *r_) {
/*	BIdivmod_basecase_c(x, y, q_, r_); */
	BIdivmod_basecase_binary(x, y, q_, r_);
}


/**
 * See "Fast Recursive Division", Christoph Burnikel, Joachim Ziegler,
 * MPI-I-98-1-022, October 1998.
 *
 * "Fast Division of Large Integers. A Comparison of Algorithms", Karl Hasselstrom, Master Thesis.
 *
 */


/* forward declaration for co-recursion */
void BIdivmod_divide_conquere_32(BI x, BI y, BI *q_, BI *r_);

void BIdivmod_divide_conquere_42(BI x, BI y, BI *q_, BI *r_) {
	BI x1, x2, x3, x4;
	BI y1, y2;
	BI x123;
	int n = y.l;
	int n2 = n/2;
	int i;
	BI t, r1, q1, q2;

/*	printf("subdziele42 %d przez %d\n", x.l, y.l); */

	if (x.l < y.l) {
		*q_ = BIzero();
		*r_ = BIclone(x);
		return;
	}

	if (!(y.l & 1) || y.l < BI_DIV_LIMIT) {
/*		printf("base case\n"); */
		BIdivmod_basecase(x, y, q_, r_);
		return;
	}


	assert((n & 1) == 0);/* even */
	assert(n > 6); /* sufficiently big */
	assert(!is_zero(y));
	assert(x.l <= 4*n2);
	if (!(x.l > 3*n2)) {
		BIdivmod_divide_conquere_32(x, y, q_, r_);
		return;
	}
	assert(x.l > 3*n2);
	assert(y.b[y.l-1] & (1 << (BIBit-1))); /* leading bit in y set */

	/* x = [x1, x2, x3, x4] */
	x4 = BIslice(x, 0, n2);
	x3 = BIslice(x, n2, n);
	x2 = BIslice(x, n, n+n2);
	x1 = BIslice(x, n+n2, x.l);

	/* y = [y1, y2] */
	y2 = BIslice(y, 0, n2);
	y1 = BIslice(y, n2, y.l);

	x123 = BIslice(x, n2, x.l); /* 3 most significant "digits", x1, x2, x3 */

	/* divide [x1, x2, x3] by [y1, y2] */
	BIdivmod_divide_conquere_32(x123, y, &q1, &r1);

	/* construct [r1;1, r1;2, x4] */
	t = BIshift(r1, -n2*BIBit);
	BIfree(r1);
	for (i = 0; i < n2; i++) {
		t.b[i] = x4.b[i];
	}

	/* divide [t1, t2, t3] = [r1, r2, x4] by [y1, y2] */
	BIdivmod_divide_conquere_32(t, y, &q2, r_);
	BIfree(t);

	/* construct [q1, q2] */
	assert(q2.l == n2);
	*q_ = BInew(q1.l + n2);
	/* copy q2 */
	for (i = 0; i < q2.l; i++) {
		q_->b[i] = q2.b[i];
	}
	/* fill in beetwen zeros */
	for (; i < n2; i++) {
		q_->b[i] = 0;
	}
	/* copy q1 */
	for (i = 0; i < q1.l; i++) {
		q_->b[n2+i] = q1.b[i];
	}
	BIfree(q1);
	BIfree(q2);

	/* r_ is correct */
}

void BIdivmod_divide_conquere_32(BI x, BI y, BI *q_, BI *r_) {
	BI x1, x2, x3;
	BI y1, y2;
	BI x12;
	BI y1_shift, temp;
	BI qq, r1;
	BI d;
	BI r;
	int i;
	int n = y.l/2;

/*	printf("subdziele32 %d przez %d\n", x.l, y.l); */

	if (x.l < y.l) {
		*q_ = BIzero();
		*r_ = BIclone(x);
		return;
	}

	if (!(y.l & 1) || y.l < DIV_LIMIT) {
/*		printf("base case\n"); */
		BIdivmod_basecase(x, y, q_, r_);
		return;
	}

	assert(!is_zero(y));

	assert(x.l > y.l);
	assert(x.l <= y.l+n);
	assert(y.b[y.l-1] & (1 << (BIBit-1))); /* leading bit in y set */
	assert(x.l > y.l);
	assert(x.l < y.l+n);

	x3 = BIslice(x, 0, n);
	x2 = BIslice(x, n, n+n);
	x1 = BIslice(x, n+n, x.l);

	y2 = BIslice(y, 0, n);
	y1 = BIslice(y, n, y.l);

	x12 = BIslice(x, n, x.l);

	if (BIcmp(x1, y1) < 0) { /* x1 < y1 */
		BIdivmod_divide_conquere_42(x12, y1, &qq, &r1);
	} else { /* x1 >= y1 */

		/* r1 = [x1,x2] - [y1, 0] + [0, y1] */
		y1_shift = BIshift(y1, -n*BIBit);
		temp = BIsub(x12, y1_shift);
		BIfree(y1_shift);

		toadd(&temp, y1);

		/* qq = BI_beta^n - 1 */
		qq = BInew(n);
		for (i = 0; i < n; i++) {
			qq.b[i] = BIMask;
		}

		/** r1 == [x1, x2] - qq*y1 */
	}

	d = BImul(qq, y2);
	r = BIshift(r1, -n*BIBit);
	toadd(&r, x3); /* x3? */

	while (BIcmp(r, d) < 0) {  /* r < d,   r-d < 0 */
		toadd(&r, y);
		tosub(&qq, BI_1);
	}

	tosub(&r, d);

	BIfree(d);

	*r_ = r;
	*q_ = qq;
}

/* O(M(n) + n log n), dla mnozenia metoda szkolna, z M(n) = O(n^2) */
/* O(2 M(n) + n log n), dla mnozenia metoda Karatsuby, z M(n) = O(n^1.585) */
void BIdivmod_divide_conquere(BI x, BI y, BI *q_, BI *r_) {
	int m, j, n, shift, t, i;
	int k;
	BI z, r, q;


	if (is_zero(y)) {
		fprintf(stderr, "Division by zero\n");
		exit(151);
	}

	assert(!is_zero(y));

/*	printf("dziele %d przez %d\n", x.l, y.l); */

	m = 1;
	while (y.l > m*BI_DIV_LIMIT) {
		m *= 2;
	}
	m /= 2;
/*	printf("m = %d\n", m); */

	j = y.l/m;
	n = j*m;
	if (n < y.l) {
		j += 1;
		n += m;
	}
	assert(y.l <= n);
/*	printf("j = %d\n", j); */
/*	printf("n = %d\n", n); */

	shift = n*BIBit - 1 - BImsb(y);
/*	BIprint(y); */
/*	printf("shift = %d\n", shift); */
	if (shift) {
		x = BIshift(x, -shift);
		y = BIshift(y, -shift);
	}
/*	printf("dziele %d przez %d\n", x.l, y.l); */
/*	BIprint(y); */

	t = x.l/n;
/*	printf("t = %d\n", t); */
	while (t*n < x.l) {
		t++;
	}
	/* */
/*	printf("t = %d\n", t); */
	if (t < 2) {
		t = 2;
	}
/*	printf("t = %d\n", t); */

	*q_ = BInew((t-1)*n);

	z = BIslice(x, (t-2)*n, x.l);

	for (i = t-2; i >= 0; i--) {
		BIdivmod_divide_conquere_42(z, y, &q, &r);
		for (k = 0; k < q.l; k++) {
			q_->b[k + n*i] = q.b[k];
		}
		for (; k < n; k++) {
			q_->b[k + n*i] = 0;
		}
		BIfree(q);
		if (i < t-2) {
			BIfree(z);
		}
		if (i > 0) {
			z = BInew(n + r.l);
			for (k = 0; k < n; k++) {
				z.b[i] = x.b[i + (i-1)*n];
			}
			for (k = 0; k < r.l; k++) {
				z.b[k + n] = r.b[k];
			}
			BIfree(r);
		}
	}

	if (shift) {
		BIfree(x);
		BIfree(y);
		*r_ = BIshift(r, shift);
		BIfree(r);
	} else {
		*r_ = r;
	}

	*q_ = BIfix(*q_, q_->l - 1);
}

/* Newton's method for finding 1/y (with sufficient precision),
 * and then multiply it by x, x*1/y, 
 * and then fix it to make it correct.
 */
void BIdivmod_barret_newton(BI x, BI y, BI *q_, BI *r_) {
	BI v, si, ti, ui, wi, zi;
	int ki;
	int n;
	int free_qi_;
	BI mu;
	BI a, qi_, qi, bqi, ri, temp;

	assert(!(y.l == 0 && y.b[0] == 0));


	/* precision, for powmod usage it would be sufficient to have 2*y.l + 1 */
	n = 2*y.l + 1;
	assert(x.l <= n);
	ki = 1;
	/* newton method for calculating inversion */
	while (ki < n) {
		si = BIsqr(zi);
		ti = BIslice(v, v.l-ki, v.l);
		ui = BImul(ti, si);
		BIfree(si);
		wi = BIshift(zi, -1); /* BImul(BI_2, zi) */
		BIfree(zi);
		zi = BIsub(wi, ui);
		BIfree(wi);
		BIfree(ui);
	}

	mu = zi;

	assert(y.l <= x.l);
	assert(x.l <= 2*y.l);

	/* barret method for using inversion for divison */
	a = BIslice(y, x.l-1, y.l);
	qi_ = BImul(a, mu);
	BIfree(mu);
	qi = BIslice(qi_, (x.l-y.l+1), qi_.l);
	free_qi_ = 1;
	bqi = BImul(qi, y);
	if (BIcmp(x, bqi) >= 0) {
		ri = BIsub(x, bqi);
		BIfree(bqi);
		while (BIcmp(ri, y) >= 0) {
			tosub(&ri, y);
			temp = BIadd(qi, BI_1);
			if (free_qi_) {
				BIfree(qi_); /* qi is slice of qi_ */
			} else {
				BIfree(qi);
			}
			qi = temp;
			free_qi_ = 0;
		}
	} else {
		ri = BIsub(bqi, x); /* -ri */
		BIfree(bqi);
		assert(0);
	}
	*r_ = ri;
	*r_ = qi;
}


void BIdivmod(BI x, BI y, BI *q_, BI *r_) {
	assert(!is_zero(y));

	if (x.l < y.l) {
		*q_ = BIzero();
		*r_ = BIclone(x);
		return;
	}

		BIdivmod_basecase(x, y, q_, r_);
		return;


	if (x.l > BI_DIV_LIMIT && y.l > BI_DIV_LIMIT) {
		BI xprim, z;
		BI q2, r2;
		BIdivmod_divide_conquere(x, y, q_, r_);

/*
		BIdivmod_basecase(x, y, &q2, &r2);
		BIprint(x);
		BIprint(y);
		BIprint(*q_);
		BIprint(*r_);
		BIprint(q2);
		BIprint(r2);
*/
/*
		assert(BIeq(*r_, r2));
		assert(BIeq(*q_, q2));
		xprim = BImul(*q_, y);
		z = BIadd(xprim, *r_);
		assert(BIeq(z, x));
		BIfree(z);
		BIfree(xprim);
*/
	} else {
		BIdivmod_basecase(x, y, q_, r_);
	}
}

BI BImod(BI x, BI n) {
	BI q, r;
	assert(x.sign == 1);
	assert(n.sign == 1);
	BIdivmod(x, n, &q, &r);
	BIfree(q);
	return r;
}

BI BIdiv(BI x, BI n) {
	BI q, r;
	assert(x.sign == 1);
	assert(n.sign == 1);
	BIdivmod(x, n, &q, &r);
	BIfree(r);
	return q;
}

BI BImulmod(BI x, BI y, BI n) {
	BI temp, z;
	assert(x.sign == 1);
	assert(y.sign == 1);
	assert(n.sign == 1);
	temp = BImul(x, y);
	z = BImod(temp, n);
	BIfree(temp);
	return z;
}

BI BIsqrmod(BI x, BI n) {
	BI temp, z;
	assert(x.sign == 1);
	assert(n.sign == 1);
	temp = BIsqr(x);
	z = BImod(temp, n);
	BIfree(temp);
	return z;
}


BI BIpowmod(BI base, BI b, BI n) {
	BI r;
	BI temp;
	int i, j;
	BIType bit, bit0;
	int last_bit_pos;
	int can_free_base;
	int tested_bit;

	assert(base.sign == 1);
	assert(b.sign == 1);
	assert(n.sign == 1);

	r = BIone();

	can_free_base = 0;
	tested_bit = 0;

	bit0 = b.b[b.l-1];
	last_bit_pos = 0;
	while (bit0) {
		last_bit_pos++;
		bit0 >>= 1;
	}
	last_bit_pos += BIBit*(b.l-1);
	for (i = 0; i < b.l; i++) {
		bit0 = b.b[i];
		for (j = 0; j < BIBit && last_bit_pos > 0; j++, last_bit_pos--) {
			bit = (bit0 >> j) & 1;
			if (bit) {
				temp = BImulmod(r, base, n);
				BIfree(r);
				r = temp;
			}
			temp = BIsqrmod(base, n);
			if (can_free_base) { /* becuase this is calling argument */
				BIfree(base);
			}
			base = temp;

			can_free_base = 1;
			tested_bit++;
		}
	}
	if (can_free_base) {
		BIfree(base);
	}
	return r;
}


void BIprint(BI x) {
	int i = x.l-1;
	printf("0x");
	while (i) {
		printf("%04x", x.b[i]);
		i--;
	}
	printf("%04x", x.b[0]);
	printf("\n");
}

void BIprint10(BI x0) {
	char* s;
	int i, k;
	BI BI_10000 = BIvsmall(10000);
	BI x, q, r;

	s = malloc(sizeof(char)*(9 + 5*x0.l)); /* exactly 4.81648 decimal digits per 16-bits, 0.30103 dec digit per 1-bit */
	if (s == NULL) {
		perror("malloc");
		exit(6);
	}
	i = 8 + 5*x0.l;
	s[i] = '\0';
	i -= 4;
	BIdivmod(x0, BI_10000, &q, &r);
	assert(r.l == 1);
	k = r.b[0];
	s[i+3] = '0'+(k%10);
	k /= 10;
	s[i+2] = '0'+(k%10);
	k /= 10;
	s[i+1] = '0'+(k%10);
	k /= 10;
	s[i] = '0'+(k%10);
	BIfree(r);
	x = q;
	while (i > 2) {
		i -= 4;
		BIdivmod(x, BI_10000, &q, &r);
		assert(r.l == 1);
		k = r.b[0];
		s[i+3] = '0'+(k%10);
		k /= 10;
		s[i+2] = '0'+(k%10);
		k /= 10;
		s[i+1] = '0'+(k%10);
		k /= 10;
		s[i] = '0'+(k%10);
		BIfree(r);
		BIfree(x);
		x = q;
		if (q.l == 1 && q.b[0] == 0) {
			break;
		}
	}
	while (i < 8 + 5*x0.l) {
		if (s[i] == '0') {
			i++;
		} else {
			break;
		}
	}
	fprintf(stdout, "%s\n", s+i);
	BIfree(BI_10000);
	free(s);
	BIfree(x);
	return;
}


/* rabin miller */
#define BI_COMPOSITE        0
#define BI_PROBABLE_PRIME   1

/** pojedynczy test millera-rabina */
/* input: n;  nm1 = n-1;   n-1 = d*2^s;  random a: 2<=a<= n-1. */
int rabin_miller_pass(BI n, BI nm1, uint32_t s, BI d, BI a) {
	BI x;
	uint32_t r;

	x = BIpowmod(a, d, n);

	if (BIeq(x, BI_1)) {
		BIfree(x);
		return PROBABLE_PRIME;
	}

	if (BIeq(x, nm1)) {
		BIfree(x);
		return PROBABLE_PRIME;
	}

	if (s > 0) {
		for (r = 1; r <= s-1; r++) {
			BI temp = BIsqrmod(x, n);
			BIfree(x);
			x = temp;

			if (BIeq(x, BI_1)) {
				/* x==1, so all next BIsqrmod(x, n) will return 1 */
				BIfree(x);
				return COMPOSITE;
			}

			if (BIeq(x, nm1)) {
				BIfree(x);
				return PROBABLE_PRIME;
			}
		}
	}
	BIfree(x);

	return COMPOSITE;
}

/** zwraca 0 jesli liczba jest zlozona,
 * 1 jesli z duzym prawdopodobienstwem pierwsza
 */
int rabin_miller(BI n, int k) {
	int i, s;
	BI d, a, nm1;

#if 0
	if (n <= 1) {
		return COMPOSITE;
	}
	if (n == 2) {
		return PROBABLE_PRIME;
	}
#endif

	/* sprowadz n-1 do postaci d * 2^s, z nieparzystym d */
	s = 0;
	nm1 = BIsub(n, BI_1);
	d = BIclone(nm1);
	/* this can be done fast by first determining shift, and then performing single shift,
	   and if it is equal 0 modulo 16, then perform just pointer magic */
	while (!(d.b[0] & 1)) {
		BI temp = BIshift(d, 1);
		BIfree(d);
		d = temp;
		s++;
	}

#define RM_EARLY 1

#define RM_TEST(aa) do { \
	DD("testing early %d\n", aa.b[0]); \
	if ( (rabin_miller_pass(n, nm1, s, d, (aa)) == COMPOSITE)) { \
		BIfree(nm1); \
		BIfree(d); \
		return COMPOSITE; \
	} } while (0)
#define RM_TEST_G(aa) do { \
	if ((rabin_miller_pass(n, nm1, s, d, (aa)) == COMPOSITE)) { \
		return COMPOSITE; \
	} } while (0)

#if 0
	if (n < 341550071728321uLL && k > 8) {
		/* for "small" n we perform deterministic test (and do not use k at all) */
		/* for < 4759123141uLL, it is sufficient to test 2, 7, 61 */
		RM_TEST(2);
		RM_TEST(7);
		RM_TEST(61);
		/* additional tests for slightly bigger number */
		if (n >=  4759123141uLL) {
			RM_TEST(3);
			RM_TEST(5);
			RM_TEST(11);
			RM_TEST(13);
			RM_TEST(17);
		}
	} else {
#endif

/* first test by small numbers */

#ifdef RM_EARLY
		RM_TEST(BI_2);
		RM_TEST(BI_7);
		RM_TEST(BI_61);
		RM_TEST(BI_3);
		RM_TEST(BI_5);
		RM_TEST(BI_11);
		RM_TEST(BI_13);
		RM_TEST(BI_17);
		RM_TEST(BI_61);
#endif


		/* randomized , probabilistic test */
		for (i = 0; i < k; i++) {
			DD("testing %d\n", i);
			do {
				a = BIrng(nm1); /* wylosuj losowa liczbe */
				if ((a.l == 1) || BIcmp(nm1, a) < 0) {
					BIfree(a);
					DD("relosowanie\n");
					continue;
				} else {
					break;
				}
			} while (1);
			if (rabin_miller_pass(n, nm1, s, d, a) == COMPOSITE) {
				BIfree(a);
				BIfree(nm1);
				BIfree(d);
				return COMPOSITE;
			}
			BIfree(a);
		}
		BIfree(nm1);
		BIfree(d);
#if 0
	}
#endif

#undef RM_TEST
#undef RM_TEST_G
#undef RM_EARLY

	return PROBABLE_PRIME; /* prime if passed all tests */
}

int divisible_2(BI x) {
	return ((x.b[0] & 1) == 0);
}

/* divisibility tests */

/* http://www-graphics.stanford.edu/~seander/bithacks.html */
/* count the number of bits set in v */
int count_bits(uint32_t v) {
	uint32_t c; /* c accumulates the total bits set in v */
	for (c = 0; v; c++) {
		v &= (v - 1); /* clear the least significant bit set */
	}
	return c;
}

/* http://www.physicsforums.com/showthread.php?t=352014 */
int divisible_3(BI x) {
	int i;
	int32_t acc = 0;
	for (i = 0; i < x.l; i++) {
		acc += count_bits(x.b[i] & 0xaaaa);  /* 2#1010 1010 1010 1010 */
		acc -= count_bits(x.b[i] & 0x5555);  /* 2#0101 0101 0101 0101 */
	}
	return ((acc % 3) == 0);
}

/* http://mathforum.org/library/drmath/view/55908.html */
int divisible_5(BI x) {
	int i;
	uint32_t acc = 0;
	for (i = 0; i < x.l; i++) {
		acc += (x.b[i] & 0xff);
		acc += (x.b[i] >> 8);
	}
	return ((acc % 5) == 0);
}
int divisible(BI x, BI n) {
	BI r = BImod(x, n);
	int i = is_zero(x);
	BIfree(r);
	return i;
}

int divisible_7(BI x) {
	return divisible(x, BI_7);
}
int divisible_17(BI x) {
	return divisible(x, BI_17);
}
int divisible_31(BI x) {
	return divisible(x, BI_31);
}

/** zwraca 0 jesli liczb jest zlozona,
 * 1 jesli z bardzo duzym prawdopodobienstwem pierwsza */
int is_prime(BI n, int proposed_k) {
#if 0
	if (n <= 3) {
		if (n == 2 || n == 3) { return PROBABLE_PRIME; }
		return COMPOSITE;
	}
#endif
	/* zakladamy ze liczba jest duza, wieksza niz 100 */

	/* najpierw proste testy, aby niepotrzebnie nie wykonywac calego testu */


	if (divisible_2(n)) {
		return COMPOSITE;
	}

/* ponizsze male testy odfiltrowuja szybko okolo 64% liczb
 * dodanie 23 da 66%
 */
/* podzielnosc przez 3, 5, 7 mozna zrobic podobnymi regulami do 11, 9 w dziesietnym
 * generalnie patrzymy na krotkie cykle 2^k mod 5, 
 */

	/* -33.3% */
	if (divisible_3(n)) {
		return COMPOSITE;
	}
	/* -13..3% */
	if (divisible_5(n)) {
		return COMPOSITE;
	}
#if 0
	/* -7.6% */
	if (divisible_7(n)) {
		return COMPOSITE;
	}
	/* -4.1% */
/*
	if (divisible(n, BI_11)) {
		return COMPOSITE;
	}
*/
	if (divisible_17(n)) {
		return COMPOSITE;
	}
	if (divisible_31(n)) {
		return COMPOSITE;
	}
	if (divisible(n, BI_13)) {
		return COMPOSITE;
	}
	if (divisible(n, BI_23)) {
		return COMPOSITE;
	}
	if (divisible(n, BI_37)) {
		return COMPOSITE;
	}
#endif

	/* wykonaj test rabin_miller (30 tutaj by wystarczylo) */
	if (rabin_miller(n, BI_RABIN_MILLER_TESTS+proposed_k) == PROBABLE_PRIME) {
		return PROBABLE_PRIME;
	} else {
		return COMPOSITE;
	}
}

int BIinit() {
	BI_0 = BIzero();
	BI_1 = BIone();
	BI_2 = BIvsmall(2);
	BI_10000 = BIvsmall(100000);

	BI_beta = BIshift(BI_1, -BIBit);
}

int BIdeinit() {
	BIfree(BI_0);
	BIfree(BI_1);
	BIfree(BI_2);
	BIfree(BI_3);
	BIfree(BI_10000);
	BIfree(BI_beta);
}


#undef BI_RABIN_MILLER_TESTS

#undef _
#undef __
#undef DD
#undef max
#undef min

#undef BI_KARATSUBA_MUL_P
#undef BI_KARATSUBA_SQR_P
#undef BI_PRESHIFT

#undef BI_DIV_LIMIT

#ifndef BI_EXPOSE_INTERNALS
#undef BIsmallconst

#undef BIBit
#undef BIByte
#undef BIType
#undef BIMask
#endif

#endif /*  __BI_C_H__ */
