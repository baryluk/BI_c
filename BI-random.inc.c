/*** najpierw generator liczb losowych */

#define BI_WELL_W 32
#define BI_WELL_R 16
#define BI_WELL_P 0
#define BI_WELL_M1 13
#define BI_WELL_M2 9
#define BI_WELL_M3 5

#define BI_WELL_MAT0POS(t,v) (v^(v>>t))
#define BI_WELL_MAT0NEG(t,v) (v^(v<<(-(t))))
#define BI_WELL_MAT3NEG(t,v) (v<<(-(t)))
#define BI_WELL_MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))

#define BI_WELL_V0            BI_WELL_STATE[BI_WELL_state_i                   ]
#define BI_WELL_VM1           BI_WELL_STATE[(BI_WELL_state_i+BI_WELL_M1) & 0x0000000fU]
#define BI_WELL_VM2           BI_WELL_STATE[(BI_WELL_state_i+BI_WELL_M2) & 0x0000000fU]
#define BI_WELL_VM3           BI_WELL_STATE[(BI_WELL_state_i+BI_WELL_M3) & 0x0000000fU]
#define BI_WELL_VRm1          BI_WELL_STATE[(BI_WELL_state_i+15) & 0x0000000fU]
#define BI_WELL_VRm2          BI_WELL_STATE[(BI_WELL_state_i+14) & 0x0000000fU]
#define BI_WELL_newV0         BI_WELL_STATE[(BI_WELL_state_i+15) & 0x0000000fU]
#define BI_WELL_newV1         BI_WELL_STATE[BI_WELL_state_i                 ]
#define BI_WELL_newVRm1       BI_WELL_STATE[(BI_WELL_state_i+14) & 0x0000000fU]

static uint32_t BI_WELL_state_i = 0;
static uint32_t BI_WELL_STATE[R];
static uint32_t BI_WELL_z0, BI_WELL_z1, BI_WELL_z2;


void InitWELLRNG512a(uint32_t *init) {
   int j;
   BI_WELL_state_i = 0;
   for (j = 0; j < R; j++)
     BI_WELL_STATE[j] = init[j];
}

/* Copyright: Francois Panneton and Pierre L'Ecuyer (University of Montreal), Makoto Matsumoto (Hiroshima University) */
/* Notice: This code can be used freely for personal, academic, or non-commercial purposes. */
uint32_t WELLRNG512a() {
	BI_WELL_z0    = BI_WELL_VRm1;
	BI_WELL_z1    = BI_WELL_MAT0NEG (-16, BI_WELL_V0)    ^ BI_WELL_MAT0NEG (-15, BI_WELL_VM1);
	BI_WELL_z2    = BI_WELL_MAT0POS (11, BI_WELL_VM2)  ;
	BI_WELL_newV1 = BI_WELL_z1                  ^ BI_WELL_z2; 
	BI_WELL_newV0 = BI_WELL_MAT0NEG (-2, BI_WELL_z0)     ^ BI_WELL_MAT0NEG(-18,BI_WELL_z1)    ^ BI_WELL_MAT3NEG(-28,BI_WELL_z2) ^ BI_WELL_MAT4NEG(-5,0xda442d24U,BI_WELL_newV1) ;
	BI_WELL_state_i = (BI_WELL_state_i + 15) & 0x0000000fU;
	return BI_WELL_STATE[BI_WELL_state_i];
}

/* poprostu przypadkowe liczby z /dev/random do inita */
uint32_t BI_WELL_moj_init[16] = {
	0xf51999c1u, 0x2b3d3410u,  0x356136e2u, 0xb18f4af4u,
	0x8837a43fu, 0xae2bacd8u,  0x5db81d32u, 0x398f9be2u,
	0xd1ee8bd0u, 0x0d865d5cu,  0x9a4d5b56u, 0x02ddc322u,
	0x003db53fu, 0x03ed956du,  0x971f5bdeu, 0x7a4e4712u
};


uint64_t rng64() {
	uint64_t a, b;
	a = (uint64_t)WELLRNG512a();
	b = (uint64_t)WELLRNG512a();
	return (a<<32) | b;
}

/* zwraca liczbe losowba [a;b] */
int rngAB(uint64_t a, uint64_t b) {
	return a + rng64() % (b-a+1uL);
}

/* random, at most x */
/* O(n) */
BI BIrng(BI x) {
	int i;
	BI z = BInew(x.l);
	for (i = 0; i < x.l-1; i++) {
		z.b[i] = WELLRNG512a();
	}
	z.b[z.l-1] = WELLRNG512a() % (1+x.b[x.l-1]);
	return BIfix(z, z.l-1); /* just for sure */
}
