#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define ROTL(x, n)			(((x) << (n)) | ((x) >> (32 - (n))))
#define ENDIAN_CHANGE(X)	((ROTL((X),  8) & 0x00ff00ff) | (ROTL((X), 24) & 0xff00ff00))
#define WORD_LEN			8

typedef unsigned long long Long;
typedef unsigned int word;
typedef unsigned char byte;

typedef struct _ECC_PT_ {
	word X[WORD_LEN];
	word Y[WORD_LEN];
	word Z[WORD_LEN];

	char flag;
}ECC_PT;

//Function.c
__int64 cpucycles();
unsigned char getHex(unsigned char ch);
void convertStr2Byte(unsigned char* from, int size, unsigned char* to);

//ECADD
int addition(word* A, word* B, word* output);
int subtract(word* A, word* B, word* output);
int compare(word* A, word* B);
int addition_mod_P(word* A, word* B, word* P, word* output);
int subtract_mod_P(word* A, word* B, word* P, word* output);

//ECMUL
int OS_64bit_version(word* A, word* B, word* output);
int Fast_Reduction_NIST_P256(word* A, word* P, word* output);

//Inverse_function

int check_value_Binary(word* A);
void Right_Shift(word* A);
int Binary_Inverse(word* P, word* A, word* output);
int Fermat_based_inversion(word* z, word* P, word* output);

//Jacobian
int Jaco2Aff(
	word* X, word* Y, word* Z,
	word* P,
	word* output_X, word* output_Y
);

int ECDBL_Jacobian(
	word* X1, word* Y1, word* Z1, int flag,
	word* P,
	word* output_X, word* output_Y, word* output_Z
);

int ECADD_Jacobian(
	word* X1, word* Y1, word* Z1, int Pflag,
	word* X2, word* Y2, int Qflag,
	word* P,
	word* output_X, word* output_Y, word* output_Z
);


int LtoR_SM_Jaco_ver(word* scalar, word* X, word* Y, word* P, word* outputX, word* outputY, word* outputZ);

//Affine
int ECADD(word* X1, word* X2, word* Y1, word* Y2, word* P, word* outputX, word* outputY, int PFlag, int QFlag);
int ECDBL(word* X1, word* Y1, word* P, word* outputX, word* outputY, int Flag);
int LtoR_SM(word* scalar, word* X, word* Y, word* P, word* outputX, word* outputY);


void wNAF_Jaco(word* X, word* Y, word* Z, char* NAF, word* RX, word* RY, word* RZ, word* P);
void NAF_Gen(word* scalar, char* NAF, word* Prime);
void wNAF_Aff(word* X, word* Y, word* P, word* outputX, word* outputY, char* NAF);