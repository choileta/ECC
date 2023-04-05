#include "Type.h"

int Jaco2Aff(
	word* X, word* Y, word* Z,
	word* P,
	word* output_X, word* output_Y
) 
{
	word mid_value[WORD_LEN];
	word MULT_value[WORD_LEN * 2];

	word Z_2[WORD_LEN];
	word Z_3[WORD_LEN];
	word Z_inv[WORD_LEN];
	word Z_value[WORD_LEN];

	//Z^2
	OS_64bit_version(Z, Z, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, Z_2);

	//Z^3
	OS_64bit_version(Z_2, Z, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, Z_3);
	
	//(Z^3)^-1
	Binary_Inverse(P, Z_3, Z_inv);

	//(Z^2)-1
	OS_64bit_version(Z_inv, Z, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, Z_value);

	OS_64bit_version(Y, Z_inv, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, output_Y);

	OS_64bit_version(X, Z_value, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, output_X);

	return 0;
}

int ECDBL_Jacobian(
	word* X1, word* Y1, word* Z1, int flag,
	word* P, 
	word* output_X, word* output_Y, word* output_Z
) 
{

	if (flag == 1) {
		return 1;
	}

	word T1[WORD_LEN];
	word T2[WORD_LEN];
	word T3[WORD_LEN];

	word X_3[WORD_LEN];
	word Y_3[WORD_LEN];
	word Z_3[WORD_LEN];

	word mid_value[WORD_LEN];
	word value[WORD_LEN];
	word MULT_value[WORD_LEN * 2];

	//T1 <- Z^2
	OS_64bit_version(Z1, Z1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T1);

	//T2 <- X1 - T1
	subtract_mod_P(X1, T1, P, T2);

	//T1 <- X1 + T1
	addition_mod_P(X1, T1, P, mid_value);
	memcpy(T1, mid_value, sizeof(word) * WORD_LEN);

	//T2 <- T2 * T1
	OS_64bit_version(T1, T2, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T2);

	//T2 <- 3T2
	addition_mod_P(T2, T2, P, mid_value);
	addition_mod_P(mid_value, T2, P, value);
	memcpy(T2, value, sizeof(word) * WORD_LEN);

	//Y3 <- 2Y1
	addition_mod_P(Y1, Y1, P, Y_3);

	//Z3 <- Y3 * Z1
	OS_64bit_version(Y_3, Z1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, Z_3);

	//Y3 <- Y3 * Y3
	OS_64bit_version(Y_3, Y_3, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, Y_3);

	//T3 <- Y3 * X1
	OS_64bit_version(Y_3, X1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T3);

	//Y3 <- Y3 * Y3
	OS_64bit_version(Y_3, Y_3, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, Y_3);

	//Y3 <- Y3/2
	char check_value = 0;
	check_value = (Y_3[7] & 0x1);
	if(check_value ==0)
		Right_Shift(Y_3);
	else {
		word carry = 0;
		carry = addition(Y_3, P, mid_value);
		memcpy(Y_3, mid_value, sizeof(word) * WORD_LEN);
		Right_Shift(Y_3);
		Y_3[0] = Y_3[0] + (carry << 31);
	}

	//X3 <- T2^2
	OS_64bit_version(T2, T2, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, X_3);

	//T1 <- 2T3
	addition_mod_P(T3, T3, P, T1);
	
	//X3 <- X3 - T1
	subtract_mod_P(X_3, T1, P, mid_value);
	memcpy(X_3, mid_value, sizeof(word) * WORD_LEN);

	//T1 <- T3 - X3
	subtract_mod_P(T3, X_3, P, T1);

	//T1 <- T1*T2
	OS_64bit_version(T1, T2, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T1);

	//Y3<- T1-Y3
	subtract_mod_P(T1, Y_3, P, mid_value);
	memcpy(Y_3, mid_value, sizeof(word) * WORD_LEN);

	//END
	memcpy(output_X, X_3, sizeof(word) * WORD_LEN);
	memcpy(output_Y, Y_3, sizeof(word) * WORD_LEN);
	memcpy(output_Z, Z_3, sizeof(word) * WORD_LEN);

	return 0;

}

int ECADD_Jacobian(
	word* X1, word* Y1, word* Z1, int Pflag,
	word* X2, word* Y2, int Qflag,
	word* P,
	word* output_X, word* output_Y, word* output_Z
) {

	word value[WORD_LEN] = { 0, }; 
	word mid_value[WORD_LEN];
	word MULT_value[WORD_LEN * 2];

	word T1[WORD_LEN];
	word T2[WORD_LEN];
	word T3[WORD_LEN];
	word T4[WORD_LEN];

	word X3[WORD_LEN];
	word Y3[WORD_LEN];
	word Z3[WORD_LEN];

	word Zvalue[WORD_LEN];

	int compare_value = 0;

	if (Qflag == 1) {
		memcpy(output_X, X1, sizeof(word) * WORD_LEN);
		memcpy(output_Y, Y1, sizeof(word) * WORD_LEN);
		memcpy(output_Z, Z1, sizeof(word) * WORD_LEN);
		return 0;
	}

	if (Pflag == 1) {
		memcpy(output_X, X2, sizeof(word) * WORD_LEN);
		memcpy(output_Y, Y2, sizeof(word) * WORD_LEN);
		memset(output_Z, 0, sizeof(word) * WORD_LEN); output_Z[7] = 1;
		return 0;
	}

	//T1 <- (Z1)^2
	OS_64bit_version(Z1, Z1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T1);

	//T2 <- T1 * Z1
	OS_64bit_version(T1, Z1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T2);

	//T1 <- T1 * X2
	OS_64bit_version(T1, X2, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T1);

	//T2 <- T2 * Y2
	OS_64bit_version(T2, Y2, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T2);

	//T1 <- T1 - X1
	subtract_mod_P(T1, X1, P, mid_value);
	memcpy(T1, mid_value, sizeof(word) * WORD_LEN);

	//T2 <- T2 - Y1
	subtract_mod_P(T2, Y1, P, mid_value);
	memcpy(T2, mid_value, sizeof(word) * WORD_LEN);

	//CHECK T1 == 0
	compare_value = compare(T1, P);
	if (compare_value == 0) {
		//CHECK T2 == 0
		if (compare(T2, P) == 0) {
			memcpy(mid_value, P, sizeof(word) * WORD_LEN);
			memset(value, 0, sizeof(word) * WORD_LEN); value[7] = 1;
			addition_mod_P(mid_value, value, P, Zvalue);
			ECDBL_Jacobian(X2, X2, Zvalue, 0, P, output_X, output_Y, output_Z);
			return 0;
		}
		else
			return 1;
	}

	//Z3 <- Z1 * T1
	OS_64bit_version(Z1, T1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, Z3);

	//T3 <- T1^2
	OS_64bit_version(T1, T1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T3);

	//T4 <- T3 * T1
	OS_64bit_version(T3, T1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T4);

	//T3 <- T3 * X1
	OS_64bit_version(T3, X1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T3);

	//T1 <- 2T3
	addition_mod_P(T3, T3, P, T1);

	//X3 <- T2^2
	OS_64bit_version(T2, T2, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, X3);

	//X3 <- X3 - T1
	subtract_mod_P(X3, T1, P, mid_value);
	memcpy(X3, mid_value, sizeof(word) * WORD_LEN);

	//X3 <- X3 - T4
	subtract_mod_P(X3, T4, P, mid_value);
	memcpy(X3, mid_value, sizeof(word) * WORD_LEN);

	//T3 <- T3 - X3
	subtract_mod_P(T3, X3, P, mid_value);
	memcpy(T3, mid_value, sizeof(word) * WORD_LEN);

	//T3 <- T3 * T2
	OS_64bit_version(T3, T2, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T3);

	//T4 <- T4 * Y1
	OS_64bit_version(T4, Y1, MULT_value);
	Fast_Reduction_NIST_P256(MULT_value, P, T4);

	//Y3 <- T3 - T4
	subtract_mod_P(T3, T4, P, Y3);

	//END
	memcpy(output_X, X3, sizeof(word)* WORD_LEN);
	memcpy(output_Y, Y3, sizeof(word)* WORD_LEN);
	memcpy(output_Z, Z3, sizeof(word)* WORD_LEN);
	return 0;
}

int LtoR_SM_Jaco_ver(word* scalar, word* X, word* Y, word* P, word* outputX, word* outputY, word* outputZ) {
	int i = 0;
	int j = 0;
	int Flag = 0;
	int check_value = 0;

	word Q_X[WORD_LEN];
	word Q_Y[WORD_LEN];
	word Q_Z[WORD_LEN];

	word mid_value_X[WORD_LEN];
	word mid_value_Y[WORD_LEN];
	word mid_value_Z[WORD_LEN];

	word Jaco_Z[WORD_LEN] = { 0,0,0,0,0,0,0,1 };
	word mid_value[WORD_LEN];
	word mid_value2[WORD_LEN];
	memcpy(Q_Z, Jaco_Z, sizeof(word) * WORD_LEN);
	memcpy(mid_value_Z, Jaco_Z, sizeof(word) * WORD_LEN);
	
	Flag = 1;
	for (i = 0; i < WORD_LEN; i++) {
		for (j = 0; j < 32; j++) {
			ECDBL_Jacobian(Q_X, Q_Y, Q_Z, Flag, P, mid_value_X, mid_value_Y, mid_value_Z);
			memcpy(Q_X, mid_value_X, WORD_LEN * sizeof(word));
			memcpy(Q_Y, mid_value_Y, WORD_LEN * sizeof(word));
			memcpy(Q_Z, mid_value_Z, sizeof(word) * WORD_LEN);

			check_value = (scalar[i] >> (31 - j)) & 0x1;
			
			if (check_value == 1) {
				
				ECADD_Jacobian(Q_X, Q_Y, Q_Z, Flag, X, Y, 0, P, mid_value_X, mid_value_Y, mid_value_Z);
				memcpy(Q_X, mid_value_X, WORD_LEN * sizeof(word));
				memcpy(Q_Y, mid_value_Y, WORD_LEN * sizeof(word));
				memcpy(Q_Z, mid_value_Z, WORD_LEN * sizeof(word));
				Flag = 0;
			}
		}
	}

	memcpy(outputX, Q_X, WORD_LEN * sizeof(word));
	memcpy(outputY, Q_Y, WORD_LEN * sizeof(word));
	memcpy(outputZ, Q_Z, WORD_LEN * sizeof(word));
	return 0;
}

void NAF_Gen(word* scalar, char* NAF, word* Prime)
{
	int i = 0;
	word One[8] = { 0, };
	word K[8] = { 0, };
	word value0[8] = { 0, };
	word value2[8] = { 0, };
	char mod = 0x00;
	memcpy(K, scalar, WORD_LEN * sizeof(word));
	One[7] = 0x1;

	while(compare(K, One) != -1)
	{
		if ((K[7] & 0x1) != 0) 
		{
			NAF[i] = K[7] & 0x0000000f;
			if (NAF[i] > 8)
			{
				NAF[i] -= 16;
				mod = ~(NAF[i]) + 1;
				value0[7] = mod;
				addition_mod_P(K, value0, Prime, value2);
				memcpy(K, value2, WORD_LEN * sizeof(word));
			}
			else
			{
				value0[7] = NAF[i];
				subtract_mod_P(K, value0, Prime, value2);
				memcpy(K, value2, WORD_LEN * sizeof(word));
			}

		}
		else
		{
			NAF[i] = 0;
		}
		Right_Shift(K);
		i += 1;
	}
}

void show(word* X) {
	for (int i = 0; i < 8; i++)
		printf("%08X ", X[i]);
	printf("\n");
}

void wNAF_Aff(word* X, word* Y, word* P, word* outputX, word* outputY, char* NAF)
{
	int i = 0; 
	int j = 0;
	int bits = 0;
	char temp = 0;

	//temp0
	word temp0X[WORD_LEN];
	word temp0Y[WORD_LEN];

	//temp1
	word temp1X[WORD_LEN];
	word temp1Y[WORD_LEN];

	//temp2
	word temp2X[WORD_LEN];
	word temp2Y[WORD_LEN];

	//buffer
	word temp3X[WORD_LEN];
	word temp3Y[WORD_LEN];

	//Pi Define
	word P0X_1[WORD_LEN];
	word P0Y_1[WORD_LEN];

	word P0X_3[WORD_LEN];
	word P0Y_3[WORD_LEN];

	word P0X_5[WORD_LEN];
	word P0Y_5[WORD_LEN];

	word P0X_7[WORD_LEN];
	word P0Y_7[WORD_LEN];

	word P1X_1[WORD_LEN];
	word P1Y_1[WORD_LEN];

	word P1X_3[WORD_LEN];
	word P1Y_3[WORD_LEN];

	word P1X_5[WORD_LEN];
	word P1Y_5[WORD_LEN];

	word P1X_7[WORD_LEN];
	word P1Y_7[WORD_LEN];


	//1
	memcpy(P0X_1, X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_1, Y, sizeof(word) * WORD_LEN);

	//3
	ECDBL(X, Y, P, temp0X, temp0Y, 0);
	ECADD(X, temp0X, Y, temp0Y, P, temp1X, temp1Y, 0, 0);
	memcpy(P0X_3, temp1X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_3, temp1Y, sizeof(word) * WORD_LEN);

	//5
	ECADD(P0X_3, temp0X, P0Y_3, temp0Y, P, temp2X, temp2Y, 0, 0);
	memcpy(P0X_5, temp2X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_5, temp2Y, sizeof(word) * WORD_LEN);

	//7
	ECADD(P0X_5, temp0X, P0Y_5, temp0Y, P, temp1X, temp1Y, 0, 0);
	memcpy(P0X_7, temp1X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_7, temp1Y, sizeof(word) * WORD_LEN);

	////////////////////////////////////////////////////////////////////////
	//1
	memcpy(P1X_1, P0X_1, sizeof(word) * WORD_LEN);
	memcpy(P1Y_1, P0Y_1, sizeof(word) * WORD_LEN);
	subtract_mod_P(P, P0Y_1, P, P1Y_1);

	//3
	memcpy(P1X_3, P0X_3, sizeof(word) * WORD_LEN);
	memcpy(P1Y_3, P0Y_3, sizeof(word) * WORD_LEN);
	subtract_mod_P(P, P0Y_3, P, P1Y_3);

	//5
	memcpy(P1X_5, P0X_5, sizeof(word) * WORD_LEN);
	memcpy(P1Y_5, P0Y_5, sizeof(word) * WORD_LEN);
	subtract_mod_P(P, P0Y_5, P, P1Y_5);

	//7
	memcpy(P1X_7, P0X_7, sizeof(word) * WORD_LEN);
	memcpy(P1Y_7, P0Y_7, sizeof(word) * WORD_LEN);
	subtract_mod_P(P, P0Y_7, P, P1Y_7);

	bits = 0;
	for (i = 256; i >= 0; i--)
	{

		if (bits == 1) {
			ECDBL(temp0X, temp0Y, P, temp1X, temp1Y, 0);
			if (NAF[i] != 0)
			{
				if (NAF[i] > 0)
				{
					switch ((NAF[i]-1)/2)
					{
					case 0:
						ECADD(temp1X, P0X_1, temp1Y, P0Y_1, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					case 1:
						ECADD(temp1X, P0X_3, temp1Y, P0Y_3, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					case 2:
						ECADD(temp1X, P0X_5, temp1Y, P0Y_5, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					case 3:
						ECADD(temp1X, P0X_7, temp1Y, P0Y_7, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					default:
						break;
					}
					continue;
				}
				else
				{
					switch ((-NAF[i]-1)/2)
					{
					case 0:
						ECADD(temp1X, P1X_1, temp1Y, P1Y_1, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					case 1:
						ECADD(temp1X, P1X_3, temp1Y, P1Y_3, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					case 2:
						ECADD(temp1X, P1X_5, temp1Y, P1Y_5, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					case 3:
						ECADD(temp1X, P1X_7, temp1Y, P1Y_7, P, temp2X, temp2Y, 0, 0);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						break;
					default:
						break;
					}
					continue;
				}
			}
			memcpy(temp0X, temp1X, sizeof(word) * WORD_LEN);
			memcpy(temp0Y, temp1Y, sizeof(word) * WORD_LEN);
		}
		if (NAF[i] == 0)
			continue;
		if (NAF[i] != 0)
		{
			switch ((NAF[i]-1)/2)
			{
			case 0:
				memcpy(temp0X, P0X_1, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_1, sizeof(word) * WORD_LEN);
				break;
			case 1:
				memcpy(temp0X, P0X_3, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_3, sizeof(word) * WORD_LEN);
				break;
			case 2:
				memcpy(temp0X, P0X_5, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_5, sizeof(word) * WORD_LEN);
				break;
			case 3:
				memcpy(temp0X, P0X_7, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_7, sizeof(word) * WORD_LEN);
				break;
			default:
				break;
			}
			bits = 1;
		}
	}
	memcpy(outputX, temp0X, sizeof(word) * WORD_LEN);
	memcpy(outputY, temp0Y, sizeof(word) * WORD_LEN);


}

void wNAF_Jaco(word* X, word* Y, word* Z, char* NAF, word* RX, word* RY, word* RZ, word* P)
{
	int i = 0;
	int j = 0;	
	int value = 0;
	char temp = 0;
	int bits = 0;
	//temp0
	word temp0X[WORD_LEN];
	word temp0Y[WORD_LEN];
	word temp0Z[WORD_LEN];
	//temp1
	word temp1X[WORD_LEN];
	word temp1Y[WORD_LEN];
	word temp1Z[WORD_LEN];
	//temp2
	word temp2X[WORD_LEN];
	word temp2Y[WORD_LEN];
	word temp2Z[WORD_LEN];
	//temp3
	word temp3X[WORD_LEN];
	word temp3Y[WORD_LEN];


	//Pi Define
	word P0X_1[WORD_LEN];
	word P0Y_1[WORD_LEN];
	word P0Z_1[WORD_LEN];

	word P0X_3[WORD_LEN];
	word P0Y_3[WORD_LEN];
	word P0Z_3[WORD_LEN];

	word P0X_5[WORD_LEN];
	word P0Y_5[WORD_LEN];
	word P0Z_5[WORD_LEN];

	word P0X_7[WORD_LEN];
	word P0Y_7[WORD_LEN];
	word P0Z_7[WORD_LEN];

	word P1X_1[WORD_LEN];
	word P1Y_1[WORD_LEN];
	word P1Z_1[WORD_LEN];

	word P1X_3[WORD_LEN];
	word P1Y_3[WORD_LEN];
	word P1Z_3[WORD_LEN];

	word P1X_5[WORD_LEN];
	word P1Y_5[WORD_LEN];
	word P1Z_5[WORD_LEN];

	word P1X_7[WORD_LEN];
	word P1Y_7[WORD_LEN];
	word P1Z_7[WORD_LEN];

	//1
	memcpy(P0X_1, X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_1, Y, sizeof(word) * WORD_LEN);
	memcpy(P0Z_1, Z, sizeof(word) * WORD_LEN);
	ECDBL_Jacobian(X, Y, Z, 0, P, temp0X, temp0Y, temp0Z);
	Jaco2Aff(temp0X, temp0Y, temp0Z, P, temp3X, temp3Y);
	ECADD_Jacobian(X, Y, Z, 0, temp3X, temp3Y, 0, P, temp1X, temp1Y, temp1Z);
	
	//3
	memcpy(P0X_3, temp1X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_3, temp1Y, sizeof(word) * WORD_LEN);
	memcpy(P0Z_3, temp1Z, sizeof(word) * WORD_LEN);

	//5
	ECADD_Jacobian(P0X_3, P0Y_3, P0Z_3, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
	memcpy(P0X_5, temp2X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_5, temp2Y, sizeof(word) * WORD_LEN);
	memcpy(P0Z_5, temp2Z, sizeof(word) * WORD_LEN);

	//7
	ECADD_Jacobian(P0X_5, P0Y_5, P0Z_5, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
	memcpy(P0X_7, temp2X, sizeof(word) * WORD_LEN);
	memcpy(P0Y_7, temp2Y, sizeof(word) * WORD_LEN);
	memcpy(P0Z_7, temp2Z, sizeof(word) * WORD_LEN);


	////////////////////////////////////////////////////////////////////////

	//1
	memcpy(P1X_1, P0X_1, sizeof(word) * WORD_LEN);
	memcpy(P1Y_1, P0Y_1, sizeof(word) * WORD_LEN);
	memcpy(P1Z_1, P0Z_1, sizeof(word) * WORD_LEN);
	subtract_mod_P(P, P0Y_1, P, P1Y_1);

	//3
	memcpy(P1X_3, P0X_3, sizeof(word) * WORD_LEN);
	memcpy(P1Y_3, P0Y_3, sizeof(word) * WORD_LEN);
	memcpy(P1Z_3, P0Z_3, sizeof(word) * WORD_LEN);
	subtract_mod_P(P, P0Y_3, P, P1Y_3);

	//5
	memcpy(P1X_5, P0X_5, sizeof(word) * WORD_LEN);
	memcpy(P1Y_5, P0Y_5, sizeof(word) * WORD_LEN);
	memcpy(P1Z_5, P0Z_5, sizeof(word) * WORD_LEN);
	subtract_mod_P(P, P0Y_5, P, P1Y_5);

	//7
	memcpy(P1X_7, P0X_7, sizeof(word)* WORD_LEN);
	memcpy(P1Y_7, P0Y_7, sizeof(word)* WORD_LEN);
	memcpy(P1Z_7, P0Z_7, sizeof(word)* WORD_LEN);
	subtract_mod_P(P, P0Y_7, P, P1Y_7);

	bits = 0;
	for (i = 256; i >=0; i--)
	{
		if (bits == 1) 
		{
			ECDBL_Jacobian(temp0X, temp0Y, temp0Z, 0, P, temp1X, temp1Y, temp1Z);
			if (NAF[i] != 0) {
				if (NAF[i] > 0) 
				{
					switch ((NAF[i]-1)/2)
					{
					case 0:
						Jaco2Aff(P0X_1, P0Y_1, P0Z_1, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					case 1:
						Jaco2Aff(P0X_3, P0Y_3, P0Z_3, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					case 2:
						Jaco2Aff(P0X_5, P0Y_5, P0Z_5, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					case 3:
						Jaco2Aff(P0X_7, P0Y_7, P0Z_7, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					default:
						break;
					}
					continue;
				}
				else
				{
					switch ((-NAF[i]-1)/2)
					{
					case 0:
						Jaco2Aff(P1X_1, P1Y_1, P1Z_1, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					case 1:
						Jaco2Aff(P1X_3, P1Y_3, P1Z_3, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					case 2:
						Jaco2Aff(P1X_5, P1Y_5, P1Z_5, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					case 3:
						Jaco2Aff(P1X_7, P1Y_7, P1Z_7, P, temp3X, temp3Y);
						ECADD_Jacobian(temp1X, temp1Y, temp1Z, 0, temp3X, temp3Y, 0, P, temp2X, temp2Y, temp2Z);
						memcpy(temp0X, temp2X, sizeof(word) * WORD_LEN);
						memcpy(temp0Y, temp2Y, sizeof(word) * WORD_LEN);
						memcpy(temp0Z, temp2Z, sizeof(word) * WORD_LEN);
						break;
					default:
						break;
					}
				}
				continue;
			}
			memcpy(temp0X, temp1X, sizeof(word) * WORD_LEN);
			memcpy(temp0Y, temp1Y, sizeof(word) * WORD_LEN);
			memcpy(temp0Z, temp1Z, sizeof(word) * WORD_LEN);
		}
		if (NAF[i] == 0)
			continue;
		if (NAF[i] != 0) {
			switch ((NAF[i] - 1) / 2)
			{
			case 0:
				memcpy(temp0X, P0X_1, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_1, sizeof(word) * WORD_LEN);
				memcpy(temp0Z, P0Z_1, sizeof(word) * WORD_LEN);
				break;
			case 1:
				memcpy(temp0X, P0X_3, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_3, sizeof(word) * WORD_LEN);
				memcpy(temp0Z, P0Z_3, sizeof(word) * WORD_LEN);
				break;
			case 2:
				memcpy(temp0X, P0X_5, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_5, sizeof(word) * WORD_LEN);
				memcpy(temp0Z, P0Z_5, sizeof(word) * WORD_LEN);
				break;
			case 3:
				memcpy(temp0X, P0X_7, sizeof(word) * WORD_LEN);
				memcpy(temp0Y, P0Y_7, sizeof(word) * WORD_LEN);
				memcpy(temp0Z, P0Z_7, sizeof(word) * WORD_LEN);
				break;
			default:
				break;
			}
			bits = 1;
		}
	}
	memcpy(RX, temp0X, sizeof(word) * WORD_LEN);
	memcpy(RY, temp0Y, sizeof(word) * WORD_LEN);
	memcpy(RZ, temp0Z, sizeof(word) * WORD_LEN);
}
