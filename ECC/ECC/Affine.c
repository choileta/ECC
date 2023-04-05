#include "Type.h"

/**
* @brief		ECC에서의 덧셈 함수
* @detail		(X3, Y3) = (X1, Y1) + (X2, Y2) in ECC
* @param		word * X1			값을 더할 첫 번째 X좌표 데이터
* @param		word * X2			값을 더할 두 번째 X좌표 데이터
* @param		word * Y1			값을 더할 첫 번째 Y좌표 데이터
* @param		word * Y2			값을 더할 두 번째 Y좌표 데이터
* @param		word * P			유한체 소수 P
* @param		word * outputX		X값을 저장할 배열
* @param		word * outputY		Y값을 저장할 배열
* @param		int PFlag			(X1, Y1) 무한원점 Flag
* @param		int QFlag			(X2, Y2) 무한원점 Flag
* @return		int					함수 성공 여부 or Flag 값
*/
int ECADD(word* X1, word* X2, word* Y1, word* Y2, word* P, word* outputX, word* outputY, int PFlag, int QFlag) {

	word mid_value[WORD_LEN];
	word mid_output[WORD_LEN * 2];
	word X2_X1[WORD_LEN];
	word X2_X1_INV[WORD_LEN];
	word Y2_Y1[WORD_LEN];
	word slope[WORD_LEN];
	word slope_SQR[WORD_LEN];

	if (PFlag == 1) {
		memcpy(outputX, X2, sizeof(word) * WORD_LEN);
		memcpy(outputY, Y2, sizeof(word) * WORD_LEN);
		return 0;
	}

	if (QFlag == 1) {
		memcpy(outputX, X1, sizeof(word) * WORD_LEN);
		memcpy(outputY, Y1, sizeof(word) * WORD_LEN);
		return 0;
	}

	//X3
	//X2-X1 & X2-X1 INVERSE
	subtract_mod_P(X2, X1, P, X2_X1);
	Fermat_based_inversion(X2_X1, P, X2_X1_INV);

	//Y2-Y1
	subtract_mod_P(Y2, Y1, P, Y2_Y1);

	//slope
	OS_64bit_version(Y2_Y1, X2_X1_INV, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, slope);

	//slope_SQR
	OS_64bit_version(slope, slope, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, slope_SQR);

	//Slope_SQR -X1
	subtract_mod_P(slope_SQR, X1, P, mid_value);

	//Last
	subtract_mod_P(mid_value, X2, P, outputX);

	//Y3
	//X1-X3
	subtract_mod_P(X1, outputX, P, mid_value);

	//slope * (x1-x3)
	OS_64bit_version(slope, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	//Last
	subtract_mod_P(mid_value, Y1, P, outputY);
	return 0;
}


/**
* @brief		ECC에서의 더블링 함수
* @detail		(X3, Y3) = 2(X1,Y1) in ECC
* @param		word * X1			값을 더블링 할 X좌표 데이터
* @param		word * Y1			값을 더블링 할 Y좌표 데이터
* @param		word * P			유한체 소수 P
* @param		word * outputX		X값을 저장할 배열
* @param		word * outputY		Y값을 저장할 배열
* @param		int Flag			(X1, Y1) 무한원점 Flag
* @return		int					함수 성공 여부 or Flag 값
*/
int ECDBL(word* X1, word* Y1, word* P, word* outputX, word* outputY, int Flag) {

	word mid_value[WORD_LEN];
	word mid_output[WORD_LEN * 2];
	word X_3[WORD_LEN];
	word X_2[WORD_LEN];
	word X_3_A[WORD_LEN];
	word slope[WORD_LEN];
	word INV_Y[WORD_LEN];
	word P_3[WORD_LEN] = { 0xFFFFFFFF, 0x00000001, 0x00000000, 0x00000000,
							0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFC };

	if (Flag == 1) {
		return 1;
	}

	//X3
	OS_64bit_version(X1, X1, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	addition_mod_P(mid_value, mid_value, P, X_2);
	addition_mod_P(mid_value, X_2, P, X_3);
	addition_mod_P(X_3, P_3, P, X_3_A);

	addition_mod_P(Y1, Y1, P, mid_value);
	Fermat_based_inversion(mid_value, P, INV_Y);

	OS_64bit_version(X_3_A, INV_Y, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, slope);

	OS_64bit_version(slope, slope, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	addition_mod_P(X1, X1, P, X_2);
	subtract_mod_P(mid_value, X_2, P, outputX);

	//Y3
	subtract_mod_P(X1, outputX, P, mid_value);

	OS_64bit_version(slope, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	subtract_mod_P(mid_value, Y1, P, outputY);

	return 0;
}


/**
* @brief		ECC에서의 스칼라 곱 LtoR 버전
* @detail		(X3, Y3) = K(X1,Y1) in ECC[K의 최상위 비트부터]
* @param		word * scalar		곱할 스칼라
* @param		word * X			Base-Point X
* @param		word * Y			Base-Point Y
* @param		word * P			유한체 소수 P
* @param		word * outputX		X값을 저장할 배열
* @param		word * outputY		Y값을 저장할 배열
* @return		int					함수 성공 여부
*/
int LtoR_SM(word* scalar, word* X, word* Y, word* P, word* outputX, word* outputY) {

	int i = 0;
	int j = 0;
	int Flag = 0;
	int check_value = 0;

	word Q_X[WORD_LEN] = { 0x00, };
	word Q_Y[WORD_LEN] = { 0x00, };

	word mid_value_X[WORD_LEN];
	word mid_value_Y[WORD_LEN];

	Flag = 1;
	for (i = 0; i < WORD_LEN; i++) {
		for (j = 0; j < 32; j++) {
			ECDBL(Q_X, Q_Y, P, mid_value_X, mid_value_Y, Flag);
			memcpy(Q_X, mid_value_X, WORD_LEN * sizeof(word));
			memcpy(Q_Y, mid_value_Y, WORD_LEN * sizeof(word));
			check_value = (scalar[i] >> (31 - j)) & 0x1;
			if (check_value == 1) {
				ECADD(Q_X, X, Q_Y, Y, P, mid_value_X, mid_value_Y, Flag, 0);
				memcpy(Q_X, mid_value_X, WORD_LEN * sizeof(word));
				memcpy(Q_Y, mid_value_Y, WORD_LEN * sizeof(word));
				Flag = 0;
			}
		}
	}

	memcpy(outputX, Q_X, WORD_LEN * sizeof(word));
	memcpy(outputY, Q_Y, WORD_LEN * sizeof(word));
}


/**
* @brief		ECC에서의 스칼라 곱 RtoL 버전
* @detail		(X3, Y3) = K(X1,Y1) in ECC[K의 최하위 비트부터]
* @param		word * scalar		곱할 스칼라
* @param		word * X			Base-Point X
* @param		word * Y			Base-Point Y
* @param		word * P			유한체 소수 P
* @param		word * outputX		X값을 저장할 배열
* @param		word * outputY		Y값을 저장할 배열
* @return		int					함수 성공 여부
*/
int RtoL_SM(word* scalar, word* X, word* Y, word* P, word* outputX, word* outputY) {

	word Q_X[WORD_LEN];
	word Q_Y[WORD_LEN];
	word P_value_X[WORD_LEN];
	word P_value_Y[WORD_LEN];
	word mid_value_X[WORD_LEN];
	word mid_value_Y[WORD_LEN];

	memcpy(P_value_X, X, sizeof(word) * WORD_LEN);
	memcpy(P_value_Y, Y, sizeof(word) * WORD_LEN);

	int i, j = 0;

	int check_value = 0;
	int Flag = 0;
	Flag = 1;
	for (i = 7; i >= 0; i--) {
		for (j = 0; j < 32; j++) {
			check_value = (scalar[i] >> (j)) & 0x1;

			if (check_value == 1) {
				ECADD(Q_X, P_value_X, Q_Y, P_value_Y, P, mid_value_X, mid_value_Y, Flag, 0);
				memcpy(Q_X, mid_value_X, sizeof(word) * WORD_LEN);
				memcpy(Q_Y, mid_value_Y, sizeof(word) * WORD_LEN);
				Flag = 0;
			}
			ECDBL(P_value_X, P_value_Y, P, mid_value_X, mid_value_Y, 0);
			memcpy(P_value_X, mid_value_X, sizeof(word) * WORD_LEN);
			memcpy(P_value_Y, mid_value_Y, sizeof(word) * WORD_LEN);
		}
	}

	memcpy(outputX, Q_X, WORD_LEN * sizeof(word));
	memcpy(outputY, Q_Y, WORD_LEN * sizeof(word));

}
