#define _CRT_SECURE_NO_WARNINGS
#include "Type.h"
#define TESTNUM 10000
int main() {

	unsigned long long cycle = 0;
	unsigned long long cycle1 = 0;
	unsigned long long cycle2 = 0;

	FILE* P = NULL;

	//SM Affine 좌표계 결과 작성 파일
	FILE* ret_SM_AFF = NULL;
	ret_SM_AFF = fopen("ret_SM_Affine.txt", "w");
	assert(ret_SM_AFF != NULL);

	//wNAF Affine 좌표계 결과 작성 파일
	FILE* ret_wNAF_AFF = NULL;
	ret_wNAF_AFF = fopen("ret_wNAF_Affine.txt", "w");

	/////////////////////////////////////////////////////////////////

	//SM Jacobian 좌표계 결과 작성 파일
	FILE* ret_SM_Jaco = NULL;
	ret_SM_Jaco = fopen("ret_SM_Jacobian.txt", "w");
	assert(ret_SM_Jaco != NULL);

	//wNAF Jacobian 좌표계 결과 작성 파일
	FILE* ret_wNAF_Jaco = NULL;
	ret_wNAF_Jaco = fopen("ret_wNAF_Jacobian.txt", "w");

	//저장되어 있는 소수 읽기
	P = fopen("P256.txt", "r");
	assert(P != NULL);

	//저장되어 있는 Scalar 읽는 파일
	FILE* TV_scalar = fopen("TV_Scalar.txt", "r");
	assert(TV_scalar != NULL);

	//구현 정확성 검증
	FILE* SM_Answer = NULL;
	SM_Answer = fopen("TV_SM_P256.txt", "r");
	assert(SM_Answer != NULL);


	//ECC Base Point 정의
	word Base_Point_X[8] = {
	0x6b17d1f2,0xe12c4247,0xf8bce6e5,0x63a440f2,
	0x77037d81,0x2deb33a0,0xf4a13945,0xd898c296
	};

	word Base_Point_Y[8] = {
		0x4fe342e2,0xfe1a7f9b,0x8ee7eb4a,0x7c0f9e16,
		0x2bce3357,0x6b315ece,0xcbb64068,0x37bf51f5
	};
	word Base_Point_Z[WORD_LEN] = { 0,0,0,0,0,0,0,1 };

	byte buf[2000];
	byte value[32];

	word P256[8];
	word outputX[8];
	word outputY[8];
	word outputZ[8];
	word value_X[8];
	word value_Y[8];
	word value_Z[8];
	word bufferX[8];
	word bufferY[8];
	word Scalar[8];
	char NAF[257] = { 0, };
	int i = 0;
	int j = 0;


	//소수 읽기
	fgets((char*)buf, sizeof(buf), P);
	convertStr2Byte(buf, 64, value);
	for (i = 0; i < 8; i++)
		P256[i] = ENDIAN_CHANGE((*(word*)(value + (4 * i))));
	fclose(P);

	//결과 파일 작성 + 성능 측정
#if 1

	//Affine ECADD
	for (i = 0; i < TESTNUM; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			bufferX[j] = i ^ j;
			bufferY[j] = i + j;
		}
		cycle1 = cpucycles();
		ECADD(Base_Point_X, bufferX, Base_Point_Y, bufferY, P256, outputX, outputY, 0, 0);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);
	}
	printf("Affine ECADD RDTSC = %10lld\n", (cycle) / TESTNUM);
	cycle = 0;

	//Jacobian ECADD
	for (i = 0; i < TESTNUM; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			bufferX[j] = i ^ j;
			bufferY[j] = i + j;
		}
		cycle1 = cpucycles();
		ECADD_Jacobian(Base_Point_X, Base_Point_Y, Base_Point_Z, 0, bufferX, bufferY, 0, P256, outputX, outputY, outputZ);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);
	}
	printf("Jacobian ECADD RDTSC = %10lld\n", (cycle) / TESTNUM);
	cycle = 0;

	//Affine ECDBL
	for (i = 0; i < TESTNUM; i++)
	{
		cycle1 = cpucycles();
		ECDBL(Base_Point_X, Base_Point_Y, P256, outputX, outputY, 0);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);
	}
	printf("Affine ECDBL RDTSC = %10lld\n", (cycle) / TESTNUM);
	cycle = 0;

	//Jacobian ECDBL
	for (i = 0; i < TESTNUM; i++)
	{
		cycle1 = cpucycles();
		ECDBL_Jacobian(Base_Point_X, Base_Point_Y, Base_Point_Z, 0, P256, outputX, outputY, outputZ);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);
	}
	printf("Jacobian ECDBL RDTSC = %10lld\n", (cycle) / TESTNUM);
	cycle = 0;

	//Affine SM
	for (i = 0; i < TESTNUM; i++)
	{
		fgets(buf, sizeof(buf), TV_scalar);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			Scalar[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));
		
		fgets(buf, sizeof(buf), TV_scalar);
		
		cycle1 = cpucycles();
		LtoR_SM(Scalar, Base_Point_X, Base_Point_Y, P256, outputX, outputY);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);
		
		for (j = 0; j < 8; j++)
			fprintf(ret_SM_AFF, "%08X", outputX[j]);
		fprintf(ret_SM_AFF, "\n");
		
		for (j = 0; j < 8; j++)
			fprintf(ret_SM_AFF, "%08X", outputY[j]);
		fprintf(ret_SM_AFF, "\n\n");
	}
	printf("Affine SM RDTSC = %10lld\n", (cycle) / TESTNUM);
	fclose(ret_SM_AFF);
	cycle = 0;
	fseek(TV_scalar, 0, SEEK_SET);

	//Affine wNAF
	for (i = 0; i < TESTNUM; i++)
	{
		fgets(buf, sizeof(buf), TV_scalar);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			Scalar[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));

		fgets(buf, sizeof(buf), TV_scalar);

		cycle1 = cpucycles();
		NAF_Gen(Scalar, NAF, P256);
		wNAF_Aff(Base_Point_X, Base_Point_Y, P256, outputX, outputY, NAF);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);

		for (j = 0; j < 8; j++)
			fprintf(ret_wNAF_AFF, "%08X", outputX[j]);
		fprintf(ret_wNAF_AFF, "\n");

		for (j = 0; j < 8; j++)
			fprintf(ret_wNAF_AFF, "%08X", outputY[j]);
		fprintf(ret_wNAF_AFF, "\n\n");
		memset(NAF, 0, 257);
	}
	printf("Affine wNAF RDTSC = %10lld\n", (cycle) / TESTNUM);
	fclose(ret_wNAF_AFF);
	cycle = 0;
	fseek(TV_scalar, 0, SEEK_SET);

	//Jacobian SM
	for (i = 0; i < TESTNUM; i++) {
		
		fgets(buf, sizeof(buf), TV_scalar);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
				Scalar[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));
		
		fgets(buf, sizeof(buf), TV_scalar);
		
		cycle1 = cpucycles();
		LtoR_SM_Jaco_ver(Scalar, Base_Point_X, Base_Point_Y, P256, outputX, outputY, outputZ);
		Jaco2Aff(outputX, outputY, outputZ, P256, value_X, value_Y);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);
		
		for (j = 0; j < 8; j++)
			fprintf(ret_SM_Jaco, "%08X", value_X[j]);
		fprintf(ret_SM_Jaco, "\n");
		
		for (j = 0; j < 8; j++)
			fprintf(ret_SM_Jaco, "%08X", value_Y[j]);
		fprintf(ret_SM_Jaco, "\n\n");
		
	}
	printf("Jacobian SM RDTSC = %10lld\n", (cycle) / TESTNUM);
	fclose(ret_SM_Jaco);
	cycle = 0;
	fseek(TV_scalar, 0, SEEK_SET);

	//Jacobian wNAF
	for (i = 0; i < TESTNUM; i++) {

		fgets(buf, sizeof(buf), TV_scalar);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			Scalar[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));

		fgets(buf, sizeof(buf), TV_scalar);

		cycle1 = cpucycles();
		NAF_Gen(Scalar, NAF, P256);
		wNAF_Jaco(Base_Point_X, Base_Point_Y, Base_Point_Z, NAF, value_X, value_Y, value_Z, P256);
		Jaco2Aff(value_X, value_Y, value_Z, P256, outputX, outputY);
		cycle2 = cpucycles();
		cycle += (cycle2 - cycle1);

		for (j = 0; j < 8; j++)
			fprintf(ret_wNAF_Jaco, "%08X", outputX[j]);
		fprintf(ret_wNAF_Jaco, "\n");

		for (j = 0; j < 8; j++)
			fprintf(ret_wNAF_Jaco, "%08X", outputY[j]);
		fprintf(ret_wNAF_Jaco, "\n\n");

	}
	printf("Jacobian wNAF RDTSC = %10lld\n", (cycle) / TESTNUM);
	fclose(ret_wNAF_Jaco);
	cycle = 0;

#endif

	//구현 정확성 검증
#if 0
	int checkSum = 1;
	//Affine wNAF 구현 정확성 검증
	for (i = 0; i < TESTNUM; i++)
	{
		//Scalar Read
		fgets(buf, sizeof(buf), TV_scalar);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			Scalar[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));
		fgets(buf, sizeof(buf), TV_scalar);

		//X Read
		fgets(buf, sizeof(buf), SM_Answer);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			bufferX[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));
		
		//Y Read
		fgets(buf, sizeof(buf), SM_Answer);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			bufferY[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));
		fgets(buf, sizeof(buf), SM_Answer);

		//Operation
		NAF_Gen(Scalar, NAF, P256);
		wNAF_Aff(Base_Point_X, Base_Point_Y, P256, outputX, outputY, NAF);

		//Check Answer
		for (j = 0; j < 8; j++)
		{
			if ((outputX[j] != bufferX[j]) || (outputY[j] != bufferY[j]))
			{
				printf("%d-th Affine wNAF Answer ERROR\n", i);
				checkSum = -1;
				goto NEXT;
			}
		}
		memset(NAF, 0, 257);
	}
	NEXT:
	if (checkSum == -1)
		printf("Affine wNAF Answer ERROR\n");
	else
		printf("Affine wNAF Answer Correct!\n");

	checkSum = 0;
	//Affine wNAF 구현 정확성 검증
	for (i = 0; i < TESTNUM; i++)
	{
		//Scalar Read
		fgets(buf, sizeof(buf), TV_scalar);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			Scalar[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));
		fgets(buf, sizeof(buf), TV_scalar);

		//X Read
		fgets(buf, sizeof(buf), SM_Answer);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			bufferX[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));

		//Y Read
		fgets(buf, sizeof(buf), SM_Answer);
		convertStr2Byte(buf, 64, value);
		for (j = 0; j < 8; j++)
			bufferY[j] = ENDIAN_CHANGE((*(word*)(value + (4 * j))));
		fgets(buf, sizeof(buf), SM_Answer);

		//Operation
		NAF_Gen(Scalar, NAF, P256);
		wNAF_Jaco(Base_Point_X, Base_Point_Y, Base_Point_Z, NAF, value_X, value_Y, value_Z, P256);
		Jaco2Aff(value_X, value_Y, value_Z, P256, outputX, outputY);

		//Check Answer
		for (j = 0; j < 8; j++)
		{
			if ((outputX[j] != bufferX[j]) || (outputY[j] != bufferY[j]))
			{
				printf("%d-th Jacobian wNAF Answer ERROR\n", i);
				checkSum = -1;
				break;
				goto NEXT2;
			}
		}
		memset(NAF, 0, 257);
	}
NEXT2:
	if (checkSum == -1)
		printf("Jacobian wNAF Answer ERROR\n");
	else
		printf("Jacobian wNAF Answer Correct!\n");
#endif

	fclose(TV_scalar);
	fclose(SM_Answer);
}

