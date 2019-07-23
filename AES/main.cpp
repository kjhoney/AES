#include <stdio.h>
#include <utility>
#define IP 0x0169					// Irreducible Polynomials  : x^8 + x^6 + x^5 + x^3 + x^1
#define CONSTANT 0x15
using namespace std;

int key_matrix[4][4][11];			// 11라운드의 key
int state[4][4];					// 매 라운드의 함수를 수행한 후 text 상태를 저장하는 행렬
int temp[4];
int RC[11] = {0, };					// RC

void RC_init() {					// RC 초기화 함수
	printf("RC: ");
	RC[1] = 1;
	for (int i = 2; i < 11; i++) {	// 1부터 left shift 하면서 배열에 저장
		RC[i] = (RC[i - 1] << 1);
		if (RC[i] > 0xff)			// 8비트를 넘어가면 mod M
			RC[i] ^= IP;
	}
	for (int i = 1; i < 11; i++)
		printf("%02X ", RC[i]);
	puts("");
}

void print_state() {				// 현재 4x4 text 행렬 상태 출력 함수
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			printf("%02X ", state[j][i]);
	}
	puts("");
}

void print_key(int round) {			// 각 라운드별 4x4 key 행렬 상태 출력 함수
	printf("ROUND %d: ", round);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			printf("%02X ", key_matrix[j][i][round]);
	}
	puts("");
}

int gf_add(int a, int b) {			// Galois Field 내 더하기 연산
	return a ^ b;
}

int gf_mul(int a, int b) {			// Galois Field 내 곱하기 연산
	int tmp = 0;
	int hi_bit_set;					// Carry 판별 비트
	for (int i = 0; i < 8; i++) {
		if (b & 1)					// b의 맨 아래 한 비트 검사
			tmp ^= a;				// 1일 때 a를 더해줌
		hi_bit_set = (a & 0x80);	// a의 맨 위 비트가 1이면 Carry 발생
		a <<= 1;					// a를 left shift
		if (hi_bit_set)				// Carry 발생시 mod M 연산
			a ^= IP;
		b >>= 1;					// b를 right shift
	}
	return tmp;
}

int gf_inv(int a, int M) {			// Galois Field 내 곱셈에 대한 역원을 찾는 함수
	if (!a) return 0;
	for (int i = 0; i < 256; i++) {	// 1부터 255 까지 a * x = 1 (mod M) 인 x를 찾으면 return
		if (gf_mul(a, i) == 1)
			return i;
	}
}

/*int Extended_Euclid(int m, int b) {
	int A[3] = { 1, 0, m };
	int B[3] = { 0, 1, b };
	while (1) {
		if (B[2] == 0)
			return 0;
		else if (B[2] == 1)
			return B[1];
		int Q = A[2] / B[2];
		printf("Q : %d\n", Q);
		int T[3] = { A[0] ^ gf_mul(Q,B[0]), A[1] ^ gf_mul(Q,B[1]), A[2] ^ gf_mul(Q ,B[2]) };
		for (int i = 0; i < 3; i++) {
			A[i] = B[i];
			B[i] = T[i];
			printf("%02X ", B[i]);
		}
		puts("");
	}
}*/

int SBOX(int input) {						// S BOX 구현 함수
	int tmp, inv, ret = 0;	
	int matrix[8][8] = { 1,0,0,0,1,1,1,1,	// Affine 변환 행렬
						1,1,0,0,0,1,1,1,
						1,1,1,0,0,0,1,1,
						1,1,1,1,0,0,0,1,
						1,1,1,1,1,0,0,0,
						0,1,1,1,1,1,0,0,
						0,0,1,1,1,1,1,0,
						0,0,0,1,1,1,1,1 };
	inv = gf_inv(input, IP);				// inv = input의 곱셈에 대한 역원
	for (int k = 0; k < 8; k++) {			// SBOX 행렬곱
		tmp = 0;
		for (int l = 0; l < 8; l++) {
			if (inv & (1 << l))
				tmp ^= matrix[k][l];
		}
		ret |= (tmp << k);					// 8비트 각 자리수를 8x1행렬의 각 자리라고 생각
	}
	ret ^= CONSTANT;						// 상수값 XOR연산
	return ret;								// 결과값 반환
}

int Inverse_SBOX(int input) {
	int tmp, tmp2, ret = 0;
	int newc = 0;
	int inv_matrix[8][8] = { 0,1,0,1,0,0,1,0,	// Affine 변환 행렬
							0,0,1,0,1,0,0,1,
							1,0,0,1,0,1,0,0,
							0,1,0,0,1,0,1,0,
							0,0,1,0,0,1,0,1,
							1,0,0,1,0,0,1,0,
							0,1,0,0,1,0,0,1,
							1,0,1,0,0,1,0,0 };
	tmp = input % 16;
	input = (tmp << 4) + (input / 16);
	for (int k = 0; k < 8; k++) {			// Inverse SBOX 행렬곱 ( 역행렬을 입력값과 상수에 모두 곱해줌 )
		tmp = 0; tmp2 = 0;
		for (int l = 0; l < 8; l++) {
			if (input & (1 << l))
				tmp ^= inv_matrix[k][l];
			if (CONSTANT & (1 << (7 - l)))
				tmp2 ^= inv_matrix[k][l];
		}
		ret |= (tmp << k);					// 8비트 각 자리수를 8x1행렬의 각 자리라고 생각
		newc |= (tmp2 << (7 - k));			// 새로 연산한 상수
	}
	ret ^= newc;							// 상수값 XOR 연산
	ret = gf_inv(ret, IP);					// 행렬연산을 끝낸 뒤에 역원을 찾아 저장
	return ret;								// 결과값 반환
}

int substitute_byte(int byte,int flag) {	// 8비트 변수를 받아 치환하는 함수   Encryption flag : 0    Decryption flag : 1
	int row, column;
	row = byte >> 4;						// SBOX 행은 변수의 앞 4비트
	column = byte % 16;						// SBOX 열은 변수의 뒤 4비트
	if (!flag)
		byte = SBOX(byte);					// Encryption시 SBOX 참조
	else 
		byte = Inverse_SBOX(byte);			// Decryption시 INVERSE SBOX 참조
	return byte;
}

void substitute_bytes(int flag) {			// Substitute byte 단계				 Encryption flag : 0    Decryption flag : 1
	for (int i = 0; i < 4; i++) {			// 4x4 행렬의 각 바이트를 위 함수를 이용해 치환
		for (int j = 0; j < 4; j++)
			state[i][j] = substitute_byte(state[i][j],flag);
	}
	printf("SB: ");							// 치환후 state 배열 출력
	print_state();
}

void R_function(int round) {				// Key Expansion 의 각 라운드 별 R Function
	temp[0] = key_matrix[1][3][round - 1];	// 전 Round의 마지막열 위로 shift
	temp[1] = key_matrix[2][3][round - 1];
	temp[2] = key_matrix[3][3][round - 1];
	temp[3] = key_matrix[0][3][round - 1];
	for (int i = 0; i < 4; i++)
		temp[i] = substitute_byte(temp[i], 0); // 전 Round의 마지막 열 치환
	temp[0] ^= RC[round];					// 맨 위 바이트를 해당 Round의 RC 값과 XOR 연산
}

void key_expansion() {						// Key Expansion
	print_key(0);							// 최초 Key 출력 ( ROUND 0 )
	for (int i = 1; i < 11; i++) {			
		R_function(i);						// 라운드 별 R Function 수행
		for (int j = 0; j < 4; j++) {		// 매 라운드 첫번째 열은 R Function 수행 값과 전 Round의 첫번째열 XOR 연산
			key_matrix[j][0][i] = temp[j] ^ key_matrix[j][0][i - 1];
			for (int k = 1; k < 4; k++)		// 나머지 열은 전 Round의 열과 해당 Round의 바로 앞 열을 XOR 연산
				key_matrix[j][k][i] = key_matrix[j][k - 1][i] ^ key_matrix[j][k][i - 1];

		}
		print_key(i);						// 해당 Round Key 출력
	}
}

void add_round_key(int round) {				// Add Round Key 단계 ( Encryption과 Decryption의 수행방법이 동일 )
	for (int i = 0; i < 4; i++) {			// 해당 Round의 Key 와 Text를 XOR 연산
		for (int j = 0; j < 4; j++) {
			state[i][j] = state[i][j] ^ key_matrix[i][j][round];
		}
	}
	printf("AR: ");
	print_state();							// Text 출력
}

void shift_row(int flag) {					// Shift Row 단계    Encryption flag : 0    Decryption flag : 1
	int tmp;
	if (!flag) {							// Encryption 시 각 행마다 left shift
		for (int i = 1; i < 4; i++) {
			for (int j = 1; j <= i; j++) {	// 2번째 행은 1번, 3번째 행은 2번, 4번째 행은 3번 shift
				tmp = state[i][0];
				state[i][0] = state[i][1];
				state[i][1] = state[i][2];
				state[i][2] = state[i][3];
				state[i][3] = tmp;
			}
		}
	}
	else {									// Decryption 시 각 행마다 right shift
		for (int i = 1; i < 4; i++) {
			for (int j = 1; j <= i; j++) {	// 2번째 행은 1번, 3번째 행은 2번, 4번째 행은 3번 shift
				tmp = state[i][0];
				state[i][0] = state[i][3];
				state[i][3] = state[i][2];
				state[i][2] = state[i][1];
				state[i][1] = tmp;
			}
		}
	}
	printf("SR: ");
	print_state();							// Text 출력
}

void mix_column(int flag) {					// Mix Column 단계    Encryption flag : 0    Decryption flag : 1
	int matrix[4][4] = {					// Encryption시 곱해줄 행렬
		0x02, 0x03, 0x01, 0x01,
		0x01, 0x02, 0x03, 0x01,
		0x01, 0x01, 0x02, 0x03,
		0x03, 0x01, 0x01, 0x02
	};
	int inv_matrix[4][4] = {				// Decryption시 곱해줄 행렬
		0x0E, 0x0B, 0x0D, 0x09,
		0x09, 0x0E, 0x0B, 0x0D,
		0x0D, 0x09, 0x0E, 0x0B,
		0x0B, 0x0D, 0x09, 0x0E,
	};
	for (int i = 0; i < 4; i++) {			// 행렬곱
		int temp[4] = { 0, };
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				if (!flag)
					temp[j] ^= gf_mul(matrix[j][k], state[k][i]);
				else
					temp[j] ^= gf_mul(inv_matrix[j][k], state[k][i]);
			}
		}
		state[0][i] = temp[0];
		state[1][i] = temp[1];
		state[2][i] = temp[2];
		state[3][i] = temp[3];
	}
	printf("MC: ");
	print_state();							// Text 출력
}

int main() {
	int read;
	FILE * plain;
	FILE * key;
	FILE * cipher;
	FILE * decrypt;
	plain = fopen("plain.bin", "rb");
	key = fopen("key.bin", "rb");			// plain.bin 과 key.bin 파일 읽기모드로 오픈
	cipher = fopen("cipher.bin", "wb");
	decrypt = fopen("decrypt.bin", "wb");

	RC_init();								// RC 초기화

	printf("PLAIN: ");
	for (int i = 0; i < 4; i++) {			// plain.bin 파일을 1바이트씩 읽어 state에 저장 및 출력
		for (int j = 0; j < 4; j++) {
			read = fgetc(plain);
			state[j][i] = read;
			printf("%02X ", read);
		}
	}

	printf("\nKEY: ");
	for (int i = 0; i < 4; i++) {			// key.bin 파일을 1바이트씩 읽어 key 행렬에 저장 및 출력
		for (int j = 0; j < 4; j++) {
			read = fgetc(key);
			key_matrix[j][i][0] = read;
			printf("%02X ", read);
		}
	}

	printf("\n\n<------ENCRYPTION------>\n\n"); // ENCRYPTION
	printf("KEY EXPANSION\n");

	key_expansion();						// Key expansion 진행

	for (int i = 0; i < 11; i++) {			// 각 Round 마다 해야되는 단계들 수행
		printf("\nROUND %d\n", i);
		if (i) {							// ROUND 0 에서는 Add round key만 수행
			substitute_bytes(0);
			shift_row(0);
			if (i != 10)					// ROUND 10 에서는 Mix column 단계를 수행하지 않음
				mix_column(0);
		}
		add_round_key(i);
	}

	printf("\nCIPHER: ");
	print_state();							// 최종 Ciphertext 출력

	for (int i = 0; i < 4; i++)				// cipher.bin 파일에 state 행렬을 1Byte씩 씀
		for (int j = 0; j < 4; j++)
			fputc(state[j][i], cipher);

	printf("\n\n<------DECRYPTION------>\n\n"); // DECRYPTION

	for (int i = 0; i < 11; i++) {			// 각 Round 마다 해야되는 단계를 수행
		printf("\nROUND %d\n", i);
		if (i) {							// ROUND 0 에서는 Add round key만 수행
			shift_row(1);
			substitute_bytes(1);
		}
		add_round_key(10 - i);				// Add round key는 Key를 ROUND 10부터 거꾸로 사용
		if (i && i != 10)					// ROUND 10 에서는 Mix column 단계를 수행하지 않음
			mix_column(1);
	}

	printf("\nDECRPYTED: ");
	print_state();							// 최종 Plaintext 출력

	for (int i = 0; i < 4; i++)				// decrypt.bin 파일에 state 행렬을 1Byte씩 씀
		for (int j = 0; j < 4; j++)
			fputc(state[j][i], decrypt);

	fclose(plain);
	fclose(key);
	fclose(cipher);
	fclose(decrypt);
}