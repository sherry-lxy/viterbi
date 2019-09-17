#define _CRT_SECURE_NO_WARNINGS // visual studio scanf error

// C
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> // 処理時間はかる
#include <malloc.h>

// C++
#include<iostream> 
#include<cstdio>
#include<iomanip>
#include <ios>     // std::left, std::right

// 読みやすいため
#define and && 
#define or ||

// 数値定義
#define N 100 // 情報長さ
#define size 2 // 桁数
#define pi 3.14159265358979 // π
#define k 100000 // 回数
#define NIL 1000000

using namespace std;

typedef struct trellis {
	int hamming; // ハミング距離
	int before; // 前のステージ
	int now; // 現在のステージ
} trellis;

void Convolution(int input, int state, int *output, int *next); // 畳み込み符号状態遷移図(入力，ステージ，出力，次のステージ)
int MetricDistance(int x, int *y0, int *y1, int s); // 軟判定から求める距離(出力信号，0の枝の軟判定値，1の枝の軟判定値，桁数)
void Receive(int send[N], int receive[N], double sigma, int soft0[N][size], int soft1[N][size]); // 受信信号の取得(送信信号，受信信号，σ，軟判定の0の枝，軟判定の1の枝)
void Receive8(int send[N], int receive[N], double sigma, int soft0[N][size], int soft1[N][size]); // 受信信号の取得(送信信号，受信信号，σ，軟判定の0の枝，軟判定の1の枝)
void Receive16(int send[N], int receive[N], double sigma, int soft0[N][size], int soft1[N][size]);// 16値軟判定 受信信号の取得(送信信号，受信信号，σ，軟判定の0の枝，軟判定の1の枝)
void Trellis(trellis trellis[N][8], int receive[N], int soft0[N][size], int soft1[N][size]); // トレリス図の作成(トレリス[情報のN番目][今のステージ][前のステージ]，受信信号，軟判定の0の枝，軟判定の1の枝)
void Viterbi(trellis trellis[N][8], int symbol[N]); // ビタビ復号

int main(void) {
	int state; // ステージ
	int input; // 入力
	int output; // 出力
	int next; // 次のステージ
	int i, j; // for
	double ran; // 乱数
	int information[N]; // 情報系列
	int send[N]; // 送信系列
	int receive[N]; // 受信系列
	double sigma; // σ
	trellis trellis[N][8];
	int symbol[N];
	int count;
	double SN;
	double pe;
	int soft0[N][size];
	int soft1[N][size];

	srand((unsigned)time(NULL)); // random初期化

	// 送信情報乱数発生
	for (i = 0; i < N - 2; i++) {
		ran = ((double)rand()) / RAND_MAX;
		if (ran < 0.5) {
			information[i] = 0;
		}
		else {
			information[i] = 1;

		}

		// cin >> information[i];
	}
	information[N - 2] = 0;
	information[N - 1] = 0;

	// 情報出力
	cout << "情報：";
	for (i = 0; i < N; i++) {
		cout << information[i] << "  ";
	}
	cout << "\n";

	state = 0;

	// ビタビ復号を適用できる形に変換する
	for (i = 0; i < N; i++) {
		Convolution(information[i], state, &output, &next); // 畳み込み

		// 状態を移す
		state = next;

		// 出力保存
		send[i] = output;
	}

	// 送信信号出力
	cout << "送信信号：";
	for (i = 0; i < N; i++) {
		cout << send[i] % 2 << send[i] / 2 << "  ";
	}
	cout << "\n\n";

	cout << "\n------------------------- 4値軟判定 ----------------------------\n";
	cout << "デシベル[dB]     σ     誤り率　　実行時間\n";
	for (SN = 1; SN <= 5; SN += 0.5) {
		cout << SN << "      ";
		count = 0; // エラーカウント初期化
		sigma = pow(0.5 * pow(10, -SN / 10), 0.5); //σを逆算する
		cout << std::setprecision(7) << sigma << "  ";

		clock_t start = clock();    // スタート時間

		for (j = 0; j < k; j++) {
			Receive(send, receive, sigma, soft0, soft1); // 受信信号の取得

			// cout << "受信信号：";
			// for (i = 0; i < N; i++) {
			// 	cout << receive[i] % 2 << receive[i] / 2 << "  ";
			// }
			// cout << "\n";

			Trellis(trellis, receive, soft0, soft1);

			// cout << "復号シンボル：";

			Viterbi(trellis, symbol);
			for (i = 0; i < N; i++) {
				// cout << symbol[i] << "   ";
				if (symbol[i] == information[i]) {
				}
				else {
					count += 1;
				}
			}
			// cout << "\n";
		}

		clock_t end = clock();     // 終了時間

		cout << count / ((double)N * (double)k) << "　　　";

		cout << std::setprecision(3) << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	}

	cout << "\n------------------------- 8値軟判定 ----------------------------\n";
	cout << "デシベル[dB]     σ     誤り率　　実行時間\n";
	for (SN = 1; SN <= 5; SN += 0.5) {
		cout << SN << "      ";
		count = 0; // エラーカウント初期化
		sigma = pow(0.5 * pow(10, -SN / 10), 0.5); //σを逆算する
		cout << std::setprecision(7) << sigma << "  ";

		clock_t start = clock();    // スタート時間

		for (j = 0; j < k; j++) {
			Receive8(send, receive, sigma, soft0, soft1); // 受信信号の取得

			// cout << "受信信号：";
			// for (i = 0; i < N; i++) {
			// 	cout << receive[i] % 2 << receive[i] / 2 << "  ";
			// }
			// cout << "\n";

			Trellis(trellis, receive, soft0, soft1);

			// cout << "復号シンボル：";

			Viterbi(trellis, symbol);
			for (i = 0; i < N; i++) {
				// cout << symbol[i] << "   ";
				if (symbol[i] == information[i]) {
				}
				else {
					count += 1;
				}
			}
			// cout << "\n";
		}

		clock_t end = clock();     // 終了時間

		cout <<count / ((double)N * (double)k) << "　　　";

		cout << std::setprecision(3) << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	}

	cout << "\n------------------------- 16値軟判定 ----------------------------\n";
	cout << "デシベル[dB]     σ     誤り率　　実行時間\n";
	for (SN = 1; SN <= 5; SN += 0.5) {
		cout << SN << "      ";
		count = 0; // エラーカウント初期化
		sigma = pow(0.5 * pow(10, -SN / 10), 0.5); //σを逆算する
		cout << std::setprecision(7) << sigma << "  ";

		clock_t start = clock();    // スタート時間

		for (j = 0; j < k; j++) {
			Receive16(send, receive, sigma, soft0, soft1); // 受信信号の取得

			// cout << "受信信号：";
			// for (i = 0; i < N; i++) {
			// 	cout << receive[i] % 2 << receive[i] / 2 << "  ";
			// }
			// cout << "\n";

			Trellis(trellis, receive, soft0, soft1);

			// cout << "復号シンボル：";

			Viterbi(trellis, symbol);
			for (i = 0; i < N; i++) {
				// cout << symbol[i] << "   ";
				if (symbol[i] == information[i]) {
				}
				else {
					count += 1;
				}
			}
			// cout << "\n";
		}

		clock_t end = clock();     // 終了時間

		cout << count / ((double)N * (double)k) << "　　　";

		cout << std::setprecision(3) << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	}

	system("pause"); // 続行するには何かキーを押してください . . .

	return 0;
}

// 畳み込み符号状態遷移図(入力，ステージ，出力，次のステージ)
void Convolution(int input, int state, int *output, int *next) {
	int s[size]; // state
	int n[size]; // next
	int o[size]; // output

	// 2進数に変換
	s[1] = state % 2;
	s[0] = state / 2;

	// 排他的論理和で出力を出す
	o[0] = input ^ s[0] ^ s[1]; // 111
	o[1] = input ^ s[1]; // 101

	*output = o[0] * 2 + o[1]; // 10進数に変換

	// 次のステージを代入
	n[1] = s[0];
	n[0] = input;

	*next = n[0] * 2 + n[1]; // 10進数に変換
}

// 軟判定から求める距離(出力信号，0の枝の軟判定値，1の枝の軟判定値，桁数)
int MetricDistance(int x, int *y0, int *y1, int s) {
	int i; // for
	int dis = 0;
	int a[2];

	// 2進数に変換
	a[1] = x % 2;
	a[0] = x / 2;

	for (i = 0; i < s; i++) {
		if (a[i] == 0) {
			dis += y0[i];
		}
		else {
			dis += y1[i];
		}
	}

	return dis;
}

// 受信信号の取得(送信信号，受信信号，σ，軟判定の0の枝，軟判定の1の枝)
void Receive(int send[N], int receive[N], double sigma, int soft0[N][size], int soft1[N][size]) {
	double x1, x2, x;
	int i, j, m; // for
	int d; // 送信信号の処理
	double ng; // n_g ガウス雑音
	int error_count = 0; // 誤りデータ数
	int s[N][2];
	int r[N][2];

	m = 0;
	for (i = 0; i < N; i++) {
		// 送信信号処理
		s[i][1] = send[i] % 2;
		s[i][0] = send[i] / 2;

		for (j = 0; j < size; j++) {
			if (s[i][j] == 0) {
				d = -1;
			}
			else {
				d = 1;
			}

			x1 = ((double)rand()) / RAND_MAX; // [0, 1]の乱数を発生する
			x2 = ((double)rand()) / RAND_MAX;

			ng = sigma * pow(-2 * log(x1 + pow(10, -21)), 0.5) * cos(2 * pi * x2); // ガウス雑音[-4, 4]

			x = d + ng; // 白色ガウス雑音の発生

			// 受信信号の判定
			if (x >= 0) { // 正
				r[i][j] = 1;
			}
			else { // 負
				r[i][j] = 0;
			}

			if (x >= 1) {
				soft0[i][j] = 3;
				soft1[i][j] = 0;
			}
			else if (x >= 0 and x < 1) {
				soft0[i][j] = 2;
				soft1[i][j] = 1;
			}
			else if (x >= -1 and x < 0) {
				soft0[i][j] = 1;
				soft1[i][j] = 2;
			}
			else if (x <= -1) {
				soft0[i][j] = 0;
				soft1[i][j] = 3;
			}
		}

		receive[i] = r[i][0] * 2 + r[i][1];
	}
}

// 8値軟判定 受信信号の取得(送信信号，受信信号，σ，軟判定の0の枝，軟判定の1の枝)
void Receive8(int send[N], int receive[N], double sigma, int soft0[N][size], int soft1[N][size]) {
	double x1, x2, x;
	int i, j, m; // for
	int d; // 送信信号の処理
	double ng; // n_g ガウス雑音
	int error_count = 0; // 誤りデータ数
	int s[N][2];
	int r[N][2];

	m = 0;
	for (i = 0; i < N; i++) {
		// 送信信号処理
		s[i][1] = send[i] % 2;
		s[i][0] = send[i] / 2;

		for (j = 0; j < size; j++) {
			if (s[i][j] == 0) {
				d = -1;
			}
			else {
				d = 1;
			}

			x1 = ((double)rand()) / RAND_MAX; // [0, 1]の乱数を発生する
			x2 = ((double)rand()) / RAND_MAX;

			ng = sigma * pow(-2 * log(x1 + pow(10, -21)), 0.5) * cos(2 * pi * x2); // ガウス雑音[-4, 4]

			x = d + ng; // 白色ガウス雑音の発生

			// 受信信号の判定
			if (x >= 0) { // 正
				r[i][j] = 1;
			}
			else { // 負
				r[i][j] = 0;
			}

			if (x >= 1.5) {
				soft0[i][j] = 7;
				soft1[i][j] = 0;
			}
			else if (x >= 1 and x < 1.5) {
				soft0[i][j] = 6;
				soft1[i][j] = 1;
			}
			else if (x >= 0.5 and x < 1) {
				soft0[i][j] = 5;
				soft1[i][j] = 2;
			}
			else if (x >= 0 and x < 0.5) {
				soft0[i][j] = 4;
				soft1[i][j] = 3;
			}
			else if (x >= -0.5 and x < 0) {
				soft0[i][j] = 3;
				soft1[i][j] = 4;
			}
			else if (x >= -1 and x < -0.5) {
				soft0[i][j] = 2;
				soft1[i][j] = 5;
			}
			else if (x >= -1.5 and x < -1) {
				soft0[i][j] = 1;
				soft1[i][j] = 6;
			}
			else if (x <= -1.5) {
				soft0[i][j] = 0;
				soft1[i][j] = 7;
			}
		}
		receive[i] = r[i][0] * 2 + r[i][1];
	}
}

// 16値軟判定 受信信号の取得(送信信号，受信信号，σ，軟判定の0の枝，軟判定の1の枝)
void Receive16(int send[N], int receive[N], double sigma, int soft0[N][size], int soft1[N][size]) {
	double x1, x2, x;
	int i, j, m; // for
	int d; // 送信信号の処理
	double ng; // n_g ガウス雑音
	int error_count = 0; // 誤りデータ数
	int s[N][2];
	int r[N][2];

	m = 0;
	for (i = 0; i < N; i++) {
		// 送信信号処理
		s[i][1] = send[i] % 2;
		s[i][0] = send[i] / 2;

		for (j = 0; j < size; j++) {
			if (s[i][j] == 0) {
				d = -1;
			}
			else {
				d = 1;
			}

			x1 = ((double)rand()) / RAND_MAX; // [0, 1]の乱数を発生する
			x2 = ((double)rand()) / RAND_MAX;

			ng = sigma * pow(-2 * log(x1 + pow(10, -21)), 0.5) * cos(2 * pi * x2); // ガウス雑音[-4, 4]

			x = d + ng; // 白色ガウス雑音の発生

			// 受信信号の判定
			if (x >= 0) { // 正
				r[i][j] = 1;
			}
			else { // 負
				r[i][j] = 0;
			}
			/*
			-0.53456 0.124639 0.523809 0.849449 1.160313 1.479969 1.66123 2.13276
			*/

			if (x >= 1.75) {
				soft0[i][j] = 15;
				soft1[i][j] = 0;
			}
			else if (x >= 1.5 and x < 1.75) {
				soft0[i][j] = 14;
				soft1[i][j] = 1;
			}
			else if (x >= 1.25 and x < 1.5) {
				soft0[i][j] = 13;
				soft1[i][j] = 2;
			}
			else if (x >= 1.0 and x < 1.25) {
				soft0[i][j] = 12;
				soft1[i][j] = 3;
			}
			else if (x >= 0.75 and x < 1.0) {
				soft0[i][j] = 11;
				soft1[i][j] = 4;
			}
			else if (x >= 0.5 and x < 0.75) {
				soft0[i][j] = 10;
				soft1[i][j] = 5;
			}
			else if (x >= 0.25 and x < 0.5) {
				soft0[i][j] = 9;
				soft1[i][j] = 6;
			}
			else if (x >= 0 and x < 0.25) {
				soft0[i][j] = 8;
				soft1[i][j] = 7;
			}
			else if (x >= -0.25 and x < 0) {
				soft0[i][j] = 7;
				soft1[i][j] = 8;
			}
			else if (x >= -0.5 and x < -0.25) {
				soft0[i][j] = 6;
				soft1[i][j] = 9;
			}
			else if (x >= -0.75 and x < -0.5) {
				soft0[i][j] = 5;
				soft1[i][j] = 10;
			}
			else if (x >= -1.0 and x < -0.75) {
				soft0[i][j] = 4;
				soft1[i][j] = 11;
			}
			else if (x >= -1.25 and x < -1.0) {
				soft0[i][j] = 3;
				soft1[i][j] = 12;
			}
			else if (x >= -1.5 and x < -1.25) {
				soft0[i][j] = 2;
				soft1[i][j] = 13;
			}
			else if (x >= -1.75 and x < -1.5) {
				soft0[i][j] = 1;
				soft1[i][j] = 14;
			}
			else if (x <= -2.0) {
				soft0[i][j] = 0;
				soft1[i][j] = 15;
			}
		}
		receive[i] = r[i][0] * 2 + r[i][1];
	}
}

// トレリス図の作成(トレリス[情報のN番目][今のステージ][前のステージ2桁目]，受信信号，軟判定の0の枝，軟判定の1の枝)
void Trellis(trellis trellis[N][8], int receive[N], int soft0[N][size], int soft1[N][size]) {
	int state0 = 0, state1 = 1, state2 = 2, state3 = 3; // ステージ
	int input0 = 0, input1 = 1; // 入力0，入力1
	int output; // 出力
	int next; // 次のステージ
	int i, j; // for
	int dis;

	for (i = 0; i < N; i++) { // ハミング距離初期化
		for (j = 0; j < 8; j++) {
			trellis[i][j].hamming = NIL;
		}
	}

	// 最初の2個を先に作る
	Convolution(input0, state0, &output, &next); // 状態遷移
	trellis[0][next].hamming = MetricDistance(output, soft0[0], soft1[0], size); // 誤り確率を計算する

	trellis[0][next].now = next; // 今のステージを保存する
	trellis[0][next].before = state0; // 前のステージを保存する
	// cout << trellis[0][0].hamming<< "   ";

	Convolution(input1, state0, &output, &next);
	trellis[0][next].hamming = MetricDistance(output, soft0[0], soft1[0], size);
	trellis[0][next].now = next;
	trellis[0][next].before = state0;
	// cout << trellis[0][2].hamming << "   \n";

	Convolution(input0, state0, &output, &next);
	trellis[1][next].hamming = trellis[0][state0].hamming + MetricDistance(output, soft0[1], soft1[1], size);
	trellis[1][next].now = next;
	trellis[1][next].before = state0;
	// cout << trellis[1][0].hamming << "   ";

	Convolution(input1, state0, &output, &next);
	trellis[1][next].hamming = trellis[0][state0].hamming + MetricDistance(output, soft0[1], soft1[1], size);
	trellis[1][next].now = next;
	trellis[1][next].before = state0;
	// cout << trellis[1][2].hamming << "   ";

	Convolution(input0, state2, &output, &next);
	trellis[1][next].hamming = trellis[0][state2].hamming + MetricDistance(output, soft0[1], soft1[1], size);
	trellis[1][next].now = next;
	trellis[1][next].before = state2;
	// cout << trellis[1][1].hamming << "   ";

	Convolution(input1, state2, &output, &next);
	trellis[1][next].hamming = trellis[0][state2].hamming + MetricDistance(output, soft0[1], soft1[1], size);
	trellis[1][next].now = next;
	trellis[1][next].before = state2;
	// cout << trellis[1][3].hamming << "   \n";

	for (i = 2; i < N; i++) {
		for (j = 0; j < 4; j++) {
			Convolution(input0, j, &output, &next);
			if (trellis[i - 1][j].hamming < trellis[i - 1][j + 4].hamming) {
				dis = trellis[i - 1][j].hamming;
			}
			else {
				dis = trellis[i - 1][j + 4].hamming;
			}

			if (trellis[i][next].hamming == NIL) {
				trellis[i][next].hamming = dis + MetricDistance(output, soft0[i], soft1[i], size);
				trellis[i][next].now = next;
				trellis[i][next].before = j;
			}
			else {
				trellis[i][next + 4].hamming = dis + MetricDistance(output, soft0[i], soft1[i], size);
				trellis[i][next + 4].now = next;
				trellis[i][next + 4].before = j;
			}

			Convolution(input1, j, &output, &next);
			if (trellis[i - 1][j].hamming < trellis[i - 1][j + 4].hamming) {
				dis = trellis[i - 1][j].hamming;
			}
			else {
				dis = trellis[i - 1][j + 4].hamming;
			}

			if (trellis[i][next].hamming == NIL) {
				trellis[i][next].hamming = dis + MetricDistance(output, soft0[i], soft1[i], size);
				trellis[i][next].now = next;
				trellis[i][next].before = j;
			}
			else {
				trellis[i][next + 4].hamming = dis + MetricDistance(output, soft0[i], soft1[i], size);
				trellis[i][next + 4].now = next;
				trellis[i][next + 4].before = j;
			}
		}
		// for (j = 0; j < 8; j++) {
		//	 cout << trellis[i][j].hamming << "   ";
		// }
		// cout << "\n";
	}
}

// ビタビ復号
void Viterbi(trellis trellis[N][8], int symbol[N]) {
	int i; // for
	int before;
	int state = 0;

	// 最後00で終わる
	symbol[N - 1] = 0;
	symbol[N - 2] = 0;


	for (i = N - 1; i >= 0; i--) {
		if (trellis[i][state].hamming < trellis[i][state + 4].hamming) {
			symbol[i] = state / 2;
			state = trellis[i][state].before;
		}
		else {
			symbol[i] = state / 2;
			state = trellis[i][state + 4].before;
		}

	}
}
