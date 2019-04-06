// BCJR.cpp

#define  _CRT_SECURE_NO_WARNINGS
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<cstring>
#include<iostream>
#include<fstream>
using namespace std;

const int message_length = 10000 + 2; //the length of message
const int codeword_length = 2 * message_length; //the length of codeword

												// code_rate
float code_rate = (float)message_length / (float)codeword_length;

// channel coefficient
#define INF 2097152
#define pi 3.1415926
double N0, sgm;

const int state_num = message_length - 2;
const int state_type = 4;
const int s_num = 4;   // for (7,5) conv. code
//double state_table[state_num][state_type];//state table, the size should be defined yourself
										  //int state_num;//the number of the state of encoder structure

int message[message_length], codeword[codeword_length];//message and codeword
int re_codeword[codeword_length];//the received codeword
int de_message[message_length];//the decoding message

double tx_symbol[codeword_length][2];//the transmitted symbols
double rx_symbol[codeword_length][2];//the received symbols

void statetable(int n);
void encoder();
void modulation();
void demodulation();
void channel();
void decoder();

void main()
{
	int i;
	float SNR, start, finish;
	long int bit_error, seq, seq_num;
	double BER;
	double progress;
	//fstream
	ofstream fout("data_BCJR");
	///////////////////////////////////////////////////

	//generate state table
	statetable(0);

	//random seed
	srand((int)time(0));

	//input the SNR and frame number
	printf("\nEnter start SNR: ");
	scanf("%f", &start);
	printf("\nEnter finish SNR: ");
	scanf("%f", &finish);
	printf("\nPlease input the number of message: ");
	scanf("%d", &seq_num);

	for (SNR = start; SNR <= finish; SNR++)
	{
		//channel noise
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2);

		bit_error = 0;

		for (seq = 1; seq <= seq_num; seq++)
		{
			//generate binary message randomly
			/****************
			Pay attention that message is appended by 0 whose number is equal to the state of encoder structure.
			****************/

			for (i = 0; i<message_length - s_num; i++)
			{
				message[i] = rand() % 2;
			}
			for (i = message_length - s_num; i<message_length; i++)
			{
				message[i] = 0;
			}
			/*
			message[0] = 1; message[1] = 1; message[2] = 0; message[3] = 1; message[4] = 1;
			message[5] = 0; message[6] = 0; message[7] = 0; message[8] = 1; message[9] = 0;
			message[10] = 0; message[11] = 0;
			*//////////

			//convolutional encoder
			encoder();

			//BPSK modulation
			modulation();

			//AWGN channel
			channel();

			//BPSK demodulation, it's needed in hard-decision Viterbi decoder
			//demodulation();

			//convolutional decoder
			decoder();

			//calculate the number of bit error
			for (i = 0; i<message_length; i++)
			{
				if (message[i] != de_message[i])
					bit_error++;
			}

			progress = (double)(seq * 100) / (double)seq_num;

			//calculate the BER
			BER = (double)bit_error / (double)(message_length*seq);

			//print the intermediate result
			printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\r", progress, SNR, bit_error, BER);
		}

		//calculate the BER
		BER = (double)bit_error / (double)(message_length*seq_num);

		//print the final result
		printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n", progress, SNR, bit_error, BER);
		fout << SNR << ',' << BER << endl;
	}
	fout.close();
	system("pause");
}

const int edge_num = 8;
struct edge
{
	edge()
	{
		P = 0.0;
		A = 0.0;
		B = 0.0;
	}
	edge(int in, int cs, int ns, int out, int id)
	{
		IN = in;
		current_state = cs;
		next_state = ns;
		OUT = out;
		ID = id;
	}
	int IN;
	int current_state;
	int next_state;
	int OUT;
	int ID;
	double P;
	double A;
	double B;
}path[state_num+5][edge_num];

void statetable(int n)
{
	path[n][0] = edge(0, 0, 0, 0, 1);
	path[n][1] = edge(1, 0, 2, 3, 2);
	path[n][2] = edge(1, 1, 2, 0, 3);
	path[n][3] = edge(0, 1, 0, 3, 4);
	path[n][4] = edge(1, 2, 3, 1, 5);
	path[n][5] = edge(0, 2, 1, 2, 6);
	path[n][6] = edge(0, 3, 1, 1, 7);
	path[n][7] = edge(1, 3, 3, 2, 8);
}

void encoder()
{
	//convolution encoder, the input is message[] and the output is codeword[]
	//(7,5)8 conv. code
	int s0 = 0;
	int s1 = 0;

	int cnt = 0;
	for (int i = 0; i<message_length; i++)
	{
		codeword[cnt++] = (message[i] + s0 + s1) % 2;
		codeword[cnt++] = (message[i] + s1) % 2;
		s1 = s0;
		s0 = message[i];
	}
}

void modulation()
{
	//BPSK modulation
	int i;

	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i<codeword_length; i++)
	{
		tx_symbol[i][0] = -1 * (2 * codeword[i] - 1);
		tx_symbol[i][1] = 0;
	}
}
void channel()
{
	//AWGN channel
	int i, j;
	double u, r, g;

	for (i = 0; i<codeword_length; i++)
	{
		for (j = 0; j<2; j++)
		{
			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			r = sgm*sqrt(2.0*log(1.0 / (1.0 - u)));

			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			g = (float)r*cos(2 * pi*u);

			rx_symbol[i][j] = tx_symbol[i][j] + g;
		}
	}
}
void demodulation()
{
	int i;
	double d1, d2;
	for (i = 0; i<codeword_length; i++)
	{
		d1 = (rx_symbol[i][0] - 1)*(rx_symbol[i][0] - 1) + rx_symbol[i][1] * rx_symbol[i][1];
		d2 = (rx_symbol[i][0] + 1)*(rx_symbol[i][0] + 1) + rx_symbol[i][1] * rx_symbol[i][1];
		if (d1<d2)
			re_codeword[i] = 0;
		else
			re_codeword[i] = 1;
	}
}
//int edge_id[state_num][state_type];
double Pa[2] = { 0.5 , 0.5 };
double state_tableA[state_num+5][state_type];
double state_tableB[state_num+5][state_type];
void decoder()
{
	//statetable();
	int BPSK_symbol[2];
	BPSK_symbol[0] = 1;
	BPSK_symbol[1] = -1;
	memset(state_tableA, 0, sizeof(state_tableA)); //Only work for double when it is 0!
	memset(state_tableB, 0, sizeof(state_tableB));
	state_tableA[0][0] = 1.0;
	int cnt = 0;
	for (int i = 0; i<codeword_length; i += 2) //i+= code_rate (rewrite)
	{
		// (rewrite)
		statetable(cnt);
		double codeword_x1 = rx_symbol[i][0];			// 2 * re_codeword[i] + re_codeword[i + 1];
		double codeword_y1 = rx_symbol[i][1];
		double codeword_x2 = rx_symbol[i + 1][0];
		double codeword_y2 = rx_symbol[i + 1][1];
		double D1[2]; double D2[2];
		D1[0] = (codeword_x1 - BPSK_symbol[0]) * (codeword_x1 - BPSK_symbol[0]) + codeword_y1 * codeword_y1;
		D1[1] = (codeword_x1 - BPSK_symbol[1]) * (codeword_x1 - BPSK_symbol[1]) + codeword_y1 * codeword_y1;
		D2[0] = (codeword_x2 - BPSK_symbol[0]) * (codeword_x2 - BPSK_symbol[0]) + codeword_y2 * codeword_y2;
		D2[1] = (codeword_x2 - BPSK_symbol[1]) * (codeword_x2 - BPSK_symbol[1]) + codeword_y2 * codeword_y2;
		double P1[2]; double P2[2];
		double sum_temp;
		P1[0] = 1 / (pi * N0) * exp(-(D1[0]) / N0);
		P1[1] = 1 / (pi * N0) * exp(-(D1[1]) / N0);
		sum_temp = P1[0] + P1[1];
		P1[0] /= sum_temp;
		P1[1] /= sum_temp;
		P2[0] = 1 / (pi * N0) * exp(-(D2[0]) / N0);
		P2[1] = 1 / (pi * N0) * exp(-(D2[1]) / N0);
		sum_temp = P2[0] + P2[1];
		P2[0] /= sum_temp;
		P2[1] /= sum_temp;
		for (int j = 0; j<edge_num; j++)
		{
			// (rewrite)
			// Pch
			double Pch1 = P1[(path[cnt][j].OUT >> 1) & 1];
			double Pch2 = P2[(path[cnt][j].OUT & 1)];
			double Pg = Pa[path[cnt][j].IN] * Pch1 * Pch2;
			path[cnt][j].P = Pg;
			//state_table[cnt][path[j].next_state] = min(state_table[cnt][path[j].next_state],state_table[cnt-1][path[j].current_state]);
		}
		cnt++;
	}


	state_tableA[0][0] = 1.0;
	for (int i = 0; i < cnt; i++)
	{
		for (int j = 0; j < edge_num; j++)
			state_tableA[i + 1][path[i][j].next_state] += state_tableA[i][path[i][j].current_state] * path[i][j].P;
		double sum_temp = 0.0;
		for (int k = 0; k < state_type; k++)
			sum_temp += state_tableA[i+1][k];
		for (int k = 0; k < state_type; k++)
			state_tableA[i + 1][k] /= sum_temp;
	}
	state_tableB[cnt][0] = 1.0; //?? end?
	for (int i = cnt - 1; i >= 0; i--)
	{
		for (int j = 0; j < edge_num; j++)
			state_tableB[i][path[i][j].current_state] += state_tableB[i + 1][path[i][j].next_state] * path[i][j].P;
		double sum_temp = 0.0;
		for (int k = 0; k < state_type; k++)
			sum_temp += state_tableB[i][k];
		for (int k = 0; k < state_type; k++)
			state_tableB[i][k] /= sum_temp;
	}
	cnt--;

	//de_message
	for (int i = 0; i < message_length; i++)
	{
		double Pp[2] = { 0.0 , 0.0 };
		for (int j = 0; j < edge_num; j++)
			Pp[path[i][j].IN] += state_tableA[i][path[i][j].current_state] * path[i][j].P * state_tableB[i + 1][path[i][j].next_state];
		de_message[i] = Pp[0] > Pp[1] ? 0 : 1;
	}
}
