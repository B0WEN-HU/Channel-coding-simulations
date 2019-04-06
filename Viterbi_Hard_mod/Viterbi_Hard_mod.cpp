// Viterbi_Hard.cpp

#define  _CRT_SECURE_NO_WARNINGS
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<fstream>
#include <bitset>
#include"statetable.hpp"

using namespace std;

const int message_length = 10000 + 2; //the length of message
const int codeword_length = 2 * message_length; //the length of codeword

												// code_rate
float code_rate = (float)message_length / (float)codeword_length;

// channel coefficient
#define INF 2097152
#define pi 3.1415926
double N0, sgm;


/************************只需要修改这个参数************************/
const int cov_g1 = 07, cov_g2 = 05;// for example:(7, 5),(15, 13),(23, 35),(171, 133)前面加0表示8进制
/************************只需要修改这个参数************************/

//the argument about statetabe build
#define Max_state_type 128//2^7
edge path[Max_state_type * 2 + 2];
int nx_id[Max_state_type * 2 + 2];// for encoder use
const int state_num = message_length;
int s_num = 0;
int state_type =0; 
int edge_num = 0;// for (7,5) conv. 
//the argument about statetabe build

int state_table[state_num+5][Max_state_type];//state table, the size should be defined yourself
									   //int state_num;//the number of the state of encoder structure

int message[message_length], codeword[codeword_length];//message and codeword
int re_codeword[codeword_length];//the received codeword
int de_message[message_length];//the decoding message

double tx_symbol[codeword_length][2];//the transmitted symbols
double rx_symbol[codeword_length][2];//the received symbols

void statetable(int, int);
void encoder();
void modulation();
void demodulation();
void channel();
void decoder();

void display(int);
int binarySearch(edge array[], int low, int high, int target);
void main()
{
	int i;
	float SNR, start, finish;
	long int bit_error, seq, seq_num;
	double BER;
	double progress;
	//fstream
	ofstream fout("data_Viterbi_Hard");
	///////////////////////////////////////////////////

	//generate state table
	//statetable();
	statetable(cov_g1, cov_g2);
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
			demodulation();

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

//void statetable()
//{
//	path[0] = edge(0, 0, 0, 0, 1);
//	path[1] = edge(1, 0, 2, 3, 2);
//	path[2] = edge(1, 1, 2, 0, 3);
//	path[3] = edge(0, 1, 0, 3, 4);
//	path[4] = edge(1, 2, 3, 1, 5);
//	path[5] = edge(0, 2, 1, 2, 6);
//	path[6] = edge(0, 3, 1, 1, 7);
//	path[7] = edge(1, 3, 3, 2, 8);
//}

void statetable(int cov_g1, int cov_g2)
{
	int in = 0, cs = 0, ns = 0, out = 0, id = 0;
	int rs_num = count_reg(cov_g1);//寄存器数目
	int out1, out2;//供计算out用
	bitset<16> g1(cov_g1);//供计算out用
	bitset<16> g2(cov_g2);//供计算out用

	//计算状态表
	for (int i = 0; i < pow(2, rs_num); i++)
	{
		cs = i;
		for (in = 0; in < 2; in++)
		{
			path[id].IN = in;
			path[id].current_state = cs;
			ns = (path[id].current_state >> 1) + in*pow(2, rs_num - 1);//寄存器状态移位
			path[id].next_state = ns;
			//****************************//计算输出
			out1 = in;
			out2 = in;
			for (int j = rs_num - 1; j >= 0; j--)
			{
				if (g1[j] == 1)
				{
					out1 = out1 ^ bitset<16>(cs)[j];

				}
				if (g2[j] == 1)
				{
					out2 = out2 ^ bitset<16>(cs)[j];
				}
			}
			path[id].OUT = out1 * 2 + out2;//output
			id++;
		}
	}

	//给id赋值
	for (int i = 0; i < pow(2, rs_num)*2; i++)
	{
		path[i].ID = i + 1;
	}

	//给s_num, state_type, edge_num 赋值
	s_num = rs_num;
	state_type = pow(2, rs_num);
	edge_num = pow(2, rs_num) * 2;

	//计算下一个id

	for (int i = 0; i <  pow(2, rs_num) * 2; i++)
	{
		id = binarySearch(path, 0, pow(2, rs_num) * 2, path[i].next_state);
		nx_id[i] = id;
	}

	display(rs_num);
}

void display(int s_num)
{	
	cout << setw(4) << "IN" << " ";
	cout << setw(4) << "CS" << " ";
	cout << setw(4) << "NS" << " ";
	cout << setw(4) << "Out" << " ";
	cout << setw(4) << "ID" << " ";
	cout << setw(4) << "ns_id" << endl;
	for (int i = 0; i < pow(2, s_num)*2; i++)
	{
		cout << setw(4) << path[i].IN << " ";
		cout << setw(4) << path[i].current_state << " ";
		cout << setw(4) << path[i].next_state << " "; 
		cout << setw(4) << path[i].OUT << " ";
		cout << setw(4) << path[i].ID << " ";
		cout << setw(4) << nx_id[i] << endl;
	}
}

void encoder()
{
	//convolution encoder, the input is message[] and the output is codeword[]
	//(7,5)8 conv. code
	int cs_id = 0;//current state id
	int cd_id = 0;//cordword id

	for (int i = 0; i<message_length; i++)
	{
		if (message[i] == 0)
		{
			codeword[cd_id + 1] = bitset<16>(path[cs_id].OUT)[0];
			codeword[cd_id] = bitset<16>(path[cs_id].OUT)[1];
			cs_id = nx_id[cs_id];
		}
		else if (message[i] == 1)
		{
			codeword[cd_id + 1] = bitset<16>(path[cs_id + 1].OUT)[0];
			codeword[cd_id] = bitset<16>(path[cs_id + 1].OUT)[1];
			cs_id = nx_id[cs_id + 1];
		}
		cd_id = cd_id + 2;
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
int edge_id[state_num+5][Max_state_type];
void decoder()
{
	//statetable();
	for (int i = 0; i<state_type; i++)
		for (int j = 0; j< state_num + 5; j++)
			state_table[j][i] = INF;
	state_table[0][0] = 0;
	int cnt = 1;
	for (int i = 0; i<codeword_length; i += 2) //i+= code_rate (rewrite)
	{
		// (rewrite)
		int codeword_temp = 2 * re_codeword[i] + re_codeword[i + 1];
		for (int j = 0; j<edge_num; j++)
		{
			// (rewrite)
			int dis = (codeword_temp ^ path[j].OUT);
			dis = (dis & 1) + ((dis >> 1) & 1);
			//state_table[cnt][path[j].next_state] = min(state_table[cnt][path[j].next_state],state_table[cnt-1][path[j].current_state]);
			if (state_table[cnt - 1][path[j].current_state] + dis < state_table[cnt][path[j].next_state])
			{
				state_table[cnt][path[j].next_state] = state_table[cnt - 1][path[j].current_state] + dis;
				edge_id[cnt][path[j].next_state] = j;
			}
		}
		cnt++;
	}
	cnt--;
	// rewind the decodeword using edge_id!
	int end_state = 0;
	int min_state = INF;

	for (int i = 0; i<state_type; i++)
		if (min_state > state_table[cnt][i])
		{
			end_state = i;
			min_state = state_table[cnt][i];
		}
	for (int i = message_length - 1; i>=0; i--)
	{
		de_message[i] = path[edge_id[i + 1][end_state]].IN;
		end_state = path[edge_id[i + 1][end_state]].current_state;
	}
}

int binarySearch(edge array[], int low, int high, int target)
{
	for (int i = low; i < high; i++)
	{
		if (target == array[i].current_state)
		{
			return i;
		}
	}
	return 0;
}
