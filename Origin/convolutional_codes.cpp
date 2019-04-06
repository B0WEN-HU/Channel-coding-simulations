

#define  _CRT_SECURE_NO_WARNINGS
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cmath>

const int message_length = 1000000; //the length of message
const int codeword_length = 2 * message_length; //the length of codeword

// code_rate
float code_rate = (float)message_length / (float)codeword_length;

// channel coefficient
#define INF 2097152
#define pi 3.1415926
double N0, sgm;

int state_num = message_length+2;
int state_type = 4;
int state_table[state_num][state_type];//state table, the size should be defined yourself
//int state_num;//the number of the state of encoder structure

int message[message_length], codeword[codeword_length];//message and codeword
int re_codeword[codeword_length];//the received codeword
int de_message[message_length];//the decoding message

double tx_symbol[codeword_length][2];//the transmitted symbols
double rx_symbol[codeword_length][2];//the received symbols

void statetable();
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

	//generate state table
	statetable();

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

		for (seq = 1; seq<=seq_num; seq++)
		{
			//generate binary message randomly
			/****************
			Pay attention that message is appended by 0 whose number is equal to the state of encoder structure.
			****************/
			for (i = 0; i<message_length - state_num; i++)
			{
				message[i] = rand() % 2;
			}
			for (i = message_length - state_num; i<message_length; i++)
			{
				message[i] = 0;
			}

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
	}
	system("pause");
}

struct edge
{
	edge(int in,int cs,int ns,int out,int id)
	{
		IN = in;
		current_state = cs;
		next_state = ns;
		OUT = out;
		ID = id;
	}
	int IN;
	int current_state;
	int next_stete;
	int OUT;
	int ID;
}path[state_type*2+2];
const int edge_num = 8;
void statetable()
{
	path[0] = edge(0,0,0,0,1);
	path[1] = edge(1,0,2,3,2);
	path[2] = edge(0,1,2,0,3);
	path[3] = edge(1,1,0,3,4);
	path[4] = edge(0,2,3,1,5);
	path[5] = edge(1,2,1,2,6);
	path[6] = edge(0,3,1,1,7);
	path[7] = edge(1,3,3,2,8);
}

void encoder()
{
	//convolution encoder, the input is message[] and the output is codeword[]
	//(7,5)8 conv. code
	int s0 = 0;
	int s1 = 1;
	
	int cnt = 0;
	for (int i = 0; i<message_length; i++)
	{
		codeword[cnt++] = (message[i] + s0 + s1) % 2;
		codeword[cnt++] = (message[i] + s1) % 2;
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
		tx_symbol[i][1]=0;
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
			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			r=sgm*sqrt(2.0*log(1.0/(1.0-u)));

			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			g=(float)r*cos(2*pi*u);

			rx_symbol[i][j]=tx_symbol[i][j]+g;
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
int edge_id[state_num][state_type];
void decoder()
{
	statetable();
	for (int i=1; i<state_type; i++)
		for (int j=0; j< state_num; j++)
			state_table[j][i] = INF;
	state_table[0][0] = 0;
	int cnt = 1;
	for (int i=0; i<codeword_length; i+=2) //i+= code_rate (rewrite)
	{
		// (rewrite)
		int codeword_temp = 2*re_codeword[i+1] + re_codeword[i];
		for (int j=0; j<edge_num; j++)
		{
			// (rewrite)
			int dis =   (codeword_temp ^ path[j].OUT);
			dis = (dis & 1) + ((dis >> 1) & 1);
			//state_table[cnt][path[j].next_state] = min(state_table[cnt][path[j].next_state],state_table[cnt-1][path[j].current_state]);
			if(state_table[cnt-1][path[j].current_state] + dis < state_table[cnt][path[j].next_state)
			{
				state_table[cnt][path[j].next_state = state_table[cnt-1][path[j].current_state] + dis;
				edge_id[cnt][path[j].next_state] = path[j].id;
			}
			cnt++;
		}
	}
	cnt--;
	// rewind the decodeword using edge_id!
	int end_state = 0;
	int min_state = INF;
	
	for (int i=0; i<state_type; i++)
		if(min_state > state_table[cnt][i])
		{
			end_state = i;
			min_state = state_table[cnt][i];
		}
	for (int i = message_length-1; i>0; i--)
	{
		de_message[i] = path[edge_id[i+1][end_state]-1].IN;
		end_state = path[edge_id[i+1][end_state]-1].current_state;
	}
}
