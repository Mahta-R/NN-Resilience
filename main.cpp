#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <iomanip>
#include <cstdlib>
#include <sstream>

using namespace std;

const float PI = 3.1415;

struct point
{
	float x;
	float y;
};

const int N = 40; 
const int N1 = N*N; //number of spins layer1

const int O = 20;
const int O1 = O*O; //number of spins layer2

const int M = 3; // number of patterns

int S1[N][N]; // spin network layer1
float W1[N1][N1]; // Weigths layer 1
int S2[O][O]; // layer 2
float W2[O1][O1]; // layer 2

float W12[O1][2]; // between 2 layers

int ***P; //array of patterns

float f(float x)  // Fermi function
{
	float a = (float) 1 / (1+exp(-x));
	return a;
}

point rotation(point a, float t) // rotate with angle t
{
	point b;

	b.x = cos(t) * a.x + sin(t) * a.y;
	b.y = -sin(t) * a.x + cos(t) * a.y;

	return b;
}
	

void calW()  // read patterns from files
{

	ifstream myFiles[M];
	
	for(int i = 0; i < M; i++)
	{
		stringstream filename;
		filename << "info" << i  << ".txt" << ends;
		myFiles[i].open(filename.str().c_str());
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				myFiles[i] >> P[j][k][i];
			}
		}
	}
}



void Wl1() // calculate weights for layer 1    .....   initialize interaction weights randomly
{
	for (int i = 0; i < N1; i++)
	{
		for (int j = 0; j < N1; j++)
		{
			int a1 = i / N;
			int a2 = i % N;
			int b1 = j / N;
			int b2 = j % N;
			for (int k = 0; k < M; k++)
			{
				
				W1[i][j] += P[a1][a2][k] * P[b1][b2][k];
			}

			W1[i][j] = (float) W1[i][j] / N1;
		}
	}

	for (int i = 0; i < O1; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			W12[i][j] = (float) (rand() % 3000000) / 10000;
			W12[i][j] = (float) W12[i][j] / (O1*2);
			//cout << W12[i][j] << "\t";
		}
	}
}

void RGs() // spins for second layer ...... change interaction weights  based on restricted boltzmann learning
{
	for (int i = 0; i < O1; i++)
	{
		int a1 = i / O;
		int a2 = i % O;

		int j1 = 2*i;
		int j2 = j1 + 1;

		int b1 = j1 / N;
		int b2 = j1 % N;

		int c1 = j2 / N;
		int c2 = j2 % N;

		float temp = f(W12[i][0]*S1[b1][b2] + W12[i][1]*S1[c1][c2]);

		if (temp > 0.5)
			S2[a1][a2] = 1;
		if (temp <= 0.5)
			S2[a1][a2] = 0;
		
		S2[a1][a2] = S1[c1][c2];

		float sn1 = f(W12[i][0]*S2[a1][a2]);
		float sn2 = f(W12[i][1]*S2[a1][a2]);
		float sn3 = f(W12[i][0]*sn1 + W12[i][1]*sn2);

		float dw = 0.1 * (S1[b1][b2]*S2[a1][a2] + S1[c1][c2]*S2[a1][a2] - sn1*sn3 - sn2*sn3);
		
		W12[i][0] += dw;
		W12[i][1] += dw;
		
	}
}  	
	
	

void RGw() // weights for second layer
{
	for (int i = 0;	i < O1; i++)
	{
		for (int j = 0; j < O1; j++)
		{
			int a = i*O + j;
			int j1 = 2*a;
			int j2 = j1 + 1;
			
			int b1 = j1 / N;
			int b2 = j1 % N;

			int c1 = j2 / N;
			int c2 = j2 % N; 
			
			W2[i][j] = W12[0][i]*W1[b1][b2] + W12[1][i]*W1[c1][c2];
		}
	}
}
	 

float hl1(int s1, int s2) // 
{
	float s;
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			s += W1[i*N+j][s1*N+j]*S1[i][j];
		}
	}

	s = -0.5*s;

	return s;
}

float hl2(float t)
{
	int St[O][O];
	point p, pn;
	
	for (int i = 0; i < O; i++)
	{
		for (int j = 0; j < O; j++)
		{
			if (S2[i][j] == 1)
			{
				p.x = i;
				p.y = j;
				pn = rotation(p,t);
				
				int i1 = pn.x;
				int j1 = pn.y;

				St[i1][j1] = 1;
			}
		}
	}

	float s1 = 0;
	float s1t = 0;
	float s2 = 0;
	float s2t = 0;

	for (int i = 0; i < O; i++)
	{
		for (int j = 0; j < O; j++)
		{
			if (St[i][j] == 1 && S2[i][j] == 0)
			{
				int k1 = i*O + j;
				
				for (int k = 0; k < O1; k++)
				{
					int i1 = k / O;
					int j1 = k % O;
					s1 += W2[k1][k]*St[i][j]*St[i1][j1];
				}
				
				s1t += s1;
				s1 = 0;
			}

			if (S2[i][j] == 1 && St[i][j] == 0)
			{
				int k1 = i*O + j;
				
				for (int k = 0; k < O1; k++)
				{
					int i1 = k / O;
					int j1 = k % O;
					s2 += W2[k1][k]*S2[i][j]*S2[i1][j1];
				}
				
				s2t += s2;
				s2 = 0;
			}
		}
	}

	float dw = s1t - s2t;
	return dw;
}
	
	
	

int main()
{
	ofstream data("r1.txt");
	
	for (int i = 0; i < N1; i++)  // set weights 0
	{
		for (int j = 0; j < N1; j++)
		{
			W1[i][j] = 0;
		}
		
	}

	
	

	P = new int**[N]; // creat 3d array for patterns
	for(int i = 0; i < N; i++)
	{
    		P[i] = new int*[N];
		
		for(int j = 0; j < N; j++)
		{
        		P[i][j] = new int[M];
		}
	}

	calW(); // read patterns

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			S1[i][j] = P[i][j][0];
		}
	}
 
	Wl1(); // set weigths in layer 1

	RGw(); // set weigths in layer 2

	ifstream num1("i2t.txt");

	for (int i = 0; i < N; i++) // input to layer 1
	{
		for (int j = 0; j < N; j++)
		{
			num1 >> S1[i][j];
		}
	}

	RGs(); // change spins layer 2 ... change interaction weigths

	float s1 = 0;
	float s2 = 0;

	for (int a = 0; a < 100; a++)
	{
		s1 = 0;
		s2 = 0;
		for (int a1 = 0; a1 < N1; a1++)
		{
			for (int a2 = 0; a2 < N1; a2++)
			{
				W1[a1][a2] = 1.005 * W1[a1][a2];
				s1 += W1[a1][a2] * W1[a1][a2];
				s2 += W1[a1][a2];
			}
		}

		for (int a1 = 0; a1 < O1; a1++)
		{
			for (int a2 = 0; a2 < O1; a2++)
			{
				W2[a1][a2] = 1.005 * W2[a1][a2];
				s1 += W2[a1][a2] * W2[a1][a2];
				s2 += W2[a1][a2];
			}
		}

		s1 = (float) s1 / (N1+O1);
		s2 = (float) s2 / (N1+O1);

	for (int z = 0; z < 100; z++)
	{
		// choose random spin and check layer 1
		int i = rand() % N;
		int j = rand() % N;

		float cw = hl1(i,j);
		
		float k = (float) (rand() % 1000) / 1000;

		if (k < exp(cw))
			S1[i][j] = 1;
		else
			S1[i][j] = 0;

		// choose an angle and check layer 2

		RGs();
		RGw();

		
		
	}
		

	float s = 0;	

	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			s += S1[i][j] * P[i][j][2] * s2 ;
		}
		
		//cout << endl;
	}

	data <<  (float) ((float)s1 / s2)*1000 - 0.1 << "\t" << (float) s / (80*s2) << endl;

	}

	return 0;
}
	


