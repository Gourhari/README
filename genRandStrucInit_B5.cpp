#include <bits/stdc++.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <string>

#define atomCount 5  //For B36
#define structCount 15
#define range 3

using namespace std;

typedef struct atom{
	char symbol[3];
	double coX;
	double coY;
	double coZ;
	double velX;
	double velY;
	double velZ;
}atom;

typedef struct particle{
	struct atom atoms[atomCount];
	struct atom pBest[atomCount];
	double cost;
	double pBestCost;
}particle;

std::string to_string(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

int main(int argc, char const *argv[])
{
	srand (static_cast <unsigned> (time(0)));
	for(int i = 1; i <= structCount; i++)
	{
		string str1 = "str";
		string str2 = to_string(i);
		string str3 = ".txt";
		string str4 = str1 + str2 + str3;

		fstream  txtFile;
		txtFile.open(str4.c_str(), ios::out | ios::trunc);
		
		for(int j = 1; j <= atomCount; j++)
		{
			float r;
			txtFile << "B\t";
			for(int k = 0; k < 3; k++)
			{
				r = (float)rand()/(float)(RAND_MAX/(2*range)) - range;
				txtFile << r << "\t";
			}
			txtFile << "\n";
		}
		
	}
	return 0;
}