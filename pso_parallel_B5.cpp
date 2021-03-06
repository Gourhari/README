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

#define atomCount 5
#define structCount 14


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



void gengjf(struct particle particles[], int particleNum, int iteration)
{
	string str1 = "str";
	string str2 = to_string(particleNum);
	string str3 = ".gjf";
	string iter = "_" + to_string(iteration);
	string str4 = str1 + str2 + iter + str3;

	fstream  gjfFile;
	gjfFile.open(str4.c_str(), ios::out | ios::trunc);
	char procShared[] = "%nprocshared=2";
	char memUsed[] = "%mem=2GB";
	char noSave[] = "%nosave";
	char titleCard[] = "B5 struct";
	char method[] = "# b3lyp/6-311+g(d,p)";
	char chargeSpin[] = "0 2";

	gjfFile << procShared << endl;
	gjfFile << memUsed << endl;
	gjfFile << noSave << endl;
	gjfFile << method << endl;
	gjfFile << endl << titleCard << endl << endl;
	gjfFile << chargeSpin << endl;

	/*
	double dbl = 1278.65626487, dbl2 = 123.31541458;
	char buffer [100];
  	int n, a=5, b=3;
  	n=snprintf (buffer, 100, "%lf plus %lf is %lf", dbl, dbl2, dbl+dbl2);
  	printf ("[%s] is a string %d chars long\n\n",buffer,n);*/

	// Commentd by SSL on 01.12.2017 cout << "Generating gjf file: " << str4 << endl;
	for(int i = 0; i < atomCount; i++)
	{
		char buff[100];
		int n;
		const char* zeroFloat = "0.00000000";
		gjfFile << " " << particles[particleNum-1].atoms[i].symbol << "\t\t\t\t\t";
		if(particles[particleNum-1].atoms[i].coX == 0.0)
			gjfFile << zeroFloat << "\t";
		else
		{
			n = sprintf(buff, "%lf\t", particles[particleNum - 1].atoms[i].coX);
			gjfFile << buff;
		}
		if(particles[particleNum-1].atoms[i].coY == 0.0)
			gjfFile << zeroFloat << "\t";
		else
		{
			n = sprintf(buff, "%lf\t", particles[particleNum - 1].atoms[i].coY);
			gjfFile << buff;
		}
		if(particles[particleNum-1].atoms[i].coZ == 0.0)
			gjfFile << zeroFloat << "\n";
		else
		{
			n = sprintf(buff, "%lf", particles[particleNum - 1].atoms[i].coZ);
			gjfFile << buff << "\n";
		}
	}

	gjfFile << "\n\n\n\n";
	gjfFile.close();
		
}

string readEnergy(string logFile)
{
	string search = "E(UB3LYP)";
	ifstream inFile;
	string line;
	inFile.open(logFile.c_str());
	string zeroRet = "0.00";
	if(!inFile)
	{
		cout << "Unable to open file" << endl;
		//return zeroRet.c_str();
		return zeroRet;
	}
	string energyLine;
	int count = 0, lastFound = 0;
	size_t pos;
	while(inFile.good())
	{
		count++;
		getline(inFile,line); 
		pos=line.find(search); 
		if(pos!=string::npos) 
		{
			lastFound = count;
			energyLine = line;
		}
	}

	
	std::string delimiter = " ";
	size_t posi = 0;
	string token;
	string energy;
	int tokenCount = 0;
	while ((posi = energyLine.find(delimiter)) != std::string::npos) 
	{
		tokenCount++;
    	token = energyLine.substr(0, posi);
	    if(tokenCount == 8)
	    {
	    	energy = token;
	    	break;
	    }
	    energyLine.erase(0, posi + delimiter.length());
	}
	
//return energy.c_str();
return energy;
}

int checkValidInd(string logFile)
{
	string searchSuccess = "Normal termination";
	string searchError = "Error termination";
	ifstream inFile;
	string line;
	inFile.open(logFile.c_str());
	
	if(!inFile)
	{
		
		return 0;
	}

	string energyLine;
	size_t posSucc, posErr;
	while(inFile.good())
	{
		getline(inFile,line); 
		posSucc = line.find(searchSuccess);
		posErr =  line.find(searchError);
		if(posSucc != string::npos)
			return 1;
		if(posErr != string::npos)
			return 1;
	}

	return 0;
}

int checkValid(int iter)
{
	string str1 = "str";
	string str3 = ".txt";
	string str4 = ".gjf";
	string str5 = ".log";
	string iterString = to_string(iter);
	for(int i=1; i <= structCount; i++)
	{
		string str2 = to_string(i);
		string iter = "_" + iterString;
		string str6 = str1 + str2 + iter + str5;
		//cout << "Checking Valid Ind: " << str6 << "\n";
		while(!checkValidInd(str6))
		{
			;
		}
	}
return 1;
}

void initializePos(struct particle particles[], int structureCount)
{
	for(int i = 1;i < structureCount+1;i++)
	{
		string str1 = "str";
		string str2 = to_string(i);
		string iter = "_0";
		string str3 = ".txt";
		string str5 = ".gjf";
		string str4 = str1 + str2 + str3;
		// Commentd by SSL on 01.12.2017 cout << "FILE to be opened: " << str4 << "\n";
		fstream valFile;

		valFile.open(str4.c_str(), ios::in);
		char data[20];
		for(int j=0; j <atomCount; j++)
		{
			valFile >> particles[i-1].atoms[j].symbol;
			strcpy(particles[i-1].pBest[j].symbol, particles[i-1].atoms[j].symbol);
			
			valFile >> data;
			string temp1(data);
			particles[i-1].atoms[j].coX = strtod(temp1.c_str(), NULL);
			particles[i-1].pBest[j].coX = particles[i-1].atoms[j].coX;
			
			valFile >> data;
			string temp2(data); 
			particles[i-1].atoms[j].coY = strtod(temp2.c_str(), NULL);
			particles[i-1].pBest[j].coY = particles[i-1].atoms[j].coY;
			
			valFile >> data; 
			string temp3(data);
			particles[i-1].atoms[j].coZ = strtod(temp3.c_str(), NULL);
			particles[i-1].pBest[j].coZ = particles[i-1].atoms[j].coZ;

			particles[i-1].atoms[j].velX = 0;
			particles[i-1].atoms[j].velY = 0;
			particles[i-1].atoms[j].velZ = 0;
		}
		
		particles[i-1].pBestCost = 100000;
		particles[i-1].cost = 100000;
		gengjf(particles, i, 0);
		string command = "g09 " + str1 + str2 + iter + str5 + " &";  // Uncomment the '&' to execute the program in the background. 
		// Commentd by SSL on 01.12.2017 cout << "command: " << command <<"\n";
		system(command.c_str());
		
		valFile.close();	
	}

	cout << "Waiting for initial log files to be generated.\n";
	while(!checkValid(0));
	for(int i = 1; i <= structureCount; i++ )
	{
		
		string str1 = "str";
		string str2 = to_string(i);
		string iter = "_0";
		string str3 = ".txt";
		string str5 = ".gjf";
		string str4 = str1 + str2 + str3;
		// Commentd by SSL on 01.12.2017 cout << "FILE to be opened: " << str4 << "\n";
		
		string str6 = ".log";
		double energyCost = atof(readEnergy(str1+str2+iter+str6).c_str());
		particles[i-1].cost = particles[i-1].pBestCost = energyCost;
		// Commentd by SSL on 01.12.2017 printf("Energy Read: %lf\n", particles[i-1].pBestCost);
	}

}

/* Commented by SSL on 01.12.2017
int relDiff(double prev, double curr)
{
	if(prev-curr <= 0)
		return 0;
	else
	{
		double temp = (prev-curr)*100/prev;
		if(temp >= 0.01)
			return 1;
		else
			return 0;
	}
	

} End of commented by SSL on 01 12 2017 */

int relDiff(double prev, double curr)
{
	if(prev-curr <= 0)
		return 0;
	else
	{
		double temp = (prev-curr)*100/prev;
		//if(temp >= 0.01) Commented by SSL
		if(temp >= 0.0000001)
			return 1;
		else
			return 0;
	}
	

}

int main(int argc, char const *argv[])
{
	
	struct particle particles[structCount];
	initializePos(particles, structCount);
	struct particle globalBest;
	globalBest.cost = 1000000;
	float wMax = 0.9;
	float wMin = 0.4;
	float w;
	int iterMax = 2000, iter = 1;
	int c1 = 2, c2 = 2;
	double r1, r2;
	for(int i = 1; i < structCount+1 ;i++)
	{
		if(globalBest.cost > particles[i-1].cost)
		{
			globalBest.cost = particles[i-1].cost;
			for(int j = 0; j < atomCount; j++)
			{
				strcpy(globalBest.atoms[j].symbol, particles[i-1].atoms[j].symbol);
				globalBest.atoms[j].coX = particles[i-1].atoms[j].coX;
				globalBest.atoms[j].coY = particles[i-1].atoms[j].coY;
				globalBest.atoms[j].coZ = particles[i-1].atoms[j].coZ;
			}
		}
	}
	
	srand (static_cast <unsigned> (time(0)));
	double prevBest = 10000;

	

	while(iter < iterMax) // && relDiff(prevBest, globalBest.cost) == 1) //Commentd the last part by SSL
	//while(iter < iterMax && relDiff(prevBest, globalBest.cost) == 1) //Commented by SSL
	{

		prevBest = globalBest.cost;
		w = wMax - ((wMax - wMin) * iter / iterMax);
		// Commentd by SSL on 01.12.2017 cout << "\n\nITERATION: " <<  iter << " \n\n";
		for(int i = 1; i < structCount+1; i++)
		{
			// Commentd by SSL on 01.12.2017 printf("Particle %d: \n", i);
			for(int j = 0; j < atomCount; j++)
			{
			
				float r;
				r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				r1 = r;
				r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				r2 = r;
				//r1 = r2 = 1;
				particles[i-1].atoms[j].velX = w * (particles[i-1].atoms[j].velX) + c1 * r1 * (particles[i-1].pBest[j].coX - particles[i-1].atoms[j].coX) + c2 * r2 * (globalBest.atoms[j].coX - particles[i-1].atoms[j].coX);	
				particles[i-1].atoms[j].velY = w * (particles[i-1].atoms[j].velY) + c1 * r1 * (particles[i-1].pBest[j].coY - particles[i-1].atoms[j].coY) + c2 * r2 * (globalBest.atoms[j].coY - particles[i-1].atoms[j].coY);	
				particles[i-1].atoms[j].velZ = w * (particles[i-1].atoms[j].velZ) + c1 * r1 * (particles[i-1].pBest[j].coZ - particles[i-1].atoms[j].coZ) + c2 * r2 * (globalBest.atoms[j].coZ - particles[i-1].atoms[j].coZ);	
				particles[i-1].atoms[j].coX = particles[i-1].atoms[j].coX + particles[i-1].atoms[j].velX;
				particles[i-1].atoms[j].coY = particles[i-1].atoms[j].coY + particles[i-1].atoms[j].velY;
				particles[i-1].atoms[j].coZ = particles[i-1].atoms[j].coZ + particles[i-1].atoms[j].velZ;
			}
			/* Commenetd by SSL on 01.12.2017
			printf("Position\n");
			for(int j = 0; j < atomCount; j++)
			{
				printf("%lf \t %lf \t %lf\n", particles[i-1].atoms[j].coX, particles[i-1].atoms[j].coY, particles[i-1].atoms[j].coZ);
			}
			printf("Velocities\n");
			for(int j = 0; j < atomCount; j++)
			{
				printf("%lf \t %lf \t %lf\n", particles[i-1].atoms[j].velX, particles[i-1].atoms[j].velY, particles[i-1].atoms[j].velZ);
			} End of Commented by SSL on 01.12.2017 */

			
			string str1 = "str";
			string str2 = to_string(i);
			string iterString = "_" + to_string(iter);
			string str3 = ".txt";
			string str5 = ".gjf";
			gengjf(particles, i, iter);
			string command = "g09 " + str1 + str2 + iterString + str5 + " &"; // Uncomment the '&' to execute the program in the background. 
			// Commentd by SSL on 01.12.2017 cout << command <<"\n";
			system(command.c_str()); 
		}
		cout << "Waiting for ITERATION:" << iter <<" log files to be generated.\n";
		while(!checkValid(iter));
		for(int i = 1; i <= structCount; i++)
		{
			string str1 = "str";
			string str2 = to_string(i);
			string iterString = "_" + to_string(iter);
			string str3 = ".txt";
			string str5 = ".gjf";
			string str6 = ".log";
			double energyCost = atof(readEnergy(str1+str2+iterString+str6).c_str());
			particles[i-1].cost = energyCost;
			// Commentd by SSL on 01.12.2017 printf("Enery read: %lf\n", energyCost);
			if(particles[i-1].cost < particles[i-1].pBestCost)
			{
				for(int j = 0; j < atomCount; j++)
				{
					particles[i-1].pBest[j].coX = particles[i-1].atoms[j].coX;
					particles[i-1].pBest[j].coY = particles[i-1].atoms[j].coY;
					particles[i-1].pBest[j].coZ = particles[i-1].atoms[j].coZ;
				}
				particles[i-1].pBestCost = particles[i-1].cost;
			}
			if(particles[i-1].cost < globalBest.cost)
			{
				for(int j = 0; j < atomCount; j++)
				{
					globalBest.atoms[j].coX = particles[i-1].atoms[j].coX;
					globalBest.atoms[j].coY = particles[i-1].atoms[j].coY;
					globalBest.atoms[j].coZ = particles[i-1].atoms[j].coZ;
				}
				globalBest.cost = particles[i-1].cost;	
			}
		}


		//Added by SSL on 29/11/2017
		/*
		printf ("==================================\n");
		printf("\\Global Best So far at the end of Iteration: %d\n", iter);
		for(int j = 0; j < atomCount; j++)
		{
	printf("%lf \t %lf \t %lf\n", globalBest.atoms[j].coX, globalBest.atoms[j].coY, globalBest.atoms[j].coZ);
		}
	
	printf("Global Best Cost: %lf\n", globalBest.cost);
		printf ("==================================\n"); */

		//End of Added by SSL on 29/11/2017



		iter++;
	}
	
		/* Added by SSL on 01.12.2017 */
		printf ("Execution Complete \n");
		printf ("==================================\n");
		printf ("Personal Bests:\n");

		for(int i = 0; i < structCount; i++)
		{
		  printf ("Particle No.: %d\n", i);
   		  for(int j = 0; j < atomCount; j++)
		     printf("%lf \t %lf \t %lf\n", particles[i].pBest[j].coX, particles[i].pBest[j].coY, particles[i].pBest[j].coZ);
	          printf ("Best Personal Cost: %lf\n",particles[i].pBestCost);
		  printf ("\n");
		}
			
			
		printf ("==================================\n");

		/* End of Added by SSL on 01.12.2017 */



		printf("Final Position\n");
		for(int j = 0; j < atomCount; j++)
		{
			printf("%lf \t %lf \t %lf\n", globalBest.atoms[j].coX, globalBest.atoms[j].coY, globalBest.atoms[j].coZ);
		}
	
	printf("Global Best Cost: %lf\n", globalBest.cost);
	return 0;
}