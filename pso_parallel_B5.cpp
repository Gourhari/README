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
#include <sys/stat.h>

#define atomCount 5  //For B5
//#define atomCount 5
//#define structCount 14
#define structCount 6
#define range 3
#define topKGBCount 10
#define iterMax 50
#define spin 2
#define topKGBdecimal 2


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

inline bool exists_test (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
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
	//char procShared[] = "%nprocshared=2";
	char procShared[] = "%nprocshared=4";

	// Commented by SSL on 04 12 2017 char memUsed[] = "%mem=2GB";
	char memUsed[] = "%mem=8GB";

	char noSave[] = "%nosave";
	char titleCard[] = "B5 struct"; //Changed by SSL from B5 to B6
	char method[] = "# b3lyp/6-311+g*";

	char chargeSpin[10];
	if(spin == 1)
	{
		string strTemp1 = "0 1";
		strcpy(chargeSpin, strTemp1.c_str());
	}
	else if (spin == 3)
	{
		string strTemp1 = "0 3";
		strcpy(chargeSpin, strTemp1.c_str());
	}

	gjfFile << procShared << endl;
	gjfFile << memUsed << endl;
	gjfFile << noSave << endl;
	gjfFile << method << endl;
	gjfFile << endl << titleCard << endl << endl;
	gjfFile << chargeSpin << endl;

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

//Changed slightly to include both singlet and triplet
string readEnergy(string logFile)
{
	string search;
	if(spin == 1)
	{
		search = "E(RB3LYP)";
	}
	else if (spin == 3)
	{
		search = "E(UB3LYP)";
	}
	
	ifstream inFile;
	string line;
	inFile.open(logFile.c_str());
	string zeroRet = "0.00";
	if(!inFile)
	{
		cout << "Unable to open file" << endl;
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
		while(!checkValidInd(str6))
		{
			;
		}
	}
return 1;
}

void initializePSO(struct particle particles[], int structureCount)
{
	string s1 = "strAll.txt";
	// Checks if there is an existing file with completed iteration and uses it if present
	if(exists_test(s1))
	{
		fstream txtFile;
		txtFile.open(s1.c_str(), ios::in);
		for(int i = 1; i <= structCount; i++)
		{
			char data[20];
			for(int j=0; j < atomCount; j++)
			{
				txtFile >> particles[i-1].atoms[j].symbol;
				strcpy(particles[i-1].pBest[j].symbol, particles[i-1].atoms[j].symbol);

				txtFile >> data;
				string temp1(data);
				particles[i-1].atoms[j].coX = strtod(temp1.c_str(), NULL);
				
				txtFile >> data;
				string temp2(data); 
				particles[i-1].atoms[j].coY = strtod(temp2.c_str(), NULL);
				
				txtFile >> data; 
				string temp3(data);
				particles[i-1].atoms[j].coZ = strtod(temp3.c_str(), NULL);
			}
			for(int j=0; j < atomCount; j++)
			{
				txtFile >> data;
				string temp4(data);
				particles[i-1].atoms[j].velX = strtod(temp4.c_str(), NULL);
				
				txtFile >> data;
				string temp5(data); 
				particles[i-1].atoms[j].velY = strtod(temp5.c_str(), NULL);
				
				txtFile >> data; 
				string temp6(data);
				particles[i-1].atoms[j].velZ = strtod(temp6.c_str(), NULL);
			}
			
			txtFile >> data;
			string temp10(data);
			particles[i-1].cost = strtod(temp10.c_str(), NULL);
			
			for(int j=0; j < atomCount; j++)
			{
				txtFile >> data;
				string temp7(data);
				particles[i-1].pBest[j].coX = strtod(temp7.c_str(), NULL);
				
				txtFile >> data;
				string temp8(data); 
				particles[i-1].pBest[j].coY = strtod(temp8.c_str(), NULL);
				
				txtFile >> data; 
				string temp9(data);
				particles[i-1].pBest[j].coZ = strtod(temp9.c_str(), NULL);
			}

			txtFile >> data;
			string temp11(data);
			particles[i-1].pBestCost = strtod(temp11.c_str(), NULL);				

			gengjf(particles, i, 0);

			string str1 = "str";
			string str2 = to_string(i);
			string iter = "_0";
			string str5 = ".gjf";
			
			string command = "g09 " + str1 + str2 + iter + str5 + " &";
			system(command.c_str());
		}	
		txtFile.close();	
	}

	//If no such file is present, it reads from the inital str1.txt,.......
	else
	{
		//cout << "---------\n\nDEBUG: Inside Initialize PSO fo no existing file.\n\n------";
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

				particles[i-1].atoms[j].velX = 0.1;
				particles[i-1].atoms[j].velY = 0.1;
				particles[i-1].atoms[j].velZ = 0.1;
			}
			
			particles[i-1].pBestCost = 100000;
			particles[i-1].cost = 100000;
			gengjf(particles, i, 0);
			string command = "g09 " + str1 + str2 + iter + str5 + " &";  // Uncomment the '&' to execute the program in the background. 
			system(command.c_str());
			
			valFile.close();	
		}
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
		cout << "Structure " << (str1 + str2) << ":\n";

		string str6 = ".log";
		double energyCost = atof(readEnergy(str1+str2+iter+str6).c_str());
		particles[i-1].cost = particles[i-1].pBestCost = energyCost;
		printf("Energy Read: %lf\n", particles[i-1].pBestCost);
	}

}

//// Added by Gour on 06/01/2018 ---- Comparing if two structures are similar
int isSimStruct(struct particle *particles, int a, int b)
{
	//cout << "Inside similar structure function for " << a << ", " << b << "\n";
	for(int i = 1; i <= atomCount; i++)
	{
		if(fabs(particles[a].atoms[i].coX - particles[b].atoms[i].coX) > 0.0000001)
			return 0;
		if(fabs(particles[a].atoms[i].coY - particles[b].atoms[i].coY) > 0.0000001)
			return 0;
		if(fabs(particles[a].atoms[i].coZ - particles[b].atoms[i].coZ) > 0.0000001)
			return 0;
	}
return 1;
}

int isSimEnergy(struct particle *particles, int a, int b)
{
	if(fabs(particles[a].cost - particles[b].cost) < 0.0001)
		return 1;
	else
		return 0;
}

//// Ending Added by Gour on 06/01/2018

//// Added by Gour on 06/01/2018 ----- Replacing the duplicate Structure with a random structure.
void randomStruct(struct particle *particles, int a)
{
	for(int k = 1; k <= atomCount; k++)
	{
		particles[a].atoms[k].coX = (float)rand()/(float)(RAND_MAX/(2*range)) - range;
		particles[a].atoms[k].coY = (float)rand()/(float)(RAND_MAX/(2*range)) - range;
		particles[a].atoms[k].coZ = (float)rand()/(float)(RAND_MAX/(2*range)) - range;
	}				
}
//// Ending Added by Gour on 06/01/2018


// Added by Gour on 10/01/2018 ------ Function used in sorting of top k global-bests
bool compareEnergy(struct particle p1, struct particle p2)
{
    return (p1.cost <= p2.cost);
}
// Ending Added by Gour on 10/01/2018


// Added by Gour on 18.01.2018 -- Function to save the positions and velocities at the end of each iteration
void saveStructures(struct particle *particles)
{
	//cout << "---------\n\nDEBUG: Inside save structures\n\n------";
	string save = "strAll.txt";
	fstream  txtFile;
	txtFile.open(save.c_str(), ios::out | ios::trunc);

	for(int i = 1; i <= structCount; i++)
	{
		for(int j = 0; j < atomCount; j++)
		{
			char buff[100];
			int n;

			txtFile << " " << particles[i-1].atoms[j].symbol << "\t\t\t\t\t";
			
			n = sprintf(buff, "%lf\t", particles[i-1].atoms[j].coX);
			txtFile << buff;
			
			n = sprintf(buff, "%lf\t", particles[i - 1].atoms[j].coY);
			txtFile << buff;
			
			n = sprintf(buff, "%lf", particles[i - 1].atoms[j].coZ);
			txtFile << buff << "\n";
			
		}

		txtFile << "\n";



		for(int j = 0; j < atomCount; j++)
		{
			char buff[100];
			int n;
			
			n = sprintf(buff, "%lf\t", particles[i-1].atoms[j].velX);
			txtFile << buff;
			
			n = sprintf(buff, "%lf\t", particles[i-1].atoms[j].velY);
			txtFile << buff;
			
			n = sprintf(buff, "%lf", particles[i-1].atoms[j].velZ);
			txtFile << buff << "\n";
			
		}

		txtFile << "\n";

		char buffTemp[100];
		int x;
		x = sprintf(buffTemp, "%lf", particles[i - 1].cost);
		txtFile << buffTemp << "\n";

		txtFile << "\n";

		for(int j = 0; j < atomCount; j++)
		{
			char buff[100];
			int n;
			
			n = sprintf(buff, "%lf\t", particles[i-1].pBest[j].coX);
			txtFile << buff;
			
			n = sprintf(buff, "%lf\t", particles[i-1].pBest[j].coY);
			txtFile << buff;
			
			n = sprintf(buff, "%lf", particles[i-1].pBest[j].coZ);
			txtFile << buff << "\n";
			
		}

		x = sprintf(buffTemp, "%lf", particles[i - 1].pBestCost);
		txtFile << buffTemp << "\n";

		txtFile << "\n";
		txtFile << "\n\n\n";
	}
	txtFile.close();
}


void savetopKGB(std::vector<struct particle> topKGB, int topKGBFilled)
{
	string topKGBFileString = "topKGB.txt";
	fstream topKGBFile;
	topKGBFile.open(topKGBFileString.c_str(), ios::out | ios::trunc);

	for(int i = topKGBFilled-1; i >= 0; i++)
	{
		char buff[100];
		int n;
		for(int j = 0; j < atomCount; j++)
		{
			topKGBFile << " " << topKGB[i].atoms[j].symbol << "\t\t\t\t\t";
			
			n = sprintf(buff, "%lf\t", topKGB[i].atoms[j].coX);
			topKGBFile << buff;
			
			n = sprintf(buff, "%lf\t", topKGB[i].atoms[j].coY);
			topKGBFile << buff;
			
			n = sprintf(buff, "%lf", topKGB[i].atoms[j].coZ);
			topKGBFile << buff << "\n";
		}
		n = sprintf(buff, "%lf\t", topKGB[i].cost);
		topKGBFile << "Energy of this structure: ";
		topKGBFile << buff;
		topKGBFile << " \n";
	}
}

int main(int argc, char const *argv[])
{
	struct particle particles[structCount];
	initializePSO(particles, structCount);
	struct particle globalBest;
	globalBest.cost = 1000000;
	// Commented by SSL on 03 12 2017 float wMax = 0.9;
	float wMax = 0.6;
	// Commented by SSL on 03.12.2017 float wMin = 0.4;
	float wMin = 0.2;
	float w;

	int iter = 1;
	
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
	
	vector <struct particle> topKGB;
	int topKGBFilled = 0;
	float compDec = 1.00;
	for(int idec = 1; idec <= topKGBdecimal; idec++)
		compDec /= 10;
	
	
	while(iter <= iterMax) 
	{
		// Added by Gour on 07.01.2018
		for(int i = 1; i <= structCount; i++)
		{
			for(int j = i+1; j <= structCount; j++)
			{
				if(isSimEnergy(particles, i, j) == 1)
				{
					cout << "Similar structures for " << i << ", " << j << "\n";
					randomStruct(particles, j);
					cout << "Replaced structure " << j << " with random structure.\n";
				}
			}
		}

		// Ending Added by Gour on 07.01.2018

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
		
		printf ("==================================\n");
		printf("\\Global Best So far at the end of Iteration: %d\n", iter);
		for(int j = 0; j < atomCount; j++)
		{
			printf("%lf \t %lf \t %lf\n", globalBest.atoms[j].coX, globalBest.atoms[j].coY, globalBest.atoms[j].coZ);
		}
	
		printf("Global Best Cost: %lf\n", globalBest.cost);
		printf ("==================================\n"); 

		//End of Added by SSL on 29/11/2017

		// Added by Gour on 06/01/2018 ---- Deleting gjf, log files of older iterations except the last 3
		if(iter >= 3)
		{
			for(int i = 1; i <= structCount; i++)
			{
				string str1 = "str";
				string str2 = to_string(i);
				string iterString = "_" + to_string(iter-3);
				string str5 = ".gjf";
				string str6 = ".log";

				string delgjf = str1 + str2 + iterString + str5;
				string dellog = str1 + str2 + iterString + str6;
				remove(delgjf.c_str());
				remove(dellog.c_str());
			}
		}
		// End of Added by Gour on 06/01/2018

		// top KGB
				
		if(topKGBFilled == 0)
		{
			//cout << "DEBUG: First push\n";
			topKGB.push_back(globalBest);
			topKGBFilled++;
		}
		else
		{
			//cout << "DEBUG: KGB Currently filling: " << topKGBFilled << ". " << globalBest.cost << " \n";
			//sort(topKGB.begin(), topKGB.end(), compareEnergy);
			//cout << "Diff: " << topKGB[topKGBFilled - 1] << "\n";// - enerVals[i] << "\n";
			if(topKGB[topKGBFilled - 1].cost - globalBest.cost > compDec){
				topKGB.push_back(globalBest);
				topKGBFilled++;
			}
			else
			{
				topKGB.pop_back();
				topKGB.push_back(globalBest);
			}
			if(topKGBFilled > topKGBCount)
			{
				//cout << "Condition Met: topKGBFilled > topKGBCount\n";
				topKGB.erase(topKGB.begin());
			}
		}
		//  Added by Gour -------Saving structures after every iteration---- on 18/01/2018
		//cout << "DEBUG: Saving topKGB\n";
		//savetopKGB(topKGB, topKGBFilled);

		//cout << "DEBUG: Done saving topKGB\n";
		saveStructures(particles);

		//Ending Added by Gour ------- on 18/01/2018
		iter++;
	}

		string topKGBFileString = "topKGB.txt";
		fstream topKGBFile;
		topKGBFile.open(topKGBFileString.c_str(), ios::out | ios::trunc);

		for(int i = topKGBFilled-1; i >= 0; i--)
		{
			char buff[100];
			int n;
			for(int j = 0; j < atomCount; j++)
			{
				topKGBFile << " " << topKGB[i].atoms[j].symbol << "\t\t\t\t\t";
				
				n = sprintf(buff, "%lf\t", topKGB[i].atoms[j].coX);
				topKGBFile << buff;
				
				n = sprintf(buff, "%lf\t", topKGB[i].atoms[j].coY);
				topKGBFile << buff;
				
				n = sprintf(buff, "%lf", topKGB[i].atoms[j].coZ);
				topKGBFile << buff << "\n";
			}
			n = sprintf(buff, "%lf\t", topKGB[i].cost);
			topKGBFile << "Energy of this structure: ";
			topKGBFile << buff;
			topKGBFile << " \n";
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