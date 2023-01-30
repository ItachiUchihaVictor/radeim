#include<stdio.h>
#include<cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]){
	std::ifstream inputfile(string(argv[1]), std::ios::in);
	std::ofstream outputfile(string(argv[2]), std::ios::out);

	int numV = 0;
	int numE = 0;

	int min = 100000000;
	int max = 0;
	inputfile >> numV >> numE;

        for(int i=0;i<numE;i++){
		int x, y;
		inputfile >> x >> y;
		max = std::max(std::max(max, x), y);
		min = std::min(std::min(min, x), y);
	}
	inputfile.close();
	
	std::ifstream inputfileagain(string(argv[1]), std::ios::in);
	inputfileagain >> numV >> numE;
	
	numV = max - min + 1;
	outputfile << numV << " " << numE << endl;

        for(int i=0;i<numE;i++){
		int x, y;
		inputfileagain >> x >> y;
		outputfile << x-min << "\t" << y-min << endl;
	}
	return 0;
}
