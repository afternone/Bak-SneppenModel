#include "bs.h"
#include <ctime>
#include <iostream>
#include <fstream>
int main(int argc, char** argv)
{
	time_t t_start = time(NULL) ;
	srand(time(NULL));
	int d=atoi(argv[1]);
	double f0=atof(argv[2]);
	double alpha=atof(argv[3]);
	int nei=atoi(argv[4]);
	int sco=atoi(argv[5]);
	int N=atoi(argv[6]);
	if(d==1) bs1d(f0,alpha,nei,sco,N,false);
	if(d==2) bs2d(f0,alpha,nei,sco,N,false);
	if(d==3) bs3d(f0,alpha,nei,sco,N,false);
	if(d==4) bs4d(f0,alpha,nei,sco,N,false);
	if(d==5) bs5d(f0,alpha,nei,sco,N,false);
	if(d==6) bs6d(f0,alpha,nei,sco,N,false);
	if(d==7) bs7d(f0,alpha,nei,sco,N,false);
	if(d==8) bs8d(f0,alpha,nei,sco,N,false);
	time_t t_end = time(NULL) ;
	std::ofstream fout("finished.txt",std::ios::app);
	fout<<std::string(argv[1])+"_"+std::string(argv[3])+"_"+std::string(argv[2])<<"\t"<<N<<"\t"<<difftime(t_end,t_start)<<" s"<<std::endl;
	fout.close();
	return 0;
}
