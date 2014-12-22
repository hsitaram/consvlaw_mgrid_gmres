#include<iostream>
#include<fstream>
#include<cmath>

int main()
{
	int np;
	double x;
	double y;
	double dx;
	double t;
	double eta;
	double dcoeff;

	std::ofstream outfile("exact.dat");

	np=100;
	dx = 1.0/(double(np)-1);
	t = 0.01;
	dcoeff=1.0;

	for(int i=0;i<np;i++)
	{
		x=i*dx;
		eta=0.5*(1.0-x)/sqrt(dcoeff*t);
		y=1-erf(eta);
		outfile<<x<<"\t"<<y<<"\n";
	}

	outfile.close();
	return(0);
}
