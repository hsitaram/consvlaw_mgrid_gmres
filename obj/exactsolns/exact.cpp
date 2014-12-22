#include<iostream>
#include<fstream>
#include<cmath>

int main()
{
	int np;
	double x;
	double y;
	double dx;
	double e;

	std::ofstream outfile("exact.dat");

	np=100;
	dx = 1.0/(double(np)-1);
	e = exp(1);

	for(int i=0;i<np;i++)
	{
		x=i*dx;
		//y=e/(e*e-1)*(exp(-x)-exp(x))+x;
		y=x*x-x;
		outfile<<x<<"\t"<<y<<"\n";
	}

	outfile.close();
	return(0);
}
