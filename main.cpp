#include "eigen-git-mirror/Eigen/Dense"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ostream>

//Choose Starting microstate
Eigen::Matrix2f get_starting_microstate()  //In reality, this would need to be Xf
{
	Eigen::Matrix2f starting_microstate<<-1, 1, 
						1, 1;
	return starting_microstate;
}


//Pick random site
class Site{
	public:
	{
	int x_Index(int Nx)
	int y_Index(int Ny)
	}
	private:
}

//member functions for site class
//get random x_index based off of Nx
double Site::x_index(int Nx) {
   	random_number=rand()%1;
	return floor(Nx*random_number);
}

//get random y_index based off of Ny
double Site::y_index(int Nx) {
   	random_number=rand()%1;
	return floor(Nx*random_number);
}

//Flip site on previous microstate
Eigen::Matrix2f flip_site(Site::Site site, int Nx, int Ny, Eigen::Matrix2f previous_microstate) //change to MatrixXf
{
	x_index=site.x_Index(Nx);
	y_index=site.y_Index(Ny);
	previous_microstate[x_index][y_index]=-previous_microstate[x_index][y_index];
	//	for (int i=0; i<previous_microstate.rows(); i++)
//	{
//		for (j=0; j<previous_microstate.cols(); j++)
//		{
//			if (i==x_index) //might be remaking index with every 
//			{
//				if (j==y_index)
//				{
//					previous_microstate[i][j]=-previous_microstate[i][j]
}

vector<double> eci(Eigen::Matrix2f microstate)
{
	//Not sure what to put here either... Somehow must calculate these. Linear reg?
}


//double Ei(Eigen::Matrix2f microstate, vector<double> ECI) //change to MatrixXf when done testing
//{
	// I don't know what to put here
//}

//double Ef(Eigen::Matrix2f microstate, vector<double> ECI)
//{
	// I don't know what to put here
//}


//Calculate change in grand canonical energy
double change_in_energy(Ei, Ef, chemical_potential, Ni, Nf)
{
	//maybe remove Ei and Ef
	//double delta_E=-2*L[i,j]*(Vp+Vn*L[i,j](L[i+1, j]+L[i, j+1]+ L[i-1, j]+L[i, j-1])) 
	//somehow do this... remember modulo functions for edges/sides!!!!
}
//Accept new microstate


//Conditionally accept new microstate
Eigen::Matrix2f which_microstate(double delta_energy, double chemical_potential, int Ni, int Nf, Eigen::Matrix2f previous_microstate, Eigen::Matrix2f testing_microstate)
{
	double k= 
	Eigen::Matrix2f ending_microstate;
	double test_number=rand()%1;
	if ((delta_energy-chemical_potential*(Nf-Ni))<0)
	{
		ending_microstate=testing_microstate;
	}
	else if (exp(-(1/(k*T))*(delta_energy-chemical_potential*(Nf-Ni))=>test_number
	{
		ending_microstate=testing_microstate;
	}
	else{
		ending_microstate=previous_microstate;
	}
	return ending_microstate;
}


//keep old microstate (may not need separate function)


//Determine next microstate


//If at sides, do this!


//Get file with lattice and basis


int main(int argc, char *argv[])
{
	if (argc<2)
	{
		std::cout<<"You;ve made a terrible mistake. Enter text file."<<'\n';
	        return 1;
	}
	std::string filename=argv[1];
}	
