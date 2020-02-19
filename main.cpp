#include "eigen-git-mirror/Eigen/Dense"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ostream>

//Choose Starting microstate
Eigen::MatrixXd get_starting_microstate()  //In reality, this would need to be Xf
{
	//Initialize Matrix with all ones
	Eigen::MatrixXd starting_microstate = Eigen::MatrixXd::Ones(10,10); 
	return starting_microstate;
}


//Pick random site
class Site{
	public:
	
	int x_Index(int Nx);
	int y_Index(int Ny);
	
	private:
};

//member functions for site class
//get random x_index based off of Nx
int Site::x_Index(int Nx) {
	double random_number=rand()%Nx;
	return floor(random_number);
}

//get random y_index based off of Ny
int Site::y_Index(int Ny) {
	double random_number=rand()%Ny;
	return floor(random_number);
}

//Flip site on previous microstate
Eigen::MatrixXd flip_site(int Nx, int Ny, Eigen::MatrixXd previous_microstate) //change to MatrixXf
{
	Site site;
	auto x_index=site.x_Index(Nx);
	auto y_index=site.y_Index(Ny);
	previous_microstate(x_index, y_index)=-previous_microstate(x_index, y_index);
//	previous_microstate.coeffRef(5, 10)=-previous_microstate.coeffRef(5,10);
	//	for (int i=0; i<previous_microstate.rows(); i++)
//	{
//		for (j=0; j<previous_microstate.cols(); j++)
//		{
//			if (i==x_index) //might be remaking index with every 
//			{
//				if (j==y_index)
//				{
//					previous_microstate[i][j]=-previous_microstate[i][j]
	return previous_microstate;
}
/*
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
*/

int main()
{

   	srand(time(NULL));
	auto my_mat_of_ones=get_starting_microstate();
	std::cout<<my_mat_of_ones;
        auto flipped=flip_site(10, 10, my_mat_of_ones);
	std::cout<<"\n\n"<<flipped;
	auto flipped_2=flip_site(10,10,flipped);
	std::cout<<"\n\n"<<flipped_2;
	auto flipped_3=flip_site(10,10,flipped_2);
	std::cout<<"\n\n"<<flipped_3;
//	if (argc<2)
//	{
//		std::cout<<"You;ve made a terrible mistake. Enter text file."<<'\n';
//	        return 1;
//	}
//	std::string filename=argv[1];
}	
