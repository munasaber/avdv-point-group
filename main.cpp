#include "eigen-git-mirror/Eigen/Dense"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ostream>

//Choose Starting microstate


//Flip site on previous microstate


//Calculate change in grand canonical energy


//Accept new microstate


//Conditionally accept new microstate


//keep old microstate (may not need separate function)


//Determine next microstate


//If at sides, do this!




int main(int argc, char *argv[])
{
	if (argc<2)
	{
		std::cout<<"You;ve made a terrible mistake. Enter text file."<<'\n';
	        return 1;
	}
	std::string filename=argv[1];
}	
