#include <iostream>
#include <stdlib.h>
#include "PipPimNpi.h"

int main(int argc, char** argv)
{
	if ( argc<2 ){
		std::cout<<"Usage: p_resolution *.root"<<std::endl;
		return -1;
	}
	TFile *f = new TFile(argv[1]);
	if (!f){
		std::cout<<"Can not open file: "<<argv[1]<<std::endl;
		return -2;
	}
	TTree *tree;
	f->GetObject("Jpsi3pi",tree);
	if (!tree){
		std::cout<<"Can not get tree: Jpsi3pi"<<std::endl;
		return -3;
	}

	PipPimNpi alg(tree);
	alg.Loop();
	std::cout<<"done"<<std::endl;
	return 0;

}
