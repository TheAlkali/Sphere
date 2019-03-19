#pragma once
#include <iostream>
#include <fstream>

//#include "Common.h"
#include "Utils.hpp"
#include "MemoryMapped.h"

class Points
{
public :
	int nP;
	int dim;
	REAL_TYPE **d;
	std::ifstream srcfile;
	enum point_type
	{
		training,
		region,
		mapping
	};


	void Initialize_MemoryMapped(std::string filename,point_type type)
	{
		
		MemoryMapped reads_file(filename);
		for (int i = 0; i < nP; ++i)
		{
			for (int j = 0; j < dim; ++j)
			{
			//	std::cout << i * (dim + 1) + j << std::endl;
				d[i][j] = ictoi_table[reads_file[i * (dim + 1) + j]];
			}
		//	std::cout << std::endl;
		}
		reads_file.close();
	}

	void Initialize(int _nP, int _dim)
	{
		nP = _nP;
		dim = _dim;
		d = new REAL_TYPE * [ nP ];
		for(int i=0;i<nP;i++)
		{
			d[i] = new REAL_TYPE [ dim ];
		}
	}

	// this function is for read point set from file
	// format:
	// (num. points) (dimensionality)
	// v0 (floats, number of elements is equal to dimensionality)
	// v1
	// ...
	int Initialize_From_File()
	{
		std::string tmp;
		int i;
		for(i=0;i<nP && srcfile.peek() != EOF;i++)
		{
			std::getline(srcfile,tmp);
			for(int k=0;k<dim;k++)
			{
				d[i][k] = ictoi_table[tmp[k]];
			//	std::cout << d[i][k];
			}
		//	std::cout << std::endl;
		}
		if (srcfile.peek() == EOF)
		{
			srcfile.close();
			return i;
		}
		return -1;
	}

/*	void Initialize_From_Bcodes(bitset<BCODE_LEN> *bcode)
	{
		for (int i = 0; i < nP; ++i)
		{
			for (int j = 0; j < dim; ++j)
			{
				d[i][j] = bcode[i][j];
			}
		}
	}

	void Initialize_For_Hashlearning()
	{
		std::string read;
		for (int i = 0; i < nP; ++i)
		{
			int rand_loc = rand() % (nP - dim + 1);
			srcfile.seekg(rand_loc * (dim + 1),std::ios_base::beg);
			getline(srcfile,read);
			for (int k = 0; k < dim; ++k)
			{
				d[i][k] = (REAL_TYPE)( read[k] - 48);
			}
		}
	}*/

	// computing center of points for zero centering
	void Compute_Center(REAL_TYPE *center)
	{
		double *tCenter = new double [dim];
		SetVector_Val<double>( tCenter , dim , 0.0 );
		for(int i=0;i<nP;i++)
		{
			for(int k=0;k<dim;k++)
			{
				tCenter[k] += d[i][k];
			}
		}
		for(int k=0;k<dim;k++)
		{
			tCenter[k] /= (double)(nP);
			center[k] = (REAL_TYPE)( tCenter[k] );
		}
		delete [] tCenter;
	}

	void ReleaseMem()
	{
		for(int i=0;i<nP;i++)
		{
			delete [] d[i];
		}

		delete [] d;
	}

	
};
