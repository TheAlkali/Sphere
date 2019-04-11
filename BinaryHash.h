#pragma once

#include "Common.h"
#include "Points.h"
#include <fstream>
#include <sstream>


// hamming distance function
__inline int Compute_HD(bitset<BCODE_LEN> &a, bitset<BCODE_LEN> &b)
{
	return ( ( a ^ b ).count() );
}

// spherical hamming distance function
__inline double Compute_SHD(bitset<BCODE_LEN> &a, bitset<BCODE_LEN> &b)
{
	return ( ( (double)( ( a ^ b ).count() ) ) / ( (double)( ( a & b ).count() ) + 0.1 ) );
}

class Index_Distance
{
public :
	bool flag;
	int index;
	REAL_TYPE dist, distSq;
	bool operator < (const Index_Distance &T) const
	{
		if( this->distSq < T.distSq )	{			return true;		}
		return false;
	}
};

class Sphere
{
public :
	REAL_TYPE *c, r, rSq;
	
	void Initialize(int _dim)
	{
		c = new REAL_TYPE [ _dim ];
		r = 0.0;		rSq = 0.0;
	}

	// function to set radius to include desired portion of training set
	void Set_Radius(Points *ps, Index_Distance *ids);
};

class SphericalHashing
{
public :
	Points *ps;

	// training set
	Points tps;

	Sphere *s;

	int code_len;
	int dim;

	Index_Distance **ids;
	bitset<NUM_TRAIN_SAMPLES> *table;

	void Initialize(Points *_p,int code_len,int dim);
	void Compute_Table();
	void Compute_Num_Overlaps(int **overlaps);
	void Set_Spheres();

	void ReleaseMem();

	template<typename ArgType>
	__inline void Compute_BCode(ArgType *x, bitset<BCODE_LEN> &y)
	{

		int start = KMER_SIZE;
	#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
	#endif
		for(int i=0;i<BCODE_LEN;i++)
		{
			ArgType dis = Compute_Distance_L2Sq<REAL_TYPE>( s[i].c , x , dim , start);
			if( dis - s[i].rSq >0.0001)
			{
				y[i] = 0;
			}
			else
			{
				y[i] = 1;
			}
			//std::cout << dis << " " << s[i].rSq << " " << y[i] << std::endl;
			start += dim;
		}
	//	std::cout << y << std::endl<< std::endl;
	}

	void Save_Sphere_Info()
	{
		std::ofstream file;
		file.open("bin/sphere_info.log");	
		for (int i = 0; i < code_len; ++i)
		{
			for (int j = 0; j < dim; ++j)
			{
				file << s[i].c[j] << " " ;
			//	std::cout << s[i].c[j] << "\t" ;
			}
			file << s[i].rSq << std::endl;
		//	std::cout << s[i].rSq << std::endl;
		}
		file.close();
	}

	void Load_Sphere_Info(int _code_len,int _seg_len)
	{
		code_len = _code_len;
		dim = _seg_len;

		std::ifstream file;
		file.open("bin/sphere_info.log");
		std::string tmp;

		s = new Sphere [ code_len ];
		for (int i = 0; i < code_len; ++i)
		{
			s[i].Initialize( dim );

			getline(file,tmp);
			std::stringstream ss(tmp);
			for (int j = 0; j < dim; ++j)
			{
				ss >> s[i].c[j];
			//	std::cout << s[i].c[j] << "\t";
			}
			ss >> s[i].rSq;
		//	std::cout << s[i].rSq << std::endl;
		}
		file.close();
	}
};
