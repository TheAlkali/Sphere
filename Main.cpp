#include "BinaryHash.h"
//#include "Evaluation.h"
#include "Mapping.hpp"
#include "fastxParser.hpp"
//#include "SAMwriter.hpp"
//#include "fastHashing.hpp"

#ifdef USE_PARALLELIZATION
#include <omp.h>
#include <iostream>
#include <sstream>
#include<unistd.h>
#endif

int main(int argc, char const *argv[])
//int main()
{	
	int gate = 1;//std::atoi(argv[1]);
//	system("rm tmp/*.log");
	Stopwatch T0("");
    T0.Reset();     T0.Start();

    Stopwatch T1("");
    T1.Reset();     T1.Start();

    int seg_len = (DIM - KMER_SIZE) / BCODE_LEN;

	Mapping map(DIM);

	if (gate == 0)
	{
	    Points ref_buff;
		ref_buff.srcfile.open(INPUT_REF_FILE_NAME,std::ifstream::in);
		if (!ref_buff.srcfile.is_open())
		{
			perror(INPUT_REF_FILE_NAME);
			exit(EXIT_FAILURE);
		}
		ref_buff.Initialize(NUM_TRAIN_SAMPLES * BCODE_LEN,seg_len);
		ref_buff.Initialize_MemoryMapped(INPUT_REF_FILE_NAME,Points::point_type::training);
		map.Learn_Spherical_Hashing(ref_buff,BCODE_LEN, seg_len);
		ref_buff.srcfile.close();

		map.Suffix_Array(seg_len);

		T1.Stop();
		printf("- Index Time Finished (%f seconds)\n\n\n\n",T1.GetTime() );
	}else if (gate == 1)
	{
    	int read_size = store_reads();
		map.Load_Spherical_Hashing(read_size,BCODE_LEN,seg_len);
		map.Load_SA(seg_len);

		Stopwatch T2("");
		T2.Reset();     T2.Start();
		std::cout <<"- Analysis of read region ..." << std::endl;
		map.Get_Read_Region();
		T2.Stop();
		std::cout << "- Analysis Finished(" << T2.GetTime() << " seconds)" << std::endl;

		map.Hash_Mapping_with_SA();
	//	map.Hash_Mapping_with_SA(PAIR_2);
		T0.Stop();
		printf("- Total Running Time (%f seconds)\n",T0.GetTime() );

		T0.Reset();     T0.Start();
	// 	SAMwriter sp;  
	//    sp.gen_SAM_file(map.rpro);
	    printf("- Generate SAM File Finished (%f seconds)\n",T0.GetTime());
	    T0.Stop();
		map.Load_Ref_Info();
		map.Output_Result(PAIR_1);
		map.Output_Result(PAIR_2);
		//map.Ref_Of_Read();
	}
		
	return 0;
}

//test the speed of mapping using hashcode and directly
//- per mapping time using hashcode(not including computing hashcode) (0.000001 seconds)
//- per mapping time using hashcode(including computing hashcode) (0.000105 seconds)
//- per mapping time directly (0.000003 seconds)
/*
void test_speed()
{
	Stopwatch T0("");
	bitset<BCODE_LEN> *qcode= new bitset<BCODE_LEN>[1000];
	bitset<BCODE_LEN> *dcode = new bitset<BCODE_LEN>[10000];
	T0.Reset();		T0.Start();
	for (int i = 0; i < 1000; ++i)
	{
		src_sh.Compute_BCode(read_buff_1.d[i],qcode[i]);
	}
	
	for (int i = 0; i < 10000; ++i)
	{
		src_sh.Compute_BCode(ref_buff.d[i],dcode[i]);
	}
	T0.Stop();
	printf("- computing hashcode (%f seconds)\n",T0.GetTime() );
	T0.Reset();		T0.Start();
	for (int j = 0; j < 1000; ++j)
	{
		for (int i = 0; i < 10000; ++i)
		{
			//src_sh.Compute_BCode(ref_buff.d[i],dcode[i]);
			Compute_HD(qcode[j],dcode[i]);
		}
	}
	
	T0.Stop();
	printf("- per mapping time using hashcode (%f seconds)\n",T0.GetTime() );

	T0.Reset();		T0.Start();
	for (int j = 0; j < 1000; ++j)
	{
		for (int i = 0; i < 10000; ++i)
		{
			Compute_Distance_L2Sq<REAL_TYPE>(read_buff_1.d[0],ref_buff.d[i],DIM);
		}
	}
	T0.Stop();
	printf("- per mapping time directly (%f seconds)\n",T0.GetTime() );
}*/