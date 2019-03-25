#include "BinaryHash.h"
//#include "Evaluation.h"
#include "Mapping.hpp"
#include "fastxParser.hpp"
//#include "fastHashing.hpp"

#ifdef USE_PARALLELIZATION
#include <omp.h>
#include <iostream>
#include <sstream>
#include<unistd.h>
#endif

// TODO pair end read result output
void output_result(std::string readfile,std::string locfile,std::string disfile,std::string resfile,std::string ref_string,std::string code_file)
{

	std::ofstream output(resfile);
	std::ifstream read(readfile);
	std::ifstream tmp_loc(locfile);
	std::ifstream tmp_dis(disfile);
	std::ifstream tmp_code(code_file);
	std::ifstream refcode_file;

	int dim = DIM;
	size_t ref_loc = 0;
//	unsigned long ref_code;
	std::string tmp_read; 
	// load mapping result from disk using cereal
	mapped_res mres;
	{
		cereal::BinaryInputArchive ar_loc(tmp_loc);
		cereal::BinaryInputArchive ar_dis(tmp_dis);
		ar_loc(CEREAL_NVP(mres.mapped_ref_loc));
		ar_dis(CEREAL_NVP(mres.min_dis));
	};

	std::vector<bitset<BCODE_LEN>> read_code;
	{
		cereal::BinaryInputArchive ar_code(tmp_code);
		ar_code(CEREAL_NVP(read_code));
	}

	Stopwatch T0("");
	T0.Reset();		T0.Start();
	std::cout << mres.mapped_ref_loc.size() << std::endl;
	for (unsigned int qIndex = 0;qIndex < mres.mapped_ref_loc.size();++qIndex)
//	for (unsigned int qIndex = 0;qIndex < 1;++qIndex)
	{
		// output results  that the distance is  1
		getline(read,tmp_read);
	//	if (mres.min_dis[qIndex] <= 20)
		{	
			output << '>' << qIndex + 1  << ":" << mres.min_dis[qIndex] << ":" ;//<< std::endl;
			output << read_code[qIndex] << std::endl;
			output << "= " ;//<< tmp_read << std::endl;

			for (unsigned int i = 0; i < DIM; ++i)
			{
				output << (char)itos_table[(int8_t)tmp_read[i]];
			//	std::cout << (char)itos_table[(int8_t)tmp_read[i]];
			}
			output << std::endl;
		//	std::cout << std::endl;
			
		//	std::cout << mres.mapped_ref_loc[qIndex].size() << std::endl;
			for (unsigned int i = 0; i < mres.mapped_ref_loc[qIndex].size(); ++i)
			{
				ref_loc = mres.mapped_ref_loc[qIndex][i];
				if (ref_loc > 0)
				{
					//------ seek with sa------
					char base;
					output << "+ ";	
		//			std::cout << "+ ";	
					for (int j = 0; j < dim; ++j)
					{
						base = ref_string[ref_loc + j];
						output << (char)itos_table[(int8_t)base];
		//				std::cout << (char)itos_table[(int8_t)base];
					}
		//			std::cout << std::endl;
					output << " " << ref_loc << std::endl;
				}
			}
		}
		
	}
	read.close();
	output.close();
//	refcode_file.close();
	T0.Stop();
	printf("- Save Mapping Results To Disk Finished (%f seconds)\n",T0.GetTime() );
}

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
		map.Get_Read_Region(PAIR_1);
		map.Get_Read_Region(PAIR_2);
		T2.Stop();
		std::cout << "- Analysis Finished(" << T2.GetTime() << " seconds)" << std::endl;

		map.Hash_Mapping_with_SA(PAIR_1);
		map.Hash_Mapping_with_SA(PAIR_2);
		T0.Stop();
		printf("- Total Running Time (%f seconds)\n",T0.GetTime() );

		T0.Reset();     T0.Start();
	 	SAMwriter sp;  
	//    sp.gen_SAM_file();
	    printf("- Generate SAM File Finished (%f seconds)\n",T0.GetTime());
	    T0.Stop();
		
		output_result(INPUT_READ_FILE_NAME_1,PAIR_1_LOC_FILE,PAIR_1_DIS_FILE,PAIR_1_RES_FILE,map.ref_string,"tmp/read_code_1.bin");
		output_result(INPUT_READ_FILE_NAME_2,PAIR_2_LOC_FILE,PAIR_2_DIS_FILE,PAIR_2_RES_FILE,map.ref_string,"tmp/read_code_2.bin");
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