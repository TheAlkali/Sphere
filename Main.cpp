#include "BinaryHash.h"
#include "Evaluation.h"
#include "Mapping.h"

#ifdef USE_PARALLELIZATION
#include <omp.h>
#include <iostream>
#include <sstream>
#include<unistd.h>
#endif

// ref_buff: data points set
// read_buff_1: query points set
Points ref_buff, read_buff_1,read_buff_2;

//Points train_bcode_P;

SphericalHashing src_sh;

int Initialize_Ref_Data(size_t buff_size)
{
#ifdef INPUT_REF_FILE_NAME
	ref_buff.Initialize(buff_size,DIM);
	int file_end = ref_buff.Initialize_From_File();
#endif			
	//if file end, return the count of unvisited data
	//if file not end, return -1
	return file_end;
}

void Learn_Spherical_Hashing(SphericalHashing &sh,Points &buff,int code_len){
	Stopwatch T0("");
	sh.Initialize(&buff,code_len);
	T0.Reset();		T0.Start();
	sh.Set_Spheres();
	T0.Stop();
	printf("- Learning Spherical Hashing Finished (%f seconds)\n",T0.GetTime());
//	sh.Save_Sphere_Info();
}

// TODO pair end read result output
void output_result(std::string readfile,std::string locfile,std::string disfile,std::string resfile){

	std::ofstream output(resfile);
	std::ifstream read(readfile);
	std::ifstream tmp_loc(locfile);
	std::ifstream tmp_dis(disfile);
//	std::ifstream refcode_file;

	int dim = read_buff_1.dim;
	size_t ref_loc = 0;
//	unsigned long ref_code;
	std::string tmp_ref,tmp_read; 
	// load mapping result from disk using cereal
	mapped_res mres_1,mres_2;
	{
		cereal::BinaryInputArchive ar_loc(tmp_loc);
		cereal::BinaryInputArchive ar_dis(tmp_dis);
		ar_loc(CEREAL_NVP(mres_1.mapped_ref_loc));
		ar_dis(CEREAL_NVP(mres_1.min_dis));
	};

	Stopwatch T0("");
	T0.Reset();		T0.Start();
	std::cout << mres_1.mapped_ref_loc.size() << std::endl;
	for (unsigned int qIndex = 0;qIndex < mres_1.mapped_ref_loc.size();++qIndex)
//	for (unsigned int qIndex = 0;qIndex < 1;++qIndex)
	{
		// output results  that the distance is  1
		getline(read,tmp_read);
	//	if (mres_1.min_dis[qIndex] <= 20)
		{	
			output << '>' << qIndex + 1  << ":" << mres_1.min_dis[qIndex] << ":" << std::endl;
		//	output << bCodeRead_SH[qIndex] << std::endl;
			output << "= " ;//<< tmp_read << std::endl;

			for (unsigned int i = 0; i < tmp_read.length(); ++i)
			{
				output << (char)itos_table[(int8_t)tmp_read[i]];
			//	std::cout << (char)itos_table[(int8_t)tmp_read[i]];
			}
			output << std::endl;
		//	std::cout << std::endl;
			
		//	std::cout << mres_1.mapped_ref_loc[qIndex].size() << std::endl;
			for (unsigned int i = 0; i < mres_1.mapped_ref_loc[qIndex].size(); ++i)
			{
				ref_loc = mres_1.mapped_ref_loc[qIndex][i];
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
	read.close();
	output.close();
//	refcode_file.close();
	T0.Stop();
	printf("- Save Mapping Results To Disk Finished (%f seconds)\n",T0.GetTime() );
}

int main()
{	
//	system("rm tmp/*.log");
	Stopwatch T0("");
    T0.Reset();     T0.Start();
//	Suffix_Array();
    Load_SA();

	srand( (unsigned int)( time(NULL) ) );

	ref_buff.srcfile.open(INPUT_REF_FILE_NAME,std::ifstream::in);
	read_buff_1.srcfile.open(INPUT_READ_FILE_NAME_1,std::ifstream::in);
	read_buff_2.srcfile.open(INPUT_READ_FILE_NAME_2,std::ifstream::in);

	//ref_buff.codefile.open(OUTPUT_REF_HASH_FILE_NAME,std::ios::out);//|std::ios::binary);
	//read_buff_1.codefile.open(OUTPUT_READ_HASH_FILE_NAME,std::ios::out);//|std::ios::binary);
	if (!ref_buff.srcfile.is_open())
	{
		perror(INPUT_REF_FILE_NAME);
		exit(EXIT_FAILURE);
	}
	if (!read_buff_1.srcfile.is_open())
	{
		perror(INPUT_READ_FILE_NAME_1);
		exit(EXIT_FAILURE);
	}
	if (!read_buff_2.srcfile.is_open())
	{
		perror(INPUT_READ_FILE_NAME_2);
		exit(EXIT_FAILURE);
	}

	ref_buff.Initialize(NUM_TRAIN_SAMPLES,DIM);
	//ref_buff.Initialize_For_Hashlearning();
	Initialize_Ref_Data(NUM_TRAIN_SAMPLES);
	Learn_Spherical_Hashing(src_sh,ref_buff,BCODE_LEN);

//	Initialize_Bcode_Data();
//	Learn_Spherical_Hashing(bcode_sh,train_bcode_P,BBCODE_LEN);

	ref_buff.srcfile.close();
	ref_buff.srcfile.open(INPUT_REF_FILE_NAME,std::ifstream::in);

//	Hash_Mapping();
	Hash_Mapping_with_SA(src_sh,read_buff_1,PAIR_1);
//	Hash_Mapping_with_SA(src_sh,read_buff_2,DIM,PAIR_2);
	T0.Stop();
	printf("- Total Running Time (%f seconds)\n",T0.GetTime() );

	output_result(INPUT_READ_FILE_NAME_1,PAIR_1_LOC_FILE,PAIR_1_DIS_FILE,PAIR_1_RES_FILE);
//	output_result(INPUT_READ_FILE_NAME_2,PAIR_2_LOC_FILE,PAIR_2_DIS_FILE,PAIR_2_RES_FILE);
		
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