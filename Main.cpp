#include "BinaryHash.h"
//#include "Evaluation.h"
#include "Mapping.hpp"
#include "fastxParser.hpp"
#include "SAMwriter.hpp"
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
	int gate = std::atoi(argv[1]);
//	system("rm tmp/*.log");
	Stopwatch T0("");
    T0.Reset();     T0.Start();

    int seg_len = (DIM - KMER_SIZE) / BCODE_LEN;


	if (gate == 0)
	{
		Mapping map(DIM);
	    Points ref_buff;
		ref_buff.srcfile.open(INPUT_REF_FILE_NAME,std::ifstream::in);
		if (!ref_buff.srcfile.is_open())
		{
			perror(INPUT_REF_FILE_NAME);
			exit(EXIT_FAILURE);
		}
		ref_buff.Initialize(NUM_TRAIN_SAMPLES * BCODE_LEN,seg_len);
		ref_buff.Initialize_From_File();
		map.Learn_Spherical_Hashing(ref_buff,BCODE_LEN, seg_len);
		ref_buff.srcfile.close();

		map.Suffix_Array(seg_len);

		T0.Stop();
		printf("- Index Time Finished (%f seconds)\n\n\n\n",T0.GetTime() );
	}else if (gate == 1)
	{
		Mapping *map = new Mapping(DIM);
    	int read_size = store_reads();
		map->Load_Spherical_Hashing(read_size,BCODE_LEN,seg_len);
		map->Load_SA(seg_len);
		std::cout << "- filter:" << filter << std::endl;
		Stopwatch T2("");
		T2.Reset();     T2.Start();
		std::cout <<"- Analysis of read region ..." << std::endl;
		map->Get_Read_Region();
		T2.Stop();
		std::cout << "- Analysis Finished(" << T2.GetTime() << " seconds)" << std::endl;

		map->Hash_Mapping_with_SA();

	/*	map->Load_Ref_Info();
		map->Output_Result(PAIR_1);
		map->Output_Result(PAIR_2);*/
		T0.Stop();

	 	SAMwriter *sp = new SAMwriter(DIM); 
	 	sp->Transfer_Info_From_Mapping(map->is_read_1_rev,map->is_read_2_rev,map->rpro,map->read_1_buff.d,map->read_2_buff.d);
	 	delete map;

	 	sp->Load_Info();
	    sp->Analyse_Result_Pair();
	    sp->Merge_Result();
		T2.Reset();     T2.Start();
	    sp->Generate_SAM();
	    T2.Stop();
	    printf("- Generate SAM File Finished (%f seconds)\n",T2.GetTime());
		printf("- Total Running Time (%f seconds)\n",T0.GetTime() );
	}
		
	return 0;
}