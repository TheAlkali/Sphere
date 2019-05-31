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
	std::string gate = "mapping";argv[1];//std::atoi(argv[1]);
	Stopwatch T0("");
    T0.Reset();     T0.Start();

    int seg_len = (DIM - KMER_SIZE - SKIP - IGNORE) - BCODE_LEN + 1; /// BCODE_LEN;
    int dim  = DIM;
    region_profile rpro;

	if (gate.compare("index") == 0)
	{
		Mapping *map = new Mapping(dim,rpro);
		std::cout << "- dim:" << DIM << std::endl;
		std::cout << "- segment length:" << seg_len << std::endl;
	    Points ref_buff;
		
		ref_buff.Initialize(NUM_TRAIN_SAMPLES * BCODE_LEN,seg_len);
		ref_buff.Initialize_From_File(INPUT_REF_FILE_NAME);
		map->Learn_Spherical_Hashing(ref_buff,BCODE_LEN, seg_len);

		map->Suffix_Array(seg_len);

		T0.Stop();
		printf("- Index Finished (%f seconds)\n\n\n\n",T0.GetTime() );
	}else if (gate.compare("mapping") == 0)
	{		
    	Points read_1_buff;
    	Points read_2_buff;
    	int read_size = store_reads(read_1_buff,read_2_buff);

    	Mapping *map = new Mapping(dim,rpro);
		map->Load_Spherical_Hashing(read_size,BCODE_LEN,seg_len);
		map->Load_SA(seg_len);
		std::cout << "- filter:" << filter << std::endl;
		Stopwatch T2("");
		T2.Reset();     T2.Start();
		std::cout <<"- Analysis of read region ..." << std::endl;
		map->Get_Read_Region(read_1_buff,read_2_buff);
		T2.Stop();
		std::cout << "- Analysis Finished(" << T2.GetTime() << " seconds)" << std::endl;

		map->Hash_Mapping_with_SA(read_1_buff,read_2_buff);

	/*	map->Load_Ref_Info();
		map->Output_Result(PAIR_1,read_1_buff);
		map->Output_Result(PAIR_2,read_2_buff);*/

	 	SAMwriter *sp = new SAMwriter(DIM,read_size); 
	 	sp->Transfer_Info_From_Mapping(map->is_read_1_rev,map->is_read_2_rev);//,map->read_1_buff,map->read_2_buff);
	 	delete map;

	    sp->Analyse_Result_Pair(rpro);

	    T2.Reset();     T2.Start();
	    sp->Generate_SAM(read_1_buff,read_2_buff);
	    T2.Stop();
	    printf("- Generate SAM File Finished (%f seconds)\n",T2.GetTime());
		T0.Stop();
		printf("- Total Running Time (%f seconds)\n",T0.GetTime() );
	}
		
	return 0;
}