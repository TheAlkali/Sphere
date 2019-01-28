#include "BinaryHash.h"
#include "Evaluation.h"
#include  <iostream>


int iteration = 0;
bool is_spherical_hash_learned = false;
Points ref_buff, read_buff_1,read_buff_2;

Points train_bcode_P;

SphericalHashing src_sh;
SphericalHashing bcode_sh;

void Initialize_Bcode_Data()
{
	train_bcode_P.Initialize(NUM_TRAIN_SAMPLES,BCODE_LEN);

	bitset<BCODE_LEN> *train_bcode_B = new bitset<BCODE_LEN> [NUM_TRAIN_SAMPLES];
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<NUM_TRAIN_SAMPLES;i++)
	{
		src_sh.Compute_BCode( ref_buff.d[i] , train_bcode_B[i] );
	}
	train_bcode_P.Initialize_From_Bcodes(train_bcode_B);
}

// initialize data and query points
int Initialize_Read_Data(size_t buff_size){
#ifdef INPUT_READ_FILE_NAME_1
	read_buff_1.Initialize(buff_size,50);
	int file_end = read_buff_1.Initialize_From_File();
#endif
	nQ = read_buff_1.nP;
	// you can control the number of queries here
	return file_end;
}

int Initialize_Ref_Data(size_t buff_size)
{
#ifdef INPUT_REF_FILE_NAME
	ref_buff.Initialize(buff_size,50);
	int file_end = ref_buff.Initialize_From_File();
#endif
	nP = ref_buff.nP;			
	//if file end, return the count of unvisited data
	//if file not end, return -1
	return file_end;
}

void Compute_Reads_BCodes(int is_readfile_end,int read_part)
{
	Stopwatch T0("");
	int qps_count = nQ;
	
	if (is_readfile_end != -1)
	{
		qps_count = is_readfile_end;
	}

	bitset<BCODE_LEN> *bCodeRead_SH = new bitset<BCODE_LEN> [qps_count];
	T0.Reset();		T0.Start();
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<qps_count;i++)
	{
		src_sh.Compute_BCode<REAL_TYPE>(read_buff_1.d[i], bCodeRead_SH[i] );
	}
	T0.Stop();
	printf("- Spherical Hashing: Computing Reads Binary Codes (%f seconds)\n",T0.GetTime() );

	unsigned long code;

	T0.Reset();		T0.Start();
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<qps_count;i++)
	{
		read_buff_1.codefile << bCodeRead_SH[i] << std::endl;
		code = bCodeRead_SH[i].to_ulong() ;
		read_buff_1.codefile.write( reinterpret_cast<const char*>(&code), sizeof(code) ) ;
	}
	T0.Stop();
	printf("- Save Reads Hashcode To Disk (%f seconds)\n",T0.GetTime() );

	delete [] bCodeRead_SH;
}

void Compute_Ref_BCodes(int is_reffile_end,int ref_part)
{
	Stopwatch T0("");
	int dps_count = nP;	
	if (is_reffile_end != -1)
	{
		dps_count = is_reffile_end;
	}
	bitset<BCODE_LEN> *bCodeRef_SH = new bitset<BCODE_LEN> [dps_count];
	bitset<BBCODE_LEN> *bCodeIdx = new bitset<BBCODE_LEN> [dps_count];

	T0.Reset();		T0.Start();
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<dps_count;i++)
	{
		src_sh.Compute_BCode<REAL_TYPE>( ref_buff.d[i] , bCodeRef_SH[i] );
		bcode_sh.Compute_BCode_for_BCode(bCodeRef_SH[i],bCodeIdx[i]);
	}
	T0.Stop();
	iteration++;
	printf("- Spherical Hashing %d : Computing Reference Binary Codes (%f seconds)\n",iteration,T0.GetTime() );

	T0.Reset();		T0.Start();
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<dps_count;i++)
	{
		std::string filename = "index/ref/" + bCodeIdx[i].to_string() + ".bin";
		std::ofstream tmp;
		// the src location
		size_t ind = ref_part + i;
		tmp.open(filename,std::ios::out|std::ios::app);
		tmp.write((char*)&ind,sizeof(ind));
		tmp.close();
	}
	T0.Stop();
	printf("- Save Index Of Reference Hashcode To Disk (%f seconds)\n",T0.GetTime() );

	unsigned long code;
	T0.Reset();		T0.Start();
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<dps_count;i++)
	{
		ref_buff.codefile << bCodeRef_SH[i] << std::endl;
		code = bCodeRef_SH[i].to_ulong();
		ref_buff.codefile.write( reinterpret_cast<const char*>(&code), sizeof(code) ) ;
	}
	T0.Stop();
	printf("- Save Reference Hashcode To Disk (%f seconds)\n",T0.GetTime() );

	delete [] bCodeRef_SH;
	delete [] bCodeIdx;
}

void Read_Hashcode_From_File(std::string filename, bitset<BCODE_LEN> *bCode)
{
	std::ifstream bcode_input;
	std::ifstream hind_input; 
	unsigned long code ;
	bcode_input.open(filename);
	
	if (!bcode_input.is_open())
	{
		std::cerr << "- can't open file " << filename << std::endl;
		exit(EXIT_FAILURE);
	}

	for(int i = 0;bcode_input.peek() != EOF;i++)
	{
		std::string tmp;
		getline(bcode_input,tmp);
		
		bcode_input.read( reinterpret_cast<char*>(&code), sizeof(code) ) ;
		bCode[i] = String2Bit<BCODE_LEN>(tmp);
	}

	bcode_input.close();
}

void Hash_Mapping()
{
	system("rm index/ref/*.bin");
	// compute hashcode
	int ref_part = 0;
	int read_part = 0;
	int is_readfile_end = -1;
	int is_reffile_end = -1;
	while(true)
	{	
		if (is_reffile_end == -1 )
		{
			is_reffile_end = Initialize_Ref_Data(REF_BUFFER_SIZE);
			Compute_Ref_BCodes(is_reffile_end,ref_part);
			ref_part += REF_BUFFER_SIZE;
		}
		if (is_readfile_end == -1)
		{
			is_readfile_end = Initialize_Read_Data(READ_BUFFER_SIZE);
			Compute_Reads_BCodes(is_readfile_end,read_part);
			read_part += READ_BUFFER_SIZE;
		}
		if (is_reffile_end != -1 && is_readfile_end != -1)
		{
		//	test_speed();
			break;
		}
	}

	ref_buff.srcfile.close();
	read_buff_1.srcfile.close();
	ref_buff.codefile.close();
	read_buff_1.codefile.close();
	ref_buff.ReleaseMem();
	read_buff_1.ReleaseMem();

	Stopwatch T0("");

	bitset<BCODE_LEN> *bCodeRead_SH = new bitset<BCODE_LEN> [nQ];
	bitset<BCODE_LEN> *bCodeRef_SH = new bitset<BCODE_LEN> [nP];

	bitset<BBCODE_LEN> bCodeRead_ind;

	T0.Reset();		T0.Start();
	Read_Hashcode_From_File(OUTPUT_REF_HASH_FILE_NAME,bCodeRef_SH);
	T0.Stop();
	printf("- Load Reference Hashcode Finished (%f seconds)\n",T0.GetTime() );

	T0.Reset();		T0.Start();
	Read_Hashcode_From_File(OUTPUT_READ_HASH_FILE_NAME,bCodeRead_SH);
	T0.Stop();
	printf("- Load Reads Hashcode Finished (%f seconds)\n",T0.GetTime() );

	T0.Reset();		T0.Start();
	printf("- Start Mapping\n");
	
	std::vector<size_t> ref_loc;
	std::ifstream ref;
	double min_dis;
	std::string filename;
	//std::ofstream tmp_codeidx;
	//tmp_codeidx.open("code_index.txt");


	for(int qIndex = 0;qIndex < nQ;qIndex++)
	{
		min_dis = 100000;

//-------------------------------------------------
//		--get search region usign hash index--		

		bcode_sh.Compute_BCode_for_BCode(bCodeRead_SH[qIndex],bCodeRead_ind);

		//tmp_codeidx << bCodeRead_ind << std::endl;

		filename = "index/ref/" + bCodeRead_ind.to_string() + ".bin";
		
//		--get mapping region using hash index--
		ref.open(filename);
		for (int i = 0; ref.peek() != EOF; ++i)
		{
			size_t ind;
			ref.read((char*)(&ind),sizeof(ind));
			int dist = Compute_HD(bCodeRead_SH[qIndex],bCodeRef_SH[ind]);
			if (min_dis > dist)
			{
				min_dis = std::move(dist);
				ref_loc.clear();
				ref_loc.push_back(std::move(ind));
			}else if (min_dis == dist)
			{
				ref_loc.push_back(std::move(ind));
			}
		}
		ref.close();
//--------------------------------------------------
//      --search in global--
/*		for (int i = 0; i < nP; ++i)
		{
			int dist = Compute_HD(bCodeRead_SH[qIndex],bCodeRef_SH[i]);
			if (min_dis > dist)
			{
				min_dis = std::move(dist);
				ref_loc.clear();
				ref_loc.push_back(i);
			}else if (min_dis == dist)
			{
				ref_loc.push_back(i);
			}
		}*/

		mres.mapped_ref_loc.push_back(ref_loc);
		mres.min_dis.push_back(min_dis);
		
	}
	T0.Stop();
	//tmp_codeidx.close();
/*	delete [] bCodeRead_SH;
	delete [] bCodeRef_SH;*/
	printf("- Mapping time: Mapping Reads Through Hash Codes Finished (%f seconds)\n",T0.GetTime() );
}