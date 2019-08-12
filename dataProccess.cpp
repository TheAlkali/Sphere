#include<iostream>
#include<string.h>
#include<time.h>
#include <algorithm>
#include <map>
//#include"kmerUtils.hpp"
#include "Points.h"
#include "Common.h"
#include"FastxParser.hpp"
//#include"SAMparser.hpp"
#include "Utils.hpp"

//#include"benchmark.hpp"

//void Analyse_Result_Spherical();
void Analyse_Result(std::string,bool flag,std::vector<std::vector<std::string>>&, std::map<std::string, int> &);
void Analyse_True_Result(std::vector<std::vector<std::string>>);
void Intersection_Of_Two_Results(std::vector<std::vector<std::string>> &,std::vector<std::vector<std::string>>  &, 
								std::map<std::string, int> &, std::map<std::string, int> &);
void Intersection_Of_Three_Results(std::vector<std::vector<std::string>> &,std::vector<std::vector<std::string>> &,std::vector<std::vector<std::string>> &,
								std::map<std::string, int> &, std::map<std::string, int> &, std::map<std::string, int> &);
void Intersection_Of_Four_Results(std::vector<std::vector<std::string>> &,std::vector<std::vector<std::string>> &,
								std::vector<std::vector<std::string>> &,std::vector<std::vector<std::string>> &,
								std::map<std::string, int> &, std::map<std::string, int> &, std::map<std::string, int> &, std::map<std::string, int> &);
void analyse_tid();
void difference_of_error(std::string,std::string);
void ref_of_art_sim(std::string);

int main(int argc, char const *argv[])
//int main()
{
//	kmerUtils kmutil;
	int klen = 50;
//	LSH lsh;
	int gate = 5;// std::atoi(argv[1]);

	char ref[] = "../reference/transcripts/Homo_sapiens.GRCh38.cdna.all.fa";
	char read1[] = "dataset/srrdata/SRR5337025_1.fastq";//"dataset/rapmap_reads_1.fastq";
	//"/home/yxt/Documents/work/RNA-seq/reads/single-cell-data/SRR5337/SRR5337026_2_1.fastq";
	//"/home/yxt/Documents/work/RNA-seq/reads/sim_reads_mismatch/sim_20M_1.fq";
	char read2[] = "dataset/srrdata/SRR5337025_2.fastq";//"dataset/rapmap_reads_2.fastq";
	//"/home/yxt/Documents/work/RNA-seq/reads/single-cell-data/SRR5337/SRR5337026_2.fastq";
	//"/home/yxt/Documents/work/RNA-seq/reads/sim_reads_mismatch/sim_20M_2.fq";
	singleSeqList *reflist = new singleSeqList();
	pairSeqList *preadlist = new pairSeqList();
//	refKmerList refklist;
			
	if (gate == 0){		
		//extract the mapping results of each aread to disk
	/*	SAMparser  sp;
		char test_sam[] = "res/res.sam";//"../rapmap_res/rapmap_res.sam";
		char test_res[] = "sphere_res.txt";
		sp.get_Sam(ref,test_sam,test_res,2922261,klen);*/
	}else if (gate == 1){		
	//	store_reads();
//		store_reads("../reads/single-cell-data/SRR5337/Simulate_GRCh38.1.fa",klen);
//		store_ref_kmers(ref,klen);
	}else if (gate == 2)
	{
		difference_of_error("analyse/unknown_errors_sphere.bin","analyse/unknown_errors_rapmap.bin");
	} else if (gate == 3){
		//merge all refs to one
	//	merge_ref_seq(ref,1);
	}else if (gate == 4){
	/*	fparser.read_fastx(read1,read2,preadlist);
		kmutil.set_kmer_length(klen);
		kmutil.genkmers_from_all_transcripts(reflist);
		delete reflist;
		refklist = kmutil.get_kmer_list();*/
	}else if (gate == 5)
	{
		std::vector<std::vector<std::string>> bowtie_ref_name_vec;
		std::vector<std::vector<std::string>> hisat_ref_name_vec;
		std::vector<std::vector<std::string>>rapmap_ref_name_vec;
		std::vector<std::vector<std::string>> sphere_ref_name_vec;
		std::map<std::string, int> read_name_vec;

	//	Analyse_Result_Spherical();

	//	Analyse_Result("../hisat2_res/hisat2_res.sam",false,hisat_ref_name_vec,2250009);
	//	Analyse_Result("../bowtie2_res/bowtie2_res.sam",false,bowtie_ref_name_vec,2522231);
	//	Analyse_Result("../rapmap_res/rapmap_res.sam",false,rapmap_ref_name_vec,2081302);

		ref_of_art_sim("dataset/srrdata/sim_read_small_150_1.fastq");
		std::cout << "analyse ref finished" << std::endl;
		
		Analyse_Result("res/sim_150_res.sam",false,bowtie_ref_name_vec, read_name_vec);
		/*bowtie_ref_name_vec.clear();
		Analyse_Result("../rapmap_res/sim_150_res.sam",false,bowtie_ref_name_vec, read_name_vec);
		bowtie_ref_name_vec.clear();

		Analyse_Result("../bowtie2_res/sim_150_res.sam",false,bowtie_ref_name_vec, read_name_vec);
		bowtie_ref_name_vec.clear();

		Analyse_Result("../hisat2_res/sim_150_res.sam",false,bowtie_ref_name_vec, read_name_vec);
		bowtie_ref_name_vec.clear();*/


		
		/*Analyse_Result("res/sim_250_res.sam",false,bowtie_ref_name_vec, read_name_vec);
		bowtie_ref_name_vec.clear();
		
		Analyse_Result("../rapmap_res/sim_250_res.sam",false,bowtie_ref_name_vec, read_name_vec);
		bowtie_ref_name_vec.clear();*/

		//Analyse_Result("../hisat2_res/sim_250_res.sam",false,bowtie_ref_name_vec, read_name_vec);
		//bowtie_ref_name_vec.clear();

		//Analyse_Result("res/250bp_57bits.sam",false,bowtie_ref_name_vec, read_name_vec);
		//bowtie_ref_name_vec.clear();
	
		//Analyse_Result("../rapmap_res/sim_250_res.sam",false,rapmap_ref_name_vec);
	}else if (gate == 6){
		std::vector<std::vector<std::string>> bowtie_ref_name_vec;
		std::vector<std::vector<std::string>> hisat_ref_name_vec;
		std::vector<std::vector<std::string>> rapmap_ref_name_vec;
		std::vector<std::vector<std::string>> sphere_ref_name_vec;

		std::map<std::string, int> bowtie_read_name_vec;
		std::map<std::string, int> hisat_read_name_vec;
		std::map<std::string, int> rapmap_read_name_vec;
		std::map<std::string, int> sphere_read_name_vec;
	//	Analyse_Result_Spherical();
		std::string sam = argv[2];
		std::cout << sam << std::endl;
		Analyse_Result("../hisat2_res/" + sam,true,hisat_ref_name_vec, hisat_read_name_vec);
		std::cout << "analyse hsiat2 finished" << std::endl;
		Analyse_Result("../bowtie2_res/" + sam,true,bowtie_ref_name_vec, bowtie_read_name_vec);
		std::cout << "analyse bowtie2 finished" << std::endl;
		Analyse_Result("../rapmap_res/" + sam,true,rapmap_ref_name_vec, rapmap_read_name_vec);
		std::cout << "analyse rapmap finished" << std::endl;
		Analyse_Result("res/" + sam,true,sphere_ref_name_vec, sphere_read_name_vec);
		std::cout << "analyse bit-mapping finished" << std::endl;

		
		std::cout << "hisat2 and rapmap" << std::endl;
		Intersection_Of_Two_Results(rapmap_ref_name_vec,hisat_ref_name_vec, rapmap_read_name_vec, hisat_read_name_vec);

		std::cout << "sphere and hisat2" << "\t" << sphere_ref_name_vec.size() << "\t" << hisat_ref_name_vec.size() << std::endl;
		Intersection_Of_Two_Results(sphere_ref_name_vec,hisat_ref_name_vec, sphere_read_name_vec, hisat_read_name_vec);

		std::cout << "hisat2 and bowtie2" << std::endl;
		Intersection_Of_Two_Results(bowtie_ref_name_vec,hisat_ref_name_vec, bowtie_read_name_vec, hisat_read_name_vec);

		std::cout << "rapmap and bowtie2" << std::endl;
		Intersection_Of_Two_Results(bowtie_ref_name_vec,rapmap_ref_name_vec, bowtie_read_name_vec, rapmap_read_name_vec);

		std::cout << "sphere and bowtie" << "\t" << sphere_ref_name_vec.size() << "\t" << bowtie_ref_name_vec.size() << std::endl;
		Intersection_Of_Two_Results(sphere_ref_name_vec,bowtie_ref_name_vec, sphere_read_name_vec, bowtie_read_name_vec);

		std::cout << "sphere and rapmap" << std::endl;
		Intersection_Of_Two_Results(sphere_ref_name_vec,rapmap_ref_name_vec, sphere_read_name_vec, rapmap_read_name_vec);

		
		std::cout << "sphere, rapmap and bowtie2" << std::endl;
		Intersection_Of_Three_Results(sphere_ref_name_vec,rapmap_ref_name_vec,bowtie_ref_name_vec, 
									sphere_read_name_vec, rapmap_read_name_vec, bowtie_read_name_vec);
		std::cout << "sphere, rapmap and hisat2" << std::endl;
		Intersection_Of_Three_Results(sphere_ref_name_vec,rapmap_ref_name_vec,hisat_ref_name_vec,
									sphere_read_name_vec, rapmap_read_name_vec, hisat_read_name_vec);
		std::cout << "sphere, hisat and bowtie2" << std::endl;
		Intersection_Of_Three_Results(sphere_ref_name_vec,hisat_ref_name_vec,bowtie_ref_name_vec,
									sphere_read_name_vec,hisat_read_name_vec,bowtie_read_name_vec);
		std::cout << "hisat, rapmap and bowtie2" << std::endl;
		Intersection_Of_Three_Results(hisat_ref_name_vec,rapmap_ref_name_vec,bowtie_ref_name_vec,
									hisat_read_name_vec,rapmap_read_name_vec,bowtie_read_name_vec);

		std::cout << "hisat, rapmap, hisat2 and bowtie2" << std::endl;
		Intersection_Of_Four_Results(hisat_ref_name_vec,rapmap_ref_name_vec,bowtie_ref_name_vec,sphere_ref_name_vec,
									hisat_read_name_vec,rapmap_read_name_vec,bowtie_read_name_vec,sphere_read_name_vec);
		
	}else if (gate == 9)
	{
		analyse_tid();
	}
	return 0;
} 

void analyse_tid(){
	std::vector<std::string> ref_name;
	{
		std::ifstream tid_file("analyse/rname_ref.bin");
		cereal::BinaryInputArchive ar(tid_file);
		ar(ref_name);
	}

	std::vector<std::string> flux_ref_name;
	{
		std::ifstream tid_file("analyse/flux_tid.bin");
		cereal::BinaryInputArchive ar(tid_file);	
		ar(flux_ref_name);
	}
	int count = 0;
	std::vector<std::string> same_tid;
	std::string ref_tmp;
	for (int i = 0; i < flux_ref_name.size(); ++i){
	//	std::cout << flux_ref_name[i] << std::endl;
		for (int j = 0; j < ref_name.size(); ++j)
		{
			ref_tmp = ref_name[j].substr(0,ref_name[j].find_first_of("."));
			if ((flux_ref_name[i]).compare(ref_tmp) == 0)
			{
				count++;
				same_tid.push_back(ref_tmp);
				break;
			}
		}
	}
	std::cout << "ref tid count :" << ref_name.size() << std::endl;
	std::cout << "flux tid count :" << flux_ref_name.size() << std::endl;
	std::cout << "same tid count :" << same_tid.size() << std::endl;
	{
		std::ofstream file("../flux_test/sam_tid.bin");
		cereal::BinaryOutputArchive ar(file);
		ar(same_tid);
	}
}

void difference_of_error(std::string file1,std::string file2)
{
	std::vector<std::string> error1,error2;
	{
		std::ifstream file(file1);
		cereal::BinaryInputArchive ar(file);
		ar(error1);
	}
	{
		std::ifstream file(file2);
		cereal::BinaryInputArchive ar(file);
		ar(error2);
	}
	std::ofstream out("analyse/sphere_error.txt");
	std::ofstream out2("analyse/same_error.txt");
	int j = 0;
	for (int i = 0; i < error1.size(); ++i)
	{
		for (j = 0; j < error2.size(); ++j)
		{
			if (error1[i].compare(error2[j]) == 0)
			{
				out2 << error1[i] << std::endl;
				break;
			}
		}
		if (j == error2.size())
		{
			out << error1[i] << std::endl;
		}
	}
	out.close();
	out2.close();
}

void ref_of_art_sim(std::string filename)
{
	std::ifstream file(filename);
	std::string line, ref;
	std::vector<std::string> ref_of_read;
	for (int i = 1; file.peek() != EOF; ++i)
	{
		getline(file, line);
		if (i % 4 == 1)
		{
			ref = line.substr(1, line.find_first_of("-") - 1);
			ref_of_read.push_back(ref);
		}
	}

	{
        std::ofstream out("../art/ref_of_mapping_read.bin");
        cereal::BinaryOutputArchive ar(out);
        ar(ref_of_read);
    }
}
/*
std::vector<std::string> Save_Ref_Of_Read()
{
	std::vector<std::vector<int> > ref_of_read_1;
	std::vector<std::vector<int> > ref_of_read_2;
	std::vector<std::string> ref_name;
	std::vector<int> loc_to_ref;
	std::vector<int> min_dis_1;
	std::vector<int> min_dis_2;
	{
        std::ifstream ref_name_file("bin/rname_ref.bin");
        cereal::BinaryInputArchive ar_ref_name(ref_name_file);
        ar_ref_name(ref_name);
    }

    {
        std::ifstream ref_of_read_file("dataset/ref_of_read_1.bin");
        cereal::BinaryInputArchive ar(ref_of_read_file);
        ar(ref_of_read_1);
    }
    {
        std::ifstream ref_of_read_file("dataset/ref_of_read_2.bin");
        cereal::BinaryInputArchive ar(ref_of_read_file);
        ar(ref_of_read_2);
    }

    {
    	std::ifstream file("tmp/tmp_dis_1.bin");
        cereal::BinaryInputArchive ar(file);
        ar(min_dis_1);
    }

    {
    	std::ifstream file("tmp/tmp_dis_2.bin");
        cereal::BinaryInputArchive ar(file);
        ar(min_dis_2);
    }

    std::vector<std::string> read_name;
    {
        std::ifstream read_name_file(PAIR_1_NAME_FILE);
        cereal::BinaryInputArchive ar(read_name_file);
        ar(read_name);
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<int> intersection;
    std::vector<int>::iterator it;
    std::vector<std::vector<std::string>> ref_of_reads;
    std::vector<std::string> ref;
    std::string tmp;
    int read_size = ref_of_read_1.size();
    int error = 0;
    int half_right = 0;
	std::ofstream error_read_half("analyse/error_read_half.txt");
	std::ofstream error_read_all("analyse/error_read_all.txt");
	std::vector<std::string> error_read_half_vec;
    for (int i = 0; i < read_size; ++i)
    {
        intersection.clear();
        sort(ref_of_read_1[i].begin(),ref_of_read_1[i].end());
        sort(ref_of_read_2[i].begin(),ref_of_read_2[i].end());
        set_intersection(ref_of_read_1[i].begin(),ref_of_read_1[i].end(),ref_of_read_2[i].begin(),ref_of_read_2[i].end(),inserter(intersection,intersection.begin()));
        if (intersection.size() == 0)
        {
        	intersection.clear();
        	if (min_dis_1[i] <= filter && min_dis_2[i] > min_dis_1[i])
        	{
        		intersection = ref_of_read_1[i];
        		half_right++;
        	//	error_read_half << read_name[i] << "\t 1" << "\n";
        		error_read_half_vec.push_back(read_name[i]);
        	}else if (min_dis_2[i] <= filter && min_dis_1[i] > min_dis_2[i])
        	{
        		intersection = ref_of_read_2[i];
        		half_right++;
        	//	error_read_half << read_name[i] << "\t 2" << "\n";
        		error_read_half_vec.push_back(read_name[i]);
        	}else if (min_dis_1[i] > filter && min_dis_2[i] > filter)
        	{
        		error++;
        		error_read_all << read_name[i] << "\n";
        	}
        }
        ref.clear();
        for (it = intersection.begin(); it != intersection.end(); ++it)
        {
        	tmp = ref_name[*it];
        	tmp = tmp.substr(0,tmp.find_first_of("."));
        //	std::cout << tmp << std::endl;
            ref.push_back(tmp);
        }
    //    std::cout << std::endl;

        if (ref.size() != 0)
        {
        	ref_of_reads.push_back(ref);
        }else
        {
        	ref.push_back("*");
        	ref_of_reads.push_back(ref);
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double fastxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    std::cout << "get the final result:" << fastxparser_time << "s" << std::endl;
    std::cout << "both distances are larger than 2:" << error << std::endl;
    std::cout << "half right:" << half_right << std::endl;
    {
        std::ofstream file("analyse/ref_of_mapping_read.bin");
        cereal::BinaryOutputArchive ar(file);
        ar(ref_of_reads);
    }  
    {
    	std::ofstream file("analyse/half_mapping_read.bin");
        cereal::BinaryOutputArchive ar(file);
        ar(error_read_half_vec);
    }
    error_read_all.close();
    error_read_half.close();
    return error_read_half_vec;
}*/

void Analyse_Sim_Result(std::vector<std::vector<std::string>> ref_of_reads,std::vector<std::string> true_ref_of_reads,bool flag,std::vector<std::string> error_read_all = std::vector<std::string>(1,"test"))
{
/*	std::vector<std::string> read_name;
    {
        std::ifstream read_name_file(PAIR_1_NAME_FILE);
        cereal::BinaryInputArchive ar(read_name_file);
        ar(read_name);
    }*/


    std::ofstream error_read_file;

    error_read_file.open("analyse/error_read_sphere.txt");
	std::vector<std::string> error_read;
	std::vector<std::string> ref;
	int total = true_ref_of_reads.size();
    int mapped = 0,error = 0,unmapped  =0;
    int j = 0;
    for (int i = 0; i < ref_of_reads.size(); ++i)
    {
        ref = ref_of_reads[i];
        //std::cout <<"------------"<< std::endl;
        //std::cout <<true_ref_of_reads[i] << std::endl;
        for (j = 0; j < ref.size(); ++j)
        {
        	//std::cout << ref[j] << std::endl;
	        if (ref[j] == "*")
	        {
	        	unmapped++;
	        	j = 0;
	        	break;
	        }
            if (ref[j].compare(true_ref_of_reads[i]) == 0)
            {
                mapped++;
                j = 0;
                break;
            }
        }
        if (j == ref.size())
        {
        	error++;
        //	error_read_file << read_name[i] << "\n";
        //	error_read.push_back(read_name[i]);
        }
    }
    std::cout << "correct mapped reads:" << (float)mapped << std::endl;
    std::cout << "incorrect mapped reads:" << (float)error << std::endl;
    std::cout << "unmapped reads:" << (float)unmapped << std::endl;
    std::cout << "total reads:" << total << std::endl;
    std::cout << "correct mapped ratio:" << (float)mapped / total << std::endl;
    std::cout << "incorrect mapped ratio:" << (float) error/ total << std::endl;
    std::cout << "unmapped ratio:" << (float)unmapped / total << std::endl << std::endl;


    if (error_read.size() > 0)
    {
    	std::ofstream unknown_error("analyse/unknown_errors_sphere.bin");
    	cereal::BinaryOutputArchive ar(unknown_error);
    	ar(error_read);
    }
}
/*
void Analyse_Result_Spherical()
{
	std::vector<std::string> error_read_half_vec = Save_Ref_Of_Read();
    std::vector<std::vector<std::string>> ref_of_reads;
    {
        std::ifstream file("analyse/ref_of_mapping_read.bin");
        cereal::BinaryInputArchive ar(file);
        ar(ref_of_reads);
    }
    std::vector<std::string> true_ref_of_reads;
    {
        std::ifstream file("analyse/true_ref_of_sim_read.bin");
        cereal::BinaryInputArchive ar(file);
        ar(true_ref_of_reads);
    } 
    Analyse_Sim_Result(ref_of_reads,true_ref_of_reads,true_ref_of_reads.size(),true,error_read_half_vec) ;
}*/

void Analyse_Result(std::string filename,bool flag,std::vector<std::vector<std::string>> &ref_of_reads, std::map<std::string, int> &read_name_vec)
{
	std::ifstream sam(filename);
//	std::ifstream sam("res/res.sam");

	std::vector<std::string> true_ref_of_reads;
/*	{
        std::ifstream rfile("../art/ref_of_mapping_read.bin");
        //std::ifstream rfile("analyse/true_ref_of_sim_read.bin");
        cereal::BinaryInputArchive ar(rfile);
        ar(true_ref_of_reads);
    }*/
    std::cout << "load ref of reads finished" << std::endl;
	char tmp;
	std::string line;
	sam.get(tmp);
	while(tmp == '@')
	{
		getline(sam,line);
		sam.get(tmp);
	}
	sam.seekg(-1,std::ios_base::cur);

	std::string ref_name;
	std::string read_name;
	std::string last_read_name;
	std::vector<std::string> ref_name_vec;
	

	std::string true_ref_name;

	std::string str;
	int count = 0;
	getline(sam,line);
	count++;
	std::stringstream ss(line);
	ss >> read_name >> str >> ref_name;
	//ref_name = ref_name.substr(0,ref_name.find_first_of("."));
//	ss >> read_name;
//	ref_name = read_name.substr(0,ref_name.find_first_of("-"));
	last_read_name = read_name;
	//read_name_vec[read_name] = 0;
	int read_size = 0;
	for(int i = 0;sam.peek() != EOF;i++)
	{
		if(read_name.compare(last_read_name) == 0)
		{
			ref_name_vec.push_back(ref_name);
			//std::cout << ref_name << std::endl;
		}else
		{
			//std::cout << read_name << std::endl;
		//	std::cout << std::endl;
			read_size++;
			ref_of_reads.push_back(ref_name_vec);
			ref_name_vec.clear();
			//read_name_vec[read_name] = read_size;
			ref_name_vec.push_back(ref_name);

			//last_read_name = last_read_name.substr(last_read_name.find_first_of(":") + 1,last_read_name.size());
			//last_read_name = last_read_name.substr(last_read_name.find_first_of(":") + 1,last_read_name.size());
			
			if (!flag)
			{
				//true_ref_name = last_read_name.substr(0,last_read_name.find_first_of(":"));
				true_ref_name = last_read_name.substr(0, last_read_name.find_first_of("-"));
				//std::cout << true_ref_name << std::endl;
				true_ref_of_reads.push_back(true_ref_name);
			}
		}
			
		last_read_name = read_name;
		line.clear();
		getline(sam,line);
		count++;
		std::stringstream ss(line);
		ss >> read_name >> str >> ref_name;
		//ref_name = ref_name.substr(0,ref_name.find_first_of("."));

	//	ss >> read_name;
	//	ref_name = read_name.substr(0,ref_name.find_first_of("-"));
	}
	std::cout << filename << std::endl;
	if (!flag)
	{
		//std::cout << "read size:" << ref_of_reads.size() << std::endl;
		Analyse_Sim_Result(ref_of_reads,true_ref_of_reads,false);
	}else
	{
		Analyse_True_Result(ref_of_reads);
	}
}

void Analyse_True_Result(std::vector<std::vector<std::string>> ref_of_reads)
{
    int mapped = 0,error = 0,unmapped  =0;
    int read_size = ref_of_reads.size();
#ifdef USE_PARALLELIZATION
    #pragma omp parallel for reduction(+:mapped,unmapped) num_threads(THREAD)
#endif
    for (int i = 0; i < read_size; ++i)
    {
        if (ref_of_reads[i][0] == "*")
        {
        	unmapped++;
        }else
        {
        	mapped++;
        }
    }
    std::cout << "total reads:" << read_size << std::endl;
    std::cout << "unmapped reads:" << unmapped  << "\t" << (float)unmapped / read_size << std::endl;
    std::cout << "mapped reads:" << mapped << "\t" << (float)mapped / read_size << std::endl << std::endl;
}

void Intersection_Of_Two_Results(std::vector<std::vector<std::string>> &ref_of_reads_1,std::vector<std::vector<std::string>> &ref_of_reads_2,
								std::map<std::string, int> &read_name_vec_1, std::map<std::string, int> &read_name_vec_2)
{
	std::vector<std::string> ref_1,ref_2;
	int count = 0;
	bool flag = true;
	int read_size = ref_of_reads_1.size();
	int j = 0, i = 0;

	std::map<std::string, int>::iterator it_1 = read_name_vec_1.begin();
	//std::map<std::string, int>::iterator it_2 = read_name_vec_2.begin();
	std::string read_name;

	while(it_1 != read_name_vec_1.end())
	{
		read_name = it_1->first;
		i = it_1->second;
		j = read_name_vec_2[read_name];
		ref_1 = ref_of_reads_1[i];
		ref_2 = ref_of_reads_2[j];
		if (ref_1.size() == ref_2.size() && ref_1[0] != "*" && ref_2[0] != "*")
		{
			/*std::cout << "------" << std::endl;
			std::cout << read_name_vec_1[read_name] << std::endl;
			for (int i = 0; i < ref_1.size(); ++i)
			{
				std::cout << ref_1[i] << std::endl;
			}
			std::cout << std::endl;
			std::cout << read_name_vec_2[read_name] << std::endl;
			for (int i = 0; i < ref_2.size(); ++i)
			{
				std::cout << ref_2[i] << std::endl;
			}*/
			sort(ref_1.begin(),ref_1.end());
			sort(ref_2.begin(),ref_2.end());
			for (int j = 0; j < ref_1.size(); ++j)
			{
				if (ref_1[j].compare(ref_2[j]) != 0)
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				count++;
			}
			flag = true;
		}
		*it_1++;
	}
	std::cout << "the same mapping ratio of them:" << count << "\t" <<(float)count / read_size << std::endl;
}

void Intersection_Of_Three_Results(std::vector<std::vector<std::string>> &ref_of_reads_1,std::vector<std::vector<std::string>> &ref_of_reads_2
								,std::vector<std::vector<std::string>> &ref_of_reads_3, std::map<std::string, int> & read_name_vec_1, 
								std::map<std::string, int> & read_name_vec_2, std::map<std::string, int> & read_name_vec_3)
{
	std::vector<std::string> ref_1,ref_2,ref_3;
	int count = 0;
	bool flag = true;
	int read_size = ref_of_reads_1.size();
	int j = 0, i = 0, k = 0;;
	std::map<std::string, int>::iterator it_1 = read_name_vec_1.begin();
	//std::map<std::string, int>::iterator it_2 = read_name_vec_2.begin();
	std::string read_name;

	while(it_1 != read_name_vec_1.end())
	{
		read_name = it_1->first;
		//std::cout << read_name << std::endl;

		i = it_1->second;
		j = read_name_vec_2[read_name];
		k = read_name_vec_3[read_name];

		ref_1 = ref_of_reads_1[i];
		ref_2 = ref_of_reads_2[j];
		ref_3 = ref_of_reads_3[k];
		if (ref_1.size() == ref_2.size() && ref_1.size() == ref_3.size() && ref_1[0] != "*" && ref_2[0] != "*" && ref_3[0] != "*")
		{
			sort(ref_1.begin(),ref_1.end());
			sort(ref_2.begin(),ref_2.end());
			sort(ref_3.begin(),ref_3.end());
			for (int i = 0; i < ref_1.size(); ++i)
			{
				if (ref_1[i].compare(ref_2[i]) == 0 && ref_1[i].compare(ref_3[i]) == 0)
				{
					continue;
				}else
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				count++;
			}
			flag = true;
		}
		*it_1++;
	}
	std::cout << "the same mapping ratio of them:" << count << "\t" <<(float)count / read_size << std::endl;
}

void Intersection_Of_Four_Results(std::vector<std::vector<std::string>> &ref_of_reads_1,std::vector<std::vector<std::string>> &ref_of_reads_2
								,std::vector<std::vector<std::string>> &ref_of_reads_3,std::vector<std::vector<std::string>> &ref_of_reads_4,
								std::map<std::string, int> & read_name_vec_1, std::map<std::string, int> & read_name_vec_2,
								std::map<std::string, int> & read_name_vec_3, std::map<std::string, int> & read_name_vec_4)
{
	std::vector<std::string> ref_1,ref_2,ref_3,ref_4;
	int count = 0;
	bool flag = true;
	int read_size = ref_of_reads_1.size();
	int j = 0, i = 0, k = 0, s = 0;
	std::map<std::string, int>::iterator it_1 = read_name_vec_1.begin();
	//std::map<std::string, int>::iterator it_2 = read_name_vec_2.begin();
	std::string read_name;

	while(it_1 != read_name_vec_1.end())
	{
		read_name = it_1->first;
		i = it_1->second;
		j = read_name_vec_2[read_name];
		k = read_name_vec_3[read_name];
		s = read_name_vec_4[read_name];

		ref_1 = ref_of_reads_1[i];
		ref_2 = ref_of_reads_2[j];
		ref_3 = ref_of_reads_3[k];
		ref_4 = ref_of_reads_4[s];
		if (ref_1.size() == ref_2.size() && ref_1.size() == ref_3.size()  && ref_1.size() == ref_4.size() && ref_1[0] != "*" && ref_2[0] != "*" && ref_3[0] != "*"
			 && ref_4[0] != "*")
		{
			sort(ref_1.begin(),ref_1.end());
			sort(ref_2.begin(),ref_2.end());
			sort(ref_3.begin(),ref_3.end());
			sort(ref_4.begin(),ref_4.end());
			for (int i = 0; i < ref_1.size(); ++i)
			{
				if (ref_1[i].compare(ref_2[i]) == 0 && ref_1[i].compare(ref_3[i]) == 0 && ref_1[i].compare(ref_4[i]) == 0)
				{
					continue;
				}else
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				count++;
			}
			flag = true;
		}
		*it_1++;
	}
	std::cout << "the same mapping ratio of them:" << count << "\t" <<(float)count / read_size << std::endl;
}