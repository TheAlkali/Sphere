#include<iostream>
#include<string.h>
#include<time.h>
#include <algorithm>
//#include"kmerUtils.hpp"
#include"fastxParser.hpp"
#include"SAMparser.hpp"
#include "Utils.hpp"

//#include"benchmark.hpp"

void gen_simulate_reads(char*,char*, singleSeqList*,int);
void gen_simulate_reads_From_TXT(int,int);
void Analyse_Result_Spherical();
void count_ref_kmer(std::string);
void Analyse_Result_Others();
void analyse_tid();
void rand_prob();

int main(int argc, char const *argv[])
//int main()
{
//	kmerUtils kmutil;
	int klen = 50;
//	LSH lsh;
	int gate = std::atoi(argv[1]);

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
		SAMparser  sp;
		char test_sam[] = "../hisat2_res/hisat2_res.sam";
		char test_res[] = "hisat2_res.txt";
		sp.get_Sam(ref,test_sam,test_res,20748,klen);
	}else if (gate == 1){		
		store_reads();
//		store_reads("../reads/single-cell-data/SRR5337/Simulate_GRCh38.1.fa",klen);
//		store_ref_kmers(ref,klen);
	}else if (gate == 2){		
		gen_simulate_reads(read1,ref,reflist,klen);
		gen_simulate_reads(read2,ref,reflist,klen);
	}else if (gate == 3){
		//merge all refs to one
		merge_ref_seq(ref,1);
	}else if (gate == 4){
	/*	fparser.read_fastx(read1,read2,preadlist);
		kmutil.set_kmer_length(klen);
		kmutil.genkmers_from_all_transcripts(reflist);
		delete reflist;
		refklist = kmutil.get_kmer_list();*/
	}else if (gate == 5){
		gen_simulate_reads_From_TXT(238956,klen);
	}else if (gate == 6){
	//	Analyse_Result_Spherical();
		Analyse_Result_Others();
	}else if (gate == 7){
		//ref kmer count:282626422
		count_ref_kmer("dataset/ref_kmer.txt");
	}else if (gate == 8){
	//	store_reads_name(read1);
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

void count_ref_kmer(std::string filename){
	std::ifstream ref_kmer;
	ref_kmer.open(filename);
	size_t count = 0;
	std::string tmp;
	while(ref_kmer.peek() != EOF){
		getline(ref_kmer,tmp);
		count++;
	}
	std::cout << "ref kmer count:" << count << std::endl;
}

//generate reads from transcriptpmes randomly
void gen_simulate_reads(char* read_name,char* ref_name,singleSeqList *reflist,int read_len){
	std::ofstream sim_read_file_1;
	sim_read_file_1.open(read_name);

	read_fastx(ref_name,reflist);
	int ref_size = reflist->seq.size();
	int gen_read_per_trans = 100;

	std::string tmp_ref;
	std::string tmp_ref_name;
	for (int i = 0; i < ref_size; ++i){
		tmp_ref.clear();
		tmp_ref_name.clear();
		tmp_ref = reflist->seq[i];
		tmp_ref_name = reflist->name[i];
		srand(time(NULL));
		for (int j = 0; j < gen_read_per_trans; ++j){
			int rand_loc = rand() % (tmp_ref.length() - read_len + 1);
			std::string aread = tmp_ref.substr(rand_loc,read_len);

			sim_read_file_1 << '@' << j + 1 << ":SIMUALATE:" << tmp_ref_name << std::endl;
			sim_read_file_1 << aread << std::endl;
			sim_read_file_1 << '+' << std::endl;
			sim_read_file_1 << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" << std::endl;

			//mismatch
			int rand_mis;
			int mis_num = 5;
			int rand_base = 0;
			const char* rand_base_str = "0";
			for (int k = 0; k < mis_num; ++k){
				rand_mis = rand() % (read_len);
				while(rand_mis < 5){
					rand_mis = rand() % (read_len);
				}
				rand_base = aread[rand_mis] - 'A' + 1;
				while((char)rand_base == aread[rand_mis] - 'A' + 1){
					rand_base = rand() % (3) + 1;

					if (rand_base == 1)
					{
						rand_base_str = "A";
					}else if (rand_base == 2)
					{
						rand_base_str = "T";
					}else if (rand_base == 3)
					{
						rand_base_str = "C";
					}else if (rand_base == 4)
					{
						rand_base_str = "G";
					}
				}
				
				aread.replace(rand_mis,1,rand_base_str);

				//std::cout << aread[rand_mis] << std::endl;
			}

			sim_read_file_1 << '@' << j + 1 << ":SIMUALATE-M:" << tmp_ref_name << std::endl;
			sim_read_file_1 << aread << std::endl;
			sim_read_file_1 << '+' << std::endl;
			sim_read_file_1 << "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"<<std::endl;
		}
	}
	std::cout << "Simulate Reads From Reference Finished." << std::endl;
	sim_read_file_1.close();
}

void gen_simulate_reads_From_TXT(int size,int dim){
	std::string aread;
	std::ofstream output;
	std::ifstream input;
	input.open("dataset/ref238956.txt");
	output.open("dataset/sim_mis_read_3.txt");
	for (int i = 0; i < 5000; ++i){
		int rand_loc = rand() % (size - dim + 1);
		input.seekg(rand_loc * (dim + 1),std::ios_base::beg);
		getline(input,aread);
		
		output << aread << std::endl;
		//mismatch
		int rand_mis;
		int mis_num = 5;
		int rand_base;
		for (int i = 0; i < mis_num; ++i){
			rand_mis = rand() % (dim);
			while(rand_mis < 5)
			{
				rand_mis = rand() % (dim);
			}
			rand_base = rand() % (3) + 1;
			while((char)rand_base == aread[rand_mis]){
				rand_base = rand() % (3) + 1;
			}
			aread.replace(rand_mis,1, std::to_string(rand_base));
		}

		//insert
	/*	aread.insert(5,"3");
		aread.pop_back();*/
		output << aread << std::endl;
	}
	input.close();
	output.close();
}

bool hamming_dis(std::string s1,std::string s2,int dim){
	int count = 0;
	for (int i = 0; i < dim; ++i){
		if (s1[i] != s2[i]){
			count++;
		}
	}
	if (count >1){
		return false;
	}else{
		return true;
	}
}


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
	std::vector<std::string> error_read_all_vec;
    for (int i = 0; i < read_size; ++i)
    {
        intersection.clear();
        sort(ref_of_read_1[i].begin(),ref_of_read_1[i].end());
        sort(ref_of_read_2[i].begin(),ref_of_read_2[i].end());
        set_intersection(ref_of_read_1[i].begin(),ref_of_read_1[i].end(),ref_of_read_2[i].begin(),ref_of_read_2[i].end(),inserter(intersection,intersection.begin()));
        if (intersection.size() == 0)
        {
        	intersection.clear();
        	if (min_dis_1[i] <= 1 && min_dis_2[i] > min_dis_1[i])
        	{
        		intersection = ref_of_read_1[i];
        		half_right++;
        		error_read_half << read_name[i] << "\t 1" << "\n";
        	}else if (min_dis_2[i] <= 1 && min_dis_1[i] > min_dis_2[i])
        	{
        		intersection = ref_of_read_2[i];
        		half_right++;
        		error_read_half << read_name[i] << "\t 2" << "\n";
        	}else if (min_dis_1[i] >= 2 && min_dis_2[i] >= 2)
        	{
        		error++;
        		error_read_all << read_name[i] << "\n";
        		error_read_all_vec.push_back(read_name[i]);
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
        	ref.push_back(" ");
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
    error_read_all.close();
    error_read_half.close();
    return error_read_all_vec;
}

void Analyse(std::vector<std::vector<std::string>> ref_of_reads,std::vector<std::string> true_ref_of_reads,int total,bool flag,std::vector<std::string> error_read_all = std::vector<std::string>(1,"test"))
{
	std::vector<std::string> read_name;
    {
        std::ifstream read_name_file(PAIR_1_NAME_FILE);
        cereal::BinaryInputArchive ar(read_name_file);
        ar(read_name);
    }

    std::ofstream error_read_file;
    if (flag)
    {
    	error_read_file.open("analyse/error_read_sphere.txt");
    }else
    {
    	error_read_file.open("analyse/error_read_bowtie.txt");
    }
	std::vector<std::string> ref;
	std::vector<std::string> error_read;
    size_t count = 0;
    int j = 0;
    for (int i = 0; i < ref_of_reads.size(); ++i)
    {
        ref = ref_of_reads[i];
    //    std::cout << true_ref_of_reads[i] << "-----" << i << std::endl;
        for (j = 0; j < ref.size(); ++j)
        {
        //	std::cout << ref[j] << std::endl;
            if (ref[j].compare(true_ref_of_reads[i]) == 0)
            {
                count++;
                break;
            }
        }
        if (j == ref.size())
        {
        	error_read_file << read_name[i] << "\n";
        	error_read.push_back(read_name[i]);
        }
    }
    std::cout << "correct mapping count:" << (float)count << std::endl;
    std::cout << "total reads:" << total << std::endl;
    std::cout << "correct mapping ratio:" << (float)count / total << std::endl;

    std::ofstream unknown_error("analyse/unknown_errors.txt");
    if (error_read_all.size() != 1)
    {
    	int j = 0;
    	for (int i = 0; i < error_read.size(); ++i)
    	{
    		for (j = 0; j < error_read_all.size(); ++j)
    		{
    			if (error_read_all[j].compare(error_read[i]) == 0)
    			{
    				break;
    			}
    		}
    		if (j == error_read_all.size())
    		{
    			unknown_error << error_read[i]<< std::endl;
    		}
    	}
    }
    unknown_error.close();
}

void Analyse_Result_Spherical()
{
	std::vector<std::string> error_read_all_vec = Save_Ref_Of_Read();
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
    Analyse(ref_of_reads,true_ref_of_reads,true_ref_of_reads.size(),true,error_read_all_vec) ;
}

void Analyse_Result_Others()
{
//	std::ifstream sam("../bowtie2_res/bowtie2_res.sam");
	std::ifstream sam("res/res.sam");
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
	std::vector<std::vector<std::string>> ref_of_reads;

	std::string true_ref_name;
	std::vector<std::string> true_ref_of_reads;

	std::string str;
	int count = 0;
	getline(sam,line);
	count++;
	std::stringstream ss(line);
	ss >> read_name >> str >> ref_name;
	ref_name = ref_name.substr(0,15);
	last_read_name = read_name;
	while(sam.peek() != EOF)
	{
		if(read_name.compare(last_read_name) == 0)
		{
			ref_name_vec.push_back(ref_name);
		//	std::cout << ref_name << std::endl;
		}else
		{
		//	std::cout << std::endl;
			ref_of_reads.push_back(ref_name_vec);
			ref_name_vec.clear();
			ref_name_vec.push_back(ref_name);
		//	std::cout << ref_name << std::endl;

			last_read_name = last_read_name.substr(last_read_name.find_first_of(":") + 1,last_read_name.size());
			last_read_name = last_read_name.substr(last_read_name.find_first_of(":") + 1,last_read_name.size());
			true_ref_name = last_read_name.substr(0,last_read_name.find_first_of(":"));
			true_ref_of_reads.push_back(true_ref_name);
		}
			
		last_read_name = read_name;
		line.clear();
		getline(sam,line);
		count++;
		std::stringstream ss(line);
		ss >> read_name >> str >> ref_name;
		ref_name = ref_name.substr(0,15);
	}
	Analyse(ref_of_reads,true_ref_of_reads,20748,false);
}