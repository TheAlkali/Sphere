#include<iostream>
#include<string.h>
#include<time.h>
//#include"kmerUtils.hpp"
#include"fastxParser.hpp"
#include"SAMparser.hpp"
#include "Utils.hpp"

//#include"benchmark.hpp"

void gen_simulate_reads(char*,char*, singleSeqList*,int);
void gen_simulate_reads_From_TXT(int,int);
void analyse_result(std::string,int);
void count_ref_kmer(std::string);

//int main(int argc, char const *argv[])
int main()
{
//	kmerUtils kmutil;
	int klen = 76;
//	LSH lsh;
	int gate = 0;//std::atoi(argv[1]);

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
		char test_sam[] = "../rapmap_res/rapmap_mapped_reads.sam";
		char test_res[] = "rapmap_res.txt";
		sp.get_Sam(ref,test_sam,test_res,20000,klen);
	}else if (gate == 1){		
		store_reads();
//		store_reads("../reads/single-cell-data/SRR5337/Simulate_GRCh38.1.fa",klen);
//		store_ref_kmers(ref,klen);
	}else if (gate == 2){		
		gen_simulate_reads(read1,ref,reflist,klen);
		gen_simulate_reads(read2,ref,reflist,klen);
	}else if (gate == 3){
		//merge all refs to one
		merge_ref_seq(ref);
	}else if (gate == 4){
	/*	fparser.read_fastx(read1,read2,preadlist);
		kmutil.set_kmer_length(klen);
		kmutil.genkmers_from_all_transcripts(reflist);
		delete reflist;
		refklist = kmutil.get_kmer_list();*/
	}else if (gate == 5){
		gen_simulate_reads_From_TXT(238956,klen);
	}else if (gate == 6){
		analyse_result("res_1.txt",klen);
	}else if (gate == 7){
		//ref kmer count:282626422
		count_ref_kmer("dataset/ref_kmer.txt");
	}else if (gate == 8){
		store_reads_name(read1);
	}

	return 0;
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

void analyse_result(std::string filename,int dim){
	std::ifstream res;
	res.open(filename);
	std::string aread;
	char tmp;
	std::string seq;
	std::string true_read;

	int all_right = 0;
	int some_error = 0;
	int no_right = 0;
	bool right = false;
	bool error = false;
	int i = 1;
	while(res.peek() != EOF){
		res.get(tmp);
		if (tmp == '='){
			aread.clear();
			if (right && error){
				some_error++;
			}else if (right && !error){
				all_right++;
			}else if (!right && error){
				no_right++;
			}
			getline(res, aread);
			aread.erase(0,aread.find_first_not_of(" "));
			right = false;
			error = false;
			if (i % 2 != 0){
				true_read = aread;
			}else{
				aread = true_read;
			}
			i++;
		}else if(tmp == '+'){
			seq.clear();
			getline(res, seq);
			std::stringstream ss(seq);
			ss >> seq;
			if (hamming_dis(true_read,seq,dim)){
				right = true;
			}else{
				error = true;
			}
		}else{
			getline(res,seq);
		}
	}

	std::cout << "all right: " << all_right << "/" << (float)all_right / i << std::endl;
	std::cout << "including error: " << some_error << "/" << (float)some_error / i << std::endl;
	std::cout << "all wrong: " << no_right << "/" << (float)no_right / i<< std::endl;

	std::cout << "right rate:" << all_right + some_error << "/" << (float)(all_right + some_error) / i << std::endl;
	res.close();
}
