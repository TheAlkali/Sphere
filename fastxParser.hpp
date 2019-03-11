#ifndef _FASTX_PARSER_
#define _FASTX_PARSER_


#include <iostream>
#include <fstream>
#include <sstream>//stringstream
#include <vector>
#include <limits>
#include <cstdint>
#include <chrono>
#include <cereal/archives/binary.hpp>

#include "Utils.hpp"

#define READ std::string //falconn::DenseVector<float>//
//#define PAIR_SEQ_LIST std::vector<std::pair<seqInfo, seqInfo> >
//#define SINGLE_SEQ_LIST std::vector<seqInfo>
#define TEST_READ 195
#define TEST_REF 195

std::ifstream seqfile__;
std::ifstream seqfile2__;
std::ofstream out_ref_file__;
std::ofstream out_read_file_1__;
std::ofstream out_read_file_2__;

int klen__;
int kmer_count__;

enum filetype__{DONE,FASTA, FASTQ,ERROR};
enum readtype__{SINGLE,PAIR};

// Adapted from
// https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
// Don't modify the qual
void reverse_complete(std::string& seq,std::string& rev_Read ) {

    rev_Read.resize(seq.length(), 'A');
    int32_t end = seq.length()-1, start = 0;
    //rev_Read[end] = '\0';
    //qualWork[end] = '\0';
    while (start < end) {
        rev_Read[start] = (char)rc_table[(int8_t)seq[end]];
        rev_Read[end] = (char)rc_table[(int8_t)seq[start]];
        ++ start;
        -- end;
    }
    // If odd # of bases, we still have to complement the middle
    if (start == end) {
        rev_Read[start] = (char)rc_table[(int8_t)seq[start]];
        // but don't need to mess with quality
        // qualWork[start] = qual[start];
    }
    //std::swap(seq, rev_Read);
    //std::swap(qual, qualWork);
}

filetype__ get_file_type(std::ifstream &file){
	
	switch(file.peek()) {
	    case EOF: return DONE;
	    case '>': return FASTA;
	    case '@': return FASTQ;
	    default: return ERROR;
	}
	
}

void open_seqfile(std::string filename){
	seqfile__.open(filename,std::ifstream::binary);//|std::fstream::in);
	if(!seqfile__.is_open()){
		std::cout<<"fail to open fasta/fastq files:" << filename << std::endl;
	}
}

void open_seqfile(std::string filename1,std::string filename2){
	seqfile__.open(filename1,std::ifstream::binary);//,std::fstream::in);
	seqfile2__.open(filename2,std::ifstream::binary);//,std::fstream::in);
	if(!seqfile__.is_open()){
		std::cerr<<"fail to open fasta/fastq files:" << filename1 << std::endl;
	}
	if(!seqfile2__.is_open()){
		std::cerr<<"fail to open fasta/fastq files:" << filename2 << std::endl;
	}
}

void read_fq_oneseq(std::ifstream &seqfile,std::string &name, READ &seq, std::string &qual){
	std::string line;
//	std::string seq;
	//ignore @
	char tmp;

	seqfile.get(tmp);
	std::getline(seqfile,line);
	//use sstream to split the name and description
	std::stringstream ss(line);
	name.clear();
	ss >> name;

//	seqint.clear();
	while(seqfile.peek() != '+' && seqfile.peek() != EOF){	
		std::getline(seqfile,line);
		seq.append(line);
	}
	
//	trans_seq2array(seq,seqint);

	//ignore the description
	seqfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	qual.clear();

	while(qual.size() < seq.length() && seqfile.good()) {
	      std::getline(seqfile, line);
	      qual.append(line);
	}
}

void  read_fa_oneseq(std::ifstream &seqfile,std::string &name, READ &seq){
	std::string line;
//	std::string seq;
	int len;
	//ignore >
	seqfile.get();
	std::getline(seqfile,line);
	std::stringstream ss(line);
	name.clear();
	ss >> name;

//	seqint.clear();
	while(seqfile.peek() != '>' && seqfile.peek() != EOF){
		std::getline(seqfile,line);
		seq.append(line);
	} 
//	trans_seq2array(seq,seqint);
}

void get_one_kmer(READ seq,int loc, bool if_out){
	READ kmer;
	int ik = 0;
	kmer.resize(klen__);
	for (int i = loc; i < loc + klen__; ++i){
		kmer[ik] = seq[i];
		ik++;
		if (if_out){
			out_ref_file__ << (char)stoic_table[(int8_t)seq[i]];
		}
		
	}
	if (if_out){
		out_ref_file__ << std::endl;	
	}
	
	kmer_count__++;																																																																																																																																																																																																																																																										
}

std::string merge_ref_seq(std::string filename,int seg_len){
	open_seqfile(filename);
	filetype__ filetype = get_file_type(seqfile__);

	size_t size;
	std::string name;
	std::string qual;
	std::string seq,all_seq;
	int len = 0;

	ref_start_name rpos;

	std::ofstream train_data(INPUT_REF_FILE_NAME);
	int train_num = 0;
	int start = KMER_SIZE;
	int cond = BCODE_LEN * seg_len + KMER_SIZE;

	auto t1 = std::chrono::high_resolution_clock::now();
	for (size  = 0;seqfile__.peek() != EOF;size++){
		seq.clear();
		read_fa_oneseq(seqfile__,name,seq);
		// generate training data
		if(seq.length() > 2 * DIM)
		{
			if (train_num < NUM_TRAIN_SAMPLES)
			{
				train_num++;
				int rand_loc = rand() % (seq.length() - DIM + 1);
				std::string train_read = seq.substr(rand_loc,DIM);
				
				while(true)
				{
					if (start < cond)
					{
						for (int j = 0; j < seg_len; ++j)
						{
							train_data << stoic_table[(int8_t)train_read[start + j]];
						//	std::cout << d[k][j];
						}
						train_data << std::endl;
						//std::cout << std::endl;
						start += seg_len;
					}else
					{
						start = KMER_SIZE;
						break;
					}		
				}
			}

			rpos.rname.push_back(name);
			rpos.ref_start.push_back(len);
			len += seq.size() + 1;
			for (int i = 0; i < seq.size(); ++i)
			{
				all_seq.push_back((char)stoic_table[(int8_t)seq[i]]);
			}
			all_seq.append("9");
		}
		
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	double fastxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	std::cout << "- Total Transcriptomes:" << size << std::endl;
	std::cout << "- Time Of Merging Ref Fasta (" << fastxparser_time << " seconds)" << std::endl;

	// save merge seq to disk
/*	{
		std::ofstream ref_seq("tmp/true_ref_seq.bin");
		cereal::BinaryOutputArchive ar(ref_seq);
		ar(CEREAL_NVP(all_seq));
	}*/
	// save ref position in merged seq to disk
	t1 = std::chrono::high_resolution_clock::now();
	{
		std::ofstream ref_start_name(MERGE_REF_POS_FILE);
		cereal::BinaryOutputArchive ar(ref_start_name);
		ar(CEREAL_NVP(rpos.ref_start),CEREAL_NVP(rpos.rname));
	}
	t2 = std::chrono::high_resolution_clock::now();
	fastxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	std::cout << "- Time Of Saving Reference Position In Merged Ref:" << fastxparser_time << std::endl;
	seqfile__.close();
	train_data.close();
	return all_seq;
}

// directly store the ref kmers to the file
void store_ref_kmers(std::string filename, int klen){
	open_seqfile(filename);
	filetype__ filetype = get_file_type(seqfile__);

	size_t size;
	std::string name;
	std::string qual;
	READ seq;
	int len;

	klen__ = klen;
	kmer_count__ = 0;

	out_ref_file__.open("data.txt",std::ofstream::binary|std::ofstream::out);
	if (!out_ref_file__.is_open())
	{
		std::cerr << "can't open out_ref_file: data.txt" << std::endl;
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	for (size  = 0;seqfile__.peek() != EOF;size++){
//	for (size = 0;size < TEST_REF;size++){
		read_fa_oneseq(seqfile__,name,seq);
		std::cout << "parse transcriptome:" << name << std::endl;
		len = seq.size();
		for (int loc = 0;loc < (len - klen + 1);loc++){	
			get_one_kmer(seq,loc,true);
		}
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	double fasxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	std::cout << "total transcriptomes:" << size << std::endl;
	std::cout << "total kmers:" << kmer_count__ << std::endl;
	std::cout << "parser time of fasta:" << fasxparser_time << std::endl;

	seqfile__.close();
	out_ref_file__.close();
}

// TODO delete
void store_reads_name(std::string filename1){
	open_seqfile(filename1);
	std::vector<std::string> name_vec_1;
	std::string name_1;
	std::string fqual;
	READ first;
	filetype__ filetype = get_file_type(seqfile__);
	if (filetype == FASTQ){
		while(seqfile__.peek() != EOF){
	//	for (;size < TEST_READ;size++){
			first.clear();
			name_1.clear();
			read_fq_oneseq(seqfile__,name_1,first,fqual);

			name_vec_1.push_back(name_1);
		}
	}
	std::cout << "- store reads name finshed "<< std::endl;
	seqfile__.close();

	//store reads name
	{
		std::ofstream reads_name_1(PAIR_1_NAME_FILE);
		cereal::BinaryOutputArchive ar(reads_name_1);
		ar(name_vec_1);
	}
}

//directly store the reads to the files
size_t store_reads(){
	open_seqfile(RAW_READ_FILE_1,RAW_READ_FILE_2);

	klen__ = DIM;

	out_read_file_1__.open(INPUT_READ_FILE_NAME_1,std::ofstream::out);
	out_read_file_2__.open(INPUT_READ_FILE_NAME_2,std::ofstream::out);
	if (!out_read_file_1__.is_open()){
		std::cerr << "can't open out_read_file:" << INPUT_READ_FILE_NAME_1 << std::endl;
	}
	if (!out_read_file_2__.is_open()){
		std::cerr << "can't open out_read_file:" << INPUT_READ_FILE_NAME_2 << std::endl;
	}

	size_t size  = 0;
	std::vector<std::string> name_vec_1;
	std::string name_1;
	std::string name_2;
	std::string fqual;
	std::string squal;
	READ first;
	READ second;
	std::pair<READ,READ> seqpair;
	std::pair<std::string,std::string> qualpair;

	filetype__ filetype = get_file_type(seqfile__);
	filetype__ filetype2 = get_file_type(seqfile2__);
	if (filetype == FASTQ && filetype2 == FASTQ){
		for (;seqfile__.peek() != EOF && seqfile2__.peek() != EOF;size++){
	//	for (;size < TEST_READ;size++){
			first.clear();
			second.clear();
			name_1.clear();
			name_2.clear();
			read_fq_oneseq(seqfile__,name_1,first,fqual);
			read_fq_oneseq(seqfile2__,name_2,second,squal);

			name_vec_1.push_back(name_1);

			std::string rc_first;
			reverse_complete(first,rc_first);
			for (int i = 0; i < klen__; ++i){
				out_read_file_1__ << rc_first[i];
			}
			out_read_file_1__ << std::endl;
		
			// pair 2 reverse complete
			
			for (int i = 0; i < klen__; ++i){
				out_read_file_2__ << stoic_table[(int8_t)second[i]];
			}
			out_read_file_2__ << std::endl;
		}
	}
	std::cout << "- total reads:" << size << std::endl;
	seqfile__.close();
	seqfile2__.close();
	out_read_file_1__.close();
	out_read_file_2__.close();

	//store reads name
	{
		std::ofstream reads_name_1(PAIR_1_NAME_FILE);
		cereal::BinaryOutputArchive ar(reads_name_1);
		ar(name_vec_1);
	}
	return size;
}

//directly store the reads to the files
void store_reads(char* filename,int klen){
	open_seqfile(filename);

	klen__ = klen;

	out_read_file_1__.open("dataset/ref_sim.txt",std::ofstream::out);
	if (!out_read_file_1__.is_open()){
		std::cerr << "can't open out_read_file: read_1.txt" << std::endl;
	}

	size_t size  = 0;
	std::vector<std::string> name_vec_1;
	std::string name_1;
	std::string fqual;
	std::string squal;
	READ first;
	std::pair<READ,READ> seqpair;
	std::pair<std::string,std::string> qualpair;

	filetype__ filetype = get_file_type(seqfile__);
	if (filetype == FASTQ){
		for (;seqfile__.peek() != EOF;size++){
	//	for (;size < TEST_READ;size++){
			first.clear();
			name_1.clear();
			read_fq_oneseq(seqfile__,name_1,first,fqual);

			name_vec_1.push_back(name_1);

			for (int i = 0; i < klen__; ++i){
				out_read_file_1__ << (char)stoic_table[(int8_t)first[i]];
			}
			out_read_file_1__ << std::endl;
		}
	}
	std::cout << "- total reads:" << size << std::endl;
	seqfile__.close();
	out_read_file_1__.close();

	//store reads name
	{
		std::ofstream reads_name_1(PAIR_1_NAME_FILE);
		cereal::BinaryOutputArchive ar(reads_name_1);
		ar(name_vec_1);
	}
}

//store data to RAM
void read_fastx(char* filename,singleSeqList *seqlist){	
	open_seqfile(filename);
	filetype__ filetype = get_file_type(seqfile__);

	size_t size;
	std::string name;
	std::string qual;
	READ seq;
	if (filetype == FASTQ){

		for (;seqfile__.peek() != EOF;size++){
	//	for(size = 0;size < TEST_READ;size++){
			read_fq_oneseq(seqfile__,name,seq,qual);
			seqlist->seq.push_back(std::move(seq));
			seqlist->name.push_back(std::move(name));
			seqlist->qual.push_back(std::move(qual));
		}
		seqfile__.close();
	}else if(filetype == FASTA){
		auto t1 = std::chrono::high_resolution_clock::now();
		for (size  = 0;seqfile__.peek() != EOF;size++){
	//	for (size = 0;size < TEST_REF;size++){
			read_fa_oneseq(seqfile__,name,seq);
			seqlist->seq.push_back(std::move(seq));
			seqlist->name.push_back(std::move(name));
	/*		for (int i = 0; i <  seq.size(); ++i){
				std::cout << seq[i];
			}
			std::cout << std::endl;*/
		}
		auto t2 = std::chrono::high_resolution_clock::now();
		double fasxparser_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
		std::cout << "Fastx Parser Time:" << fasxparser_time << std::endl;
	}
	seqfile__.close();
}

//store data to RAM
void read_fastx(char* filename1, char* filename2,pairSeqList *seqlist){	
	open_seqfile(filename1,filename2);

	size_t size  = 0;
	std::string name;
	std::string fqual;
	std::string squal;
	READ first;
	READ second;
	std::pair<READ,READ> seqpair;
	std::pair<std::string,std::string> qualpair;

	filetype__ filetype = get_file_type(seqfile__);
	filetype__ filetype2 = get_file_type(seqfile2__);
	if (filetype == FASTQ && filetype2 == FASTQ){
		for (;seqfile__.peek() != EOF && seqfile2__.peek() != EOF;size++){
	//	for (;size < TEST_READ;size++){
			read_fq_oneseq(seqfile__,name,first,fqual);
			read_fq_oneseq(seqfile2__,name,second,squal);
			seqpair.first = first;
			seqpair.second = second;
			qualpair.first = fqual;
			qualpair.second = squal;
			seqlist->seq.push_back(std::move(seqpair));
			seqlist->name.push_back(std::move(name));
			seqlist->qual.push_back(std::move(qualpair));

		}
	}else if(filetype == FASTA && filetype2 == FASTA){
		for (size  = 0;seqfile__.peek() != EOF && seqfile2__.peek() != EOF;size++){
	//	for (size = 0;size < TEST_REF;size++){
			read_fa_oneseq(seqfile__,name,seqpair.first);
			read_fa_oneseq(seqfile2__,name,seqpair.second);
			seqlist->name.push_back(std::move(name));
			seqlist->seq.push_back(std::move(seqpair));
		}
	}else{
		std::cerr << "the types of two files are not the same. " << std::endl;
		seqfile__.close();
		seqfile2__.close();
		return ;
	}
	seqfile__.close();
	seqfile2__.close();
}

#endif