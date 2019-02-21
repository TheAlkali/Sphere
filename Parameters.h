#pragma once

#define	REAL_TYPE				float


// target number of nearest neighbors
#define KNN						20

// binary code length
#define BCODE_LEN				16		

// number of training samples for spherical hashing
#define NUM_TRAIN_SAMPLES		1000

// desired portion of training set inside of one hyper-sphere
#define INCLUDING_RATIO			0.5
// desired portion of training set inside of two hyper-spheres
#define OVERLAP_RATIO			0.20

// e_m and e_s
#define EPSILON_MEAN			0.10
#define EPSILON_STDDEV			0.15

#define MAX_NUM_ITERATIONS		50

#define INPUT_REF_FILE_NAME		"dataset/ref238956.txt"	//
#define TRANSCRIPTS_FILE_NAME	"../reference/transcripts/Homo_sapiens.GRCh38.cdna.all.fa"//"dataset/rapmap_transcripts.fasta"
#define TRANSCRIPTS_STRING_FILE "dataset/srrdata/true_ref_seq.txt"//"dataset/all_ref_seq.txt"

#define RAW_READ_FILE_1			"dataset/srrdata/small_SRR5337025_1.fastq"
#define RAW_READ_FILE_2			"dataset/srrdata/small_SRR5337025_2.fastq"

#define INPUT_READ_FILE_NAME_1	"dataset/srrdata/small_srr25_1.txt"//"dataset/sim_mis_read_5_1.txt"
#define INPUT_READ_FILE_NAME_2	"dataset/srrdata/small_srr25_2.txt"//"dataset/sim_mis_read_5_2.txt"

#define MERGE_REF_POS_FILE		"dataset/merge_ref_pos.bin"
#define DIM						50
#define SAM_FILE_LOC			"res/res.sam"

#define PAIR_1 					1
#define PAIR_2 					2
#define PAIR_1_LOC_FILE			"tmp/tmp_loc_1.bin"	
#define PAIR_2_LOC_FILE			"tmp/tmp_loc_2.bin"
#define PAIR_1_DIS_FILE			"tmp/tmp_dis_1.bin"
#define PAIR_2_DIS_FILE			"tmp/tmp_dis_2.bin"
#define PAIR_1_RES_FILE			"res/res_1.txt"
#define PAIR_2_RES_FILE			"res/res_2.txt"
#define PAIR_1_NAME_FILE		"dataset/tmp_name_1.bin"
//#define OUTPUT_REF_HASH_FILE_NAME	"bin/ref_code_tmp.txt"
//#define OUTPUT_READ_HASH_FILE_NAME  "bin/read_code_tmp.txt"

#define REF_BUFFER_SIZE			200000	//
#define READ_BUFFER_SIZE		1000
#define KMER_SIZE				10	
#define TOLERANCE				7

// to disable parallelization, comment out this
#define USE_PARALLELIZATION		
