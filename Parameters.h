#pragma once

#define	REAL_TYPE				double

// binary code length
#define BCODE_LEN				28		

// number of training samples for spherical hashing
#define NUM_TRAIN_SAMPLES		10000

// desired portion of training set inside of one hyper-sphere
#define INCLUDING_RATIO			0.5
// desired portion of training set inside of two hyper-spheres
#define OVERLAP_RATIO			0.20

// e_m and e_s
#define EPSILON_MEAN			0.10
#define EPSILON_STDDEV			0.15

#define MAX_NUM_ITERATIONS		50

#define INPUT_REF_FILE_NAME		"dataset/training_data.txt"	//

#define TRANSCRIPTS_FILE_NAME	"../reference/transcripts/Homo_sapiens.GRCh38.cdna.all.fa"//"dataset/rapmap_transcripts.fasta"
#define TRANSCRIPTS_STRING_FILE "dataset/true_ref_seq.txt"

#define RAW_READ_FILE_1			"dataset/srrdata/SRR534302_1.fastq"
#define RAW_READ_FILE_2			"dataset/srrdata/SRR534302_2.fastq"

#define USED_READ_FILE_NAME_1	"tmp/tmp_read_used_1.txt"
#define USED_READ_FILE_NAME_2	"tmp/tmp_read_used_2.txt"

#define UNUSED_READ_FILE_NAME_1	"tmp/tmp_read_unused_1.txt"
#define UNUSED_READ_FILE_NAME_2	"tmp/tmp_read_unused_2.txt"

#define MERGE_REF_POS_FILE		"bin/loc_to_ref.bin"
#define MERGE_REF_SEQ_FILE		"bin/merged_ref.bin"
#define MERGE_REF_START_FILE	"bin/merged_ref_start.bin"
#define REF_NAME_FILE			"bin/rname_ref.bin"

#define REDUCED_REGION_CODE_BUCKET_FILE	"bin/code_bucket.bin"
#define REDUCED_REGION_INFO_FILE		"bin/reduced_region_info.bin"

#define DIM						50
#define SAM_FILE_LOC			"res/SRR534302_res.sam"

#define PAIR_1 					1
#define PAIR_2 					2
#define PAIR_1_REGION_FILE		"tmp/tmp_region_1.bin"
#define PAIR_2_REGION_FILE		"tmp/tmp_region_2.bin"
#define PAIR_1_RC_REGION_FILE	"tmp/tmp_rc_region_1.bin"
#define PAIR_2_RC_REGION_FILE	"tmp/tmp_rc_region_2.bin"
#define PAIR_1_RES_REGION_FILE	"res/res_read_region_1.bin"
#define PAIR_2_RES_REGION_FILE	"res/res_read_region_2.bin"

#define PAIR_1_LOC_FILE			"tmp/tmp_loc_1.bin"	
#define PAIR_2_LOC_FILE			"tmp/tmp_loc_2.bin"
#define PAIR_1_DIS_FILE			"tmp/tmp_dis_1.bin"
#define PAIR_2_DIS_FILE			"tmp/tmp_dis_2.bin"

#define PAIR_1_RES_FILE			"res/res_1.txt"
#define PAIR_2_RES_FILE			"res/res_2.txt"
#define PAIR_1_NAME_FILE		"dataset/tmp_name_1.bin"
#define REF_HASH_FILE			"bin/ref_code.bin"


#define CODE_BUFFER_SIZE		200000	//
#define READ_BUFFER_SIZE		50000
#define KMER_SIZE				10
#define SKIP					6
#define TOLERANCE				2
#define IGNORE					5

#define THREAD					4

// to disable parallelization, comment out this
#define USE_PARALLELIZATION		
