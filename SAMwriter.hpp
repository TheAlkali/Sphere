#include <iostream>
#include <fstream>

#include "algorithm"
#include "Mapping.hpp"
/*#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/basic_file_sink.h"*/

struct res_analysis
{
    int ref_of_read;
    int pos_of_ref_of_read_1;
};

struct res_analysis_intersection:res_analysis
{
    int pos_of_ref_of_read_2;

    void equal_to_pair_1(res_analysis & rai)
    {
    	ref_of_read = std::move(rai.ref_of_read);
    	pos_of_ref_of_read_1 = std::move(rai.pos_of_ref_of_read_1);
    	pos_of_ref_of_read_2 = 0;
    }
    void equal_to_pair_2(res_analysis & rai)
    {
    	ref_of_read = std::move(rai.ref_of_read);
    	pos_of_ref_of_read_2 = std::move(rai.pos_of_ref_of_read_1);
    	pos_of_ref_of_read_1 = 0;
    }
};

struct SAM_format{
	std::string qname;
	std::string rname;
	std::string cigar;
	std::string rnext;
	std::string seq;
	std::string qual;
	int flag;
	int pos;
	int mapq;
	int pnext;
	int tlen;
	template<typename OStream>
    friend OStream &operator<<(OStream &os, const SAM_format &s){	
    	os << s.qname << '\t' // QNAME
			<< s.flag << '\t' // FLAGS
			<< s.rname << '\t' // RNAME
			<< s.pos << '\t' // pos (1-based)
			<< 255 << '\t' // MAPQ
			<< "50M" << '\t' // CIGAR
			<< '=' << '\t' // MATE NAME
			<< s.pnext << '\t' // MATE pos
			<< s.tlen << '\t' // TLEN
			<< s.seq << '\t' // SEQ
			<< "*\t" << '\n';// QSTR
    }

};

class SAMwriter
{
private:
//	std::shared_ptr<spdlog::logger> samlog;
	std::vector<int> loc_to_ref;
	std::vector<size_t> ref_start;

	std::vector<std::vector<res_analysis>> ra_of_read_1;
	std::vector<std::vector<res_analysis>> ra_of_read_2;
	std::vector<int> mres_1_min_dis,mres_2_min_dis;

    region_profile rpro;
	std::vector<std::vector<int>> ref_of_merged_res;
	std::vector<std::vector<int>> pos_1_of_ref_of_merged_res;
	std::vector<std::vector<int>> pos_2_of_ref_of_merged_res;
	std::vector<int> res_from_which_pair;

	REAL_TYPE **read_1_buff;
    REAL_TYPE **read_2_buff;

    std::vector<bool> is_read_1_rev;
    std::vector<bool> is_read_2_rev;
	int dim;

	int flag_1[4] = {4,8,8,2};
	int flag_2[2] = {32,16};

public:
	SAMwriter(int dim)
	{
		this->dim = dim;
	};

	friend std::string Mapping::Get_Read_1_Buff(int idx);
	friend std::string Mapping::Get_Read_2_Buff(int idx);
	friend bool Mapping::Get_Is_Read_1_Rev(int idx);
	friend bool Mapping::Get_Is_Read_2_Rev(int idx);

	void Transfer_Info_From_Mapping(std::vector<bool> &is_read_1_rev,std::vector<bool> &is_read_2_rev,region_profile &rpro,REAL_TYPE **read_1_buff,REAL_TYPE **read_2_buff)
	{
		this->is_read_1_rev = std::move(is_read_1_rev);
		this->is_read_2_rev = std::move(is_read_2_rev);
		this->rpro = std::move(rpro);
		this->read_1_buff = &read_1_buff[0];
		this->read_2_buff = &read_2_buff[0];
	}

	static bool greater_comp(res_analysis a,res_analysis b){
    	return a.ref_of_read < b.ref_of_read;
	}

	void Load_Info()
	{
		Stopwatch T0("");
        T0.Reset();     T0.Start();

        {
            std::ifstream loc_to_ref_file(MERGE_REF_POS_FILE);
            cereal::BinaryInputArchive ar_ref_pos(loc_to_ref_file);
            ar_ref_pos(loc_to_ref);
        }

        {
        	std::ifstream ref_start_file(MERGE_REF_START_FILE);
        	cereal::BinaryInputArchive ar(ref_start_file);
        	ar(ref_start);
        }
        T0.Stop();
        printf("- Load Ref Info Finished (%f seconds)\n",T0.GetTime() );
	}

	void Analyse_Result(int pair,std::vector<size_t> &code_bucket)
    {
        std::ifstream tmp_loc;
        std::ifstream tmp_dis;

        std::ifstream read_region_file;
        if (pair == PAIR_1)
        {
            tmp_loc.open(PAIR_1_LOC_FILE);
            tmp_dis.open(PAIR_1_DIS_FILE);
	        
            read_region_file.open(PAIR_1_RES_REGION_FILE);
        }else if (pair == PAIR_2)
        {
            tmp_loc.open(PAIR_2_LOC_FILE);
            tmp_dis.open(PAIR_2_DIS_FILE);
            read_region_file.open(PAIR_2_RES_REGION_FILE);
        }

        mapped_res mres;
        // load mapping result from disk using cereal
        {
            cereal::BinaryInputArchive ar_loc(tmp_loc);
            cereal::BinaryInputArchive ar_dis(tmp_dis);
            ar_loc(CEREAL_NVP(mres.min_code_idx));
            ar_dis(CEREAL_NVP(mres.min_dis));
        };

    	std::vector<uint64_t> read_region;
        {
            cereal::BinaryInputArchive ar_read(read_region_file);
            ar_read(read_region);
        }

        Stopwatch T0("");
        T0.Reset();     T0.Start();

        size_t ref_loc = 0;
        int loc_size = 0;
        int read_size = mres.min_code_idx.size();
        std::vector<std::vector<res_analysis>> ra_of_read(read_size);
        std::vector<res_analysis> ra_vec;
        res_analysis ra;

    #ifdef USE_PARALLELIZATION
        #pragma omp parallel for
    #endif
        for (unsigned int qIndex = 0;qIndex < read_size;++qIndex)
        {      
            if (mres.min_code_idx[qIndex][0] >= 0)
            {
            	for (int j = 0; j < mres.min_code_idx[qIndex].size(); ++j)
            	{
	                if (mres.min_code_idx[qIndex][j] == rpro.code_bucket_idx[read_region[qIndex]].size() - 1)
	                {
	                	if (read_region[qIndex] + 1 < rpro.code_bucket_idx.size() - 1)
	                    {
	                    	loc_size = rpro.code_bucket_idx[read_region[qIndex] + 1][0] 
	                        	- rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j]];
	                    }else
	                    {
	                    	loc_size  = 1;
	                    }
	                }else
	                {
	                    loc_size = rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j] + 1] 
	                        - rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j]];
	                }
	                for (unsigned int i = 0; i < loc_size; ++i)
	                {
	                    ref_loc = code_bucket[rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex][j]] + i];
	                    if (ref_loc > 0)
	                    {
	                    	ra.ref_of_read = loc_to_ref[ref_loc];
	                    	ra.pos_of_ref_of_read_1 = ref_loc - ref_start[loc_to_ref[ref_loc]] + 1;
	                        ra_vec.push_back(ra);
	                    }
	                }
            	}
                ra_of_read[qIndex] = ra_vec;
                ra_vec.clear();
            }
        }

        T0.Stop();
        if (pair == PAIR_1)
        {
            ra_of_read_1 = std::move(ra_of_read);
            mres_1_min_dis = std::move(mres.min_dis);
        	printf("- Analyse Result Of Pair 1 Finished (%f seconds)\n",T0.GetTime() );
        }else if (pair == PAIR_2)
        {
            ra_of_read_2 = std::move(ra_of_read);
            mres_2_min_dis = std::move(mres.min_dis);
            printf("- Analyse Result Of Pair 2 Finished (%f seconds)\n",T0.GetTime() );
        }
    }

    void Analyse_Result_Pair()
    {
    	std::vector<size_t> code_bucket;

    	Stopwatch T0("");
    	T0.Reset();     T0.Start();
        std::ifstream bucket_file("bin/code_bucket.bin");
        bucket_file.seekg(0,bucket_file.end);
        int bucket_size = bucket_file.tellg() / sizeof(size_t);
        code_bucket.resize(bucket_size);
        bucket_file.seekg(0,bucket_file.beg);
        bucket_file.read(reinterpret_cast<char*>(&code_bucket[0]),bucket_size * sizeof(size_t));
        T0.Stop();
		printf("- Load Code Bucket Info Finished (%f seconds)\n",T0.GetTime() );

		Analyse_Result(PAIR_1,code_bucket);
		Analyse_Result(PAIR_2,code_bucket);
    }

    void res_intersection(std::vector<res_analysis> ra_of_read_1,std::vector<res_analysis> ra_of_read_2,std::vector<res_analysis_intersection> &intersection)
    {
    	res_analysis_intersection rai;
    	int idx = 0;
    	for (int i = 0; i < ra_of_read_1.size(); ++i)
    	{
    		for (; idx < ra_of_read_2.size(); ++idx)
    		{
    			if (ra_of_read_1[i].ref_of_read == ra_of_read_2[idx].ref_of_read)
    			{
    				rai.ref_of_read = ra_of_read_1[i].ref_of_read;
    				rai.pos_of_ref_of_read_1 = ra_of_read_1[i].pos_of_ref_of_read_1;
    				rai.pos_of_ref_of_read_2 = ra_of_read_2[idx].pos_of_ref_of_read_1;
    				intersection.push_back(rai);
    			}else if (ra_of_read_1[i].ref_of_read < ra_of_read_2[idx].ref_of_read)
    			{
    				break;
    			}
    		}
    	}
    }
	void Merge_Result()
	{	
		Stopwatch T0("");
        T0.Reset();     T0.Start();

	    std::vector<res_analysis_intersection> intersection;
	    std::vector<res_analysis_intersection>::iterator it;
	    std::vector<int> pos_1;
	    std::vector<int> pos_2;
	    std::vector<int> ref;

	    int read_size = ra_of_read_1.size();

	    int half_right = 0;
	    int error = 0;

	#ifdef USE_PARALLELIZATION
        #pragma omp parallel for
    #endif
	    for (int i = 0; i < read_size; ++i)
	    {
	    	ref.clear();
	        pos_1.clear();
	        pos_2.clear();
	        intersection.clear();

	        sort(ra_of_read_1[i].begin(),ra_of_read_1[i].end(),greater_comp);
	        sort(ra_of_read_2[i].begin(),ra_of_read_2[i].end(),greater_comp);
	        res_intersection(ra_of_read_1[i],ra_of_read_2[i],intersection);

	        if (intersection.size() == 0)
	        {
	        	intersection.clear();
	        	if (mres_1_min_dis[i] <= filter && mres_2_min_dis[i] >= mres_1_min_dis[i])
	        	{
	        		intersection.resize(ra_of_read_1[i].size());
	        		for (int j = 0; j < ra_of_read_1[i].size(); ++j)
	        		{
	        			intersection[j].equal_to_pair_1(ra_of_read_1[i][j]);

	        		}
	        		half_right++;
	        		res_from_which_pair.push_back(1);
	        	}else if (mres_2_min_dis[i] <= filter && mres_1_min_dis[i] >= mres_2_min_dis[i])
	        	{
	        		intersection.resize(ra_of_read_2[i].size());
	        		for (int j = 0; j < ra_of_read_2[i].size(); ++j)
	        		{
	        			intersection[j].equal_to_pair_2(ra_of_read_2[i][j]);
	        		}
	        		half_right++;
	        		res_from_which_pair.push_back(2);
	        	}else if (mres_1_min_dis[i] > filter && mres_2_min_dis[i] > filter)
	        	{
	        		error++;
	        		res_from_which_pair.push_back(0);
	        		pos_1.push_back(0);
	        		pos_2.push_back(0);
	        		// TODO
	        		ref.push_back(-1);
	        	}
	        }else
	        {
	        	res_from_which_pair.push_back(3);
	        }
	        for (it = intersection.begin(); it != intersection.end(); ++it)
	        {
	            ref.push_back(std::move((*it).ref_of_read));
	        	pos_1.push_back(std::move((*it).pos_of_ref_of_read_1));
	        	pos_2.push_back(std::move((*it).pos_of_ref_of_read_2));
	        }
	        

        	ref_of_merged_res.push_back(std::move(ref));
        	pos_1_of_ref_of_merged_res.push_back(std::move(pos_1));
        	pos_2_of_ref_of_merged_res.push_back(std::move(pos_2));
	    }
	    T0.Stop();

	    printf("- Merge Result Finished (%f seconds)\n",T0.GetTime() );

	    std::cout << "- both  hamming distances are larger than "<< filter << ":" << error << std::endl;
	    std::cout << "- half mapping:" << half_right << std::endl;
	}


	void Set_Seq_Of_SAM(bool is_read_rev,SAM_format &read,REAL_TYPE *seq)
	{
		std::string rev_read;
		std::string seq_str;
		seq_str.resize(dim);
		for (int i = 0; i < dim; ++i)
		{
			seq_str[i] = itos_table[(int8_t)seq[i]];
		}
		if (is_read_rev)
        {
            reverse_complete(seq_str,rev_read);
            read.seq = rev_read;
        }else
        {
        	read.seq.clear();
        	read.seq = seq_str;
        }	
	}
	
	void Set_Flag(SAM_format &read,int idx,int is_rc)
	{	
		read.flag = 0;
		read.flag += flag_1[res_from_which_pair[idx]];

		if (read.tlen < 0)
		{
			read.flag += 64;
		}else if (read.tlen > 0)
		{
			read.flag += 128;
		}

		read.flag += flag_2[is_rc];
		if (is_rc)
		{
			read.flag += 16;
		}else
		{
			read.flag += 32;
		}
	}

	void Generate_SAM()
	{
		std::ofstream samfile(SAM_FILE_LOC);

		std::vector<std::string> read_name;
	    {
	        std::ifstream read_name_file(PAIR_1_NAME_FILE);
	        cereal::BinaryInputArchive ar(read_name_file);
	        ar(read_name);
	    }

	    std::vector<std::string> ref_name;
		{
	        std::ifstream ref_name_file(REF_NAME_FILE);
	        cereal::BinaryInputArchive ar_ref_name(ref_name_file);
	        ar_ref_name(ref_name);
	    }

	    SAM_format read,next_read;

	    int read_size = read_name.size();
	    std::vector<int> ref;
	    std::string dim_str = std::to_string(dim);
	    std::string rev_read;
	#ifdef USE_PARALLELIZATION
        #pragma omp parallel for
    #endif
	    for (int i = 0; i < read_size; ++i)
	    {
	    	read.qname = std::move(read_name[i]);
	    	read.tlen = dim;
	    //	read.cigar = dim_str + "M";

	    	if (res_from_which_pair[i] == 1)
	    	{
	    		read.tlen = 0;
	    		read.pnext = 0;
		    	Set_Seq_Of_SAM(is_read_1_rev[i],read,read_1_buff[i]);
		    	Set_Flag(read,i,is_read_1_rev[i]);
		    	ref = std::move(ref_of_merged_res[i]);
		    	for (int j = 0; j < ref.size(); ++j)
		    	{
	    			read.pos = pos_1_of_ref_of_merged_res[i][j];
		    		if (ref[j] >= 0)
		    		{
		    			read.rname = ref_name[ref[j]];
		    		}else
		    		{
		    			read.rname = "*";
		    		}
		    		samfile << read;
		    	}
	    	}else if (res_from_which_pair[i] == 2)
	    	{
	    		read.tlen = 0;
	    		read.pnext = 0;
	    		Set_Seq_Of_SAM(is_read_2_rev[i],read,read_2_buff[i]);
	    		Set_Flag(read,i,is_read_2_rev[i]);
	    		ref = std::move(ref_of_merged_res[i]);
		    	for (int j = 0; j < ref.size(); ++j)
		    	{
	    			read.pos = pos_2_of_ref_of_merged_res[i][j];
		    		if (ref[j] >= 0)
		    		{
		    			read.rname = ref_name[ref[j]];
		    		}else
		    		{
		    			read.rname = "*";
		    		}
		    		samfile << read;
		    	}
	    	}else if (res_from_which_pair[i] == 3 || res_from_which_pair[i] == 0)
	    	{
	    		
	    		next_read.qname = read.qname;
		    	next_read.cigar = std::to_string(dim) + "M";

	    		Set_Seq_Of_SAM(is_read_1_rev[i],read,read_1_buff[i]);
	    		Set_Seq_Of_SAM(is_read_2_rev[i],next_read,read_2_buff[i]);
	    		ref = std::move(ref_of_merged_res[i]);
		    	for (int j = 0; j < ref.size(); ++j)
		    	{
		    		read.pos = pos_1_of_ref_of_merged_res[i][j];
		    		next_read.pos = pos_2_of_ref_of_merged_res[i][j];
		    		read.pnext = next_read.pos;
		    		next_read.pnext = read.pos;

		    		read.tlen = read.pos - next_read.pos + dim;
		    		next_read.tlen = -read.tlen;

		    		Set_Flag(read,i,is_read_1_rev[i]);
		    		Set_Flag(next_read,i,is_read_2_rev[i]);

		    		if (ref[j] >= 0)
		    		{
		    			read.rname = ref_name[ref[j]];
		    			next_read.rname = ref_name[ref[j]];
		    		}else
		    		{
		    			read.rname = "*";
		    			next_read.rname = "*";
		    		}
		    		samfile << read;
		    		samfile << next_read;
		    	}
	    	}
	    }
	    samfile.close();
	}
};