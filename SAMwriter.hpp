#include <iostream>
#include <fstream>
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/basic_file_sink.h"

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
    	os << '*' << '\t' // QNAME
			<< '*' << '\t' // FLAGS
			<< s.rname << '\t' // RNAME
			<< s.pos + 1 << '\t' // POS (1-based)
			<< 255 << '\t' // MAPQ
			<< '*' << '\t' // CIGAR
			<< '*' << '\t' // MATE NAME
			<< 0 << '\t' // MATE POS
			<< s.tlen << '\t' // TLEN
			<< '*' << '\t' // SEQ
			<< "*\t" << '\n';// QSTR
    }

};

class SAMwriter{
private:
	std::shared_ptr<spdlog::logger> samlog;
//	mapped_res mres;
	ref_pos rpos;
	std::vector<std::string> read_name;
public:
	SAMwriter(){
		samlog = spdlog::basic_logger_mt("basic_logger",SAM_FILE_LOC);
		samlog->set_pattern("..");
	}

	void read_mapped_info(){
		{
			std::ifstream read_name_file(PAIR_1_NAME_FILE);
			cereal::BinaryInputArchive ar_read_name(read_name_file);
			ar_read_name(read_name);
		}

		{
			std::ifstream ref_pos_file(MERGE_REF_POS_FILE);
			cereal::BinaryInputArchive ar_ref_pos(ref_pos_file);
			ar_ref_pos(rpos.ref_start,rpos.rname);
		}
	/*	{
			std::ifstream locfile(PAIR_1_LOC_FILE);
			cereal::BinaryInputArchive ar_loc(locfile);
			ar_loc(mres.mapped_ref_loc);
		}
		{
			std::ifstream disfile(PAIR_1_DIS_FILE);
			cereal::BinaryInputArchive ar_dis(disfile);
			ar_dis(mres.min_dis);
		}*/
		
	}

// TODO		
/*	SAM_format gen_SAM_format(int read_idx,int map_loc){
		
		// find the exact mapping position according to the size of transcriptomes
		size_t ridx = 0;
		for (int i = 0; i < rpos.ref_start.size(); ++i){
			if (rpos.ref_start[i] > mres.mapped_ref_loc[read_idx][map_loc]){
				ridx = i - 1;
				break;
			}else if(rpos.ref_start[i] == mres.mapped_ref_loc[read_idx][map_loc]){
				ridx = i;
				break;
			}
		}

		SAM_format sf;
		sf.rname = std::move(rpos.rname[ridx]);
		sf.pos = std::move(mres.mapped_ref_loc[read_idx][map_loc] - rpos.ref_start[ridx]);
		sf.tlen = DIM;
		if (mres.min_dis[read_idx] == 0){
			sf.cigar = "=";
		}else{
			sf.cigar = "X";
		}
		
		sf.qname = read_name[read_idx];
		
		
		sf.rnext = "*";
		sf.seq = "*";
		sf.qual = "*";
		

		sf.flag = 0;
		
		sf.mapq = 0;
		sf.pnext = 0;
		sf.tlen = 0;

		return sf;
	}*/
// TODO
	void gen_SAM_file(mapped_res &mres){
	//	mres = std::move(res);
		read_mapped_info();
		SAM_format sf;
		std::ofstream samfile("res/res.sam");

		for (int read_idx = 0; read_idx < mres.mapped_ref_loc.size(); ++read_idx)
		{
			for (int loc = 0; loc < mres.mapped_ref_loc[read_idx].size(); ++loc)
			{
			//	sam = gen_SAM_format(read_idx,loc);
			//	samlog->info(SAM_format{sam.qname,sam.rname,sam.cigar,sam.rnext,sam.seq,sam.qual,sam.flag,sam.pos,sam.mapq,sam.pnext,sam.tlen});

				// find the exact mapping position according to the size of transcriptomes
				size_t ridx = 0;
				for (int i = 0; i < rpos.ref_start.size(); ++i){
					if (rpos.ref_start[i] > mres.mapped_ref_loc[read_idx][loc]){
						ridx = i - 1;
						break;
					}else if(rpos.ref_start[i] == mres.mapped_ref_loc[read_idx][loc]){
						ridx = i;
						break;
					}
				}
				sf.rname = std::move(rpos.rname[ridx]);
				sf.pos = std::move(mres.mapped_ref_loc[read_idx][loc] - rpos.ref_start[ridx]);
				sf.tlen = DIM;
				if (mres.min_dis[read_idx] == 0){
					sf.cigar = "=";
				}else{
					sf.cigar = "X";
				}
				
				sf.qname = read_name[read_idx];
			/*	
				sf.rnext = "*";
				sf.seq = "*";
				sf.qual = "*";
				

				sf.flag = 0;
				
				sf.mapq = 0;
				sf.pnext = 0;
				sf.tlen = 0;*/

				// TODO change to fmt os lib
				samfile << sf.qname << '\t' // QNAME
					<< '*' << '\t' // FLAGS
					<< sf.rname << '\t' // RNAME
					<< sf.pos + 1 << '\t' // POS (1-based)
					<< 255 << '\t' // MAPQ
					<< sf.cigar << '\t' // CIGAR
					<< '*' << '\t' // MATE NAME
					<< 0 << '\t' // MATE POS
					<< sf.tlen << '\t' // TLEN
					<< '*' << '\t' // SEQ
					<< "*\t" << '\n';// QSTR
			}
		}
		
	}
};

