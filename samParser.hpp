#include <iostream>
#include <fstream>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "fastxParser.hpp"
#include "Utils.h"

/*#define REAL_TYPE int
#define REF_LOCATION std::std::vector<int>

struct read_mapped_ref{
	std::vector<REF_LOCATION> locs_of_ref ;
	std::string read_name;
	std::vector<std::string> ref_name;
};

struct read_info{
	REAL_TYPE **read_seq;
	std::string *read_name;
	int nP;
	int dim;
};*/

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

class samParser
{
public:
	mapped_res mres;
	ref_pos rpos;
	singleSeqList *reflist ;
	int read_index;
	int dim;
	int nP;


	samParser(){
		reflist = new singleSeqList();
		read_index = 0;
		dim = 0;
		nP = 0;
	}

	void initialize(int _nP){
		nP = _nP;
		dim = DIM;
	/*	rinfo.nP = nP;
		rinfo.dim_ = dim;
		rinfo.read_seq = new REAL_TYPE * [ nP ];
		for(int i=0;i<nP;i++)
		{
			rinfo.read_seq[i] = new REAL_TYPE [ dim ];
		}
		rinfo.read_name = new std::string [nP];
		rmref = new read_mapped_ref[nP];*/
	}

	void read_mapped_info(){
		{
			std::ifstream refposfile(MERGE_REF_POS_FILE);
			cereal::BinaryInputArchive ar_ref_pos(refposfile);
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
	SAM_format gen_SAM_format(int read_idx,int map_loc){
		
		size_t ridx = 0;
		for (int i = 0; i < rpos.ref_start.size(); ++i){
			if (rpos.ref_start[i] > mres.mapped_ref_loc[read_idx][map_loc]){
				ridx = i - 1;
				break;
			}else if(rpos.ref_start[i] == mres.mapped_ref_loc[read_idx][map_loc]){
				ridx = i;
			}
		}

		SAM_format sf;
		sf.qname = "";
		sf.rname = rpos.rname[ridx];
		sf.cigar = "";
		sf.rnext = "";
		sf.seq = "";
		sf.qual = "";
		sf.tlen = DIM;

		sf.flag = 0;
		sf.pos = mres.mapped_ref_loc[read_idx][map_loc] - rpos.ref_start[ridx];
		sf.mapq = 0;
		sf.pnext = 0;
		sf.tlen = 0;

		return sf;
	}
// TODO
	void gen_SAM_file(mapped_res &res){
		mres = std::move(res);
		read_mapped_info();
		SAM_format sam;
		auto write_sam_res = spdlog::basic_logger_mt("basic_logger",SAM_FILE_LOC);
		for (int read_idx = 0; read_idx < mres.mapped_ref_loc.size(); ++read_idx)
		{
			for (int loc = 0; loc < mres.mapped_ref_loc[read_idx].size(); ++loc)
			{
				sam = gen_SAM_format(read_idx,loc);
				write_sam_res->info(SAM_format{sam.qname,sam.rname,sam.cigar,sam.rnext,sam.seq,sam.qual,sam.flag,sam.pos,sam.mapq,sam.pnext,sam.tlen});
			}
		}
		
	}
//-------------parse sam files-----------------------------
	int get_refseq(std::string ref_name){
		for (int i = 0; i < reflist->name.size(); ++i){
			if (reflist->name[i] == ref_name){
				return i;
			}
		}
		return -1;
	}

	std::string get_One_Of_Sam(std::ofstream &out_res,std::ifstream &sam,std::string last_read_name){
		std::string line;
		char tmp;
		if (!sam.is_open()){
			std::cerr << "can't open sam file." <<std::endl;
		}
		//skip head
		sam.get(tmp);
		while(tmp == '@'){
			std::getline(sam,line);
			sam.get(tmp);
		}
		sam.seekg(-1,std::ios_base::cur);

		std::getline(sam,line);
		std::stringstream ss(line);
		std::string read_name;
		std::string ref_name;
		int ref_loc;
		int sam_index = 0;
		std::string read_seq;

		ss >> read_name;
		sam_index++;
		std::cout << read_name << std::endl;
		for (int i = 0; i < 2; ++i){
			//ref name that one read_idx mapped to
			ss >> ref_name;
			sam_index++;
		}
		//the loc of ref that one read_idx mapped to
		ss >> ref_loc;
		//the index in sam start at 1
		ref_loc--;
		sam_index++;
		while(ss.peek() != EOF){
			ss >> read_seq;
			sam_index++;
			if (sam_index == 10){
				break;
			}
		}

		if (read_name != last_read_name){
			out_res << ">" << read_name << std::endl;
			if (sam_index < 12){
				if (sam_index == 10 && read_seq != "*"){
					out_res << "=" << read_seq << std::endl;		
				}else{
					std::getline(sam,line);
					ss << line;
					while(ss.peek() != EOF){
						ss >> read_seq;
						sam_index++;
						if (sam_index == 10 && read_seq != "*"){
							out_res << "=" << read_seq << std::endl;
							break;
						}
					}
				}
			}
		}	

		
		//print the ref seq that this read_idx mapped to
		out_res << "+";
		int rind = get_refseq(ref_name);
		for (int i = 0; i < dim; ++i){
			out_res << reflist->seq[rind][ref_loc + i];
		}
		out_res << std::endl;

		return read_name;

	/*	//add info to read_mapped_ref
		rmref[read_index].read_name = read_name;
		rmref[read_index].ref_name.push_back(ref_name);
		rmref[read_index].locs_of_ref.push_back(ref_loc);


		//add new read_idx to structure
		rinfo.read_name[read_index] = read_name;
		if (sam_index < 12){
			if (sam_index == 10 && read_seq != '*'){
				//TODO
				for (int i = 0; i < rinfo.dim; ++i){
					rinfo.read_seq[read_index][i] = (REAL_TYPE)( read_seq[k] - 48)
				}
				read_index++;

			}else{
				std::getline(sam,line);
				ss << line;
				while(ss.peek() != EOF){
					ss >> read_seq;
					sam_index++;
					if (sam_index == 10) && read_seq != '*'{
						for (int i = 0; i < rinfo.dim; ++i){
							rinfo.read_seq[read_index][i] = (REAL_TYPE)( read_seq[i] - 48)
						}
						read_index++;
						break;
					}
				}
			}
			
		}*/
		
	}

	void get_Sam(char *ref_file,char *sam_file,char *sam_res_file,int _nP,int _dim){
		//extract ref seq to RAM
		read_fastx(ref_file,reflist);

		std::ifstream sam;
		std::ofstream out_res;
		sam.open(sam_file);
		out_res.open(sam_res_file);
		if (!sam.is_open())
		{
			std::cerr << "can't open sam file" << std::endl;
		}

		std::string last_read_name;
		read_index = 0;
		nP = _nP;
		dim = _dim;
		while(sam.peek() != EOF){
			last_read_name = get_One_Of_Sam(out_res,sam,last_read_name);
		}
		sam.close();
		out_res.close();
		std::cout << "Parsing Sam Finished ." << std::endl;
	}

};