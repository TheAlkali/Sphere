#include <iostream>
#include <fstream>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include "fastxParser.hpp"
#include "Utils.hpp"

/*#define REAL_TYPE int
#define REF_LOCATION std::std::vector<int>

*/

class SAMparser
{
public:
	singleSeqList *reflist ;
	int read_index;
	int dim;
	int nP;


	SAMparser(){
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