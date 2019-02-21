#pragma once

#include <sstream>
#include <iostream>
#include <cassert>
#include <memory>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

#include "string.h"
#include "Points.h"
#include "BooPHF.h"
#include "Utils.hpp"
#include "fastxParser.hpp"
#include "SAMwriter.hpp"


#include <divsufsort.h> 


std::string refseq_file = "dataset/srrdata/true_ref_seq.txt";

// custom hash for BBhash
class Custom_string_Hasher
{
public:
    // the class should have operator () with this signature :
    uint64_t operator ()   (std::string key, uint64_t seed=0) const
    {
        

        uint64_t hash  =  hash_fn(key);
        hash ^= seed;
        
        return hash;
    }
    
     std::hash<std::string> hash_fn;
};



// then tell BBhash to use this custom hash : (also appears below, line 104)
typedef boomphf::mphf< std::string, Custom_string_Hasher > boophf_t;

// BBhash perfect hash function
boophf_t * bphf = NULL;

region_profile rpro;

//int *sa_ptr ;
SArray sarry;
std::string ref_string;

void Load_SA()
{
    // load SA 
    Stopwatch T0("");
    T0.Reset();     T0.Start();
    {
        std::unique_ptr<int> sa_uniptr;
        std::ifstream saidx("sa/sa_idx.bin");
        cereal::BinaryInputArchive ar(saidx);
        ar(sarry.con);
    }
    sarry.size = sarry.con.size();  
    T0.Stop();
    printf("- Load Suffix Array Finished (%f seconds)\n",T0.GetTime());

    // ;pad SA region
    T0.Reset();     T0.Start();
    {
        std::ifstream region_file("sa/ref_suffix_region.bin");
        cereal::BinaryInputArchive ar(region_file);
        ar(rpro.region_start_idx,rpro.region_end_idx,rpro.rkmer_idx,rpro.region_visited);
    }
    T0.Stop();
    printf("- Load Suffix Array Region Finished (%f seconds)\n",T0.GetTime());   
//-----------------------------------
    std::string tmp = TRANSCRIPTS_FILE_NAME;
    ref_string = merge_ref_seq(strdup(tmp.c_str()));
 
    std::vector<std::string> region_seq;
    std::string last_kmer = "";
    std::string kmer = "";
    bool start = false;
    for (size_t i = 0; i < sarry.size; ++i) 
    {
//        saidx_txt << sa_ptr[i] << std::endl;
        kmer.clear();
/*        ref_seq.clear();
        ref_seq.seekg(sa_ptr[i],std::ios_base::beg);
        saidx_txt << sa_ptr[i] << std::endl;
        char base;
        for (int j = 0; j < KMER_SIZE; ++j)
        {
            ref_seq.get(base);
            kmer.push_back(base);
        }*/

       kmer = ref_string.substr(sarry.con[i],KMER_SIZE);

        if (last_kmer != kmer && !start)
        {
            region_seq.push_back(kmer);
            start = true;
        }else if (last_kmer != kmer && start)
        {
            region_seq.push_back(kmer);
        }
        last_kmer = kmer;
    }
//    saidx_txt.close();

    Stopwatch T1("");
    T1.Reset();     T1.Start();

    //-----compute the hashcode for ref_kmer using BBhash------
    // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
    // gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )
    double gammaFactor = 2.0; 
    int nthreads = 1;
    u_int64_t kmer_size = region_seq.size();
//    auto data_iterator = boomphf::range(static_cast<const REAL_TYPE*>(data), static_cast<const u_int64_t*>(data+nelem));

    //build the mphf
    bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(kmer_size,region_seq,nthreads,gammaFactor);
    //---------------------------------------------------------

    T1.Stop();
    printf("- Computing Hashcode Of Ref Kmer Using BBhash Finished (%f seconds)\n",T1.GetTime());
}

// change the library of constructing sa_ptr, then set idx_name to "tmp/all_ref_seq.bin"
void Suffix_Array()
{   
    system("rm sa/*");
    system("rm tmp/*");
    Stopwatch T0("");
    T0.Reset();     T0.Start();
    //------------ generate suffix array------------
    // get the whole transcripts from fasta file
    std::string tmp = TRANSCRIPTS_FILE_NAME;
    ref_string = merge_ref_seq(strdup(tmp.c_str()));
 
    sarry.size = ref_string.size();
    std::cout << "- Total Ref Length:" << sarry.size << std::endl;

    // allocate sa_ptr
    sarry.con_ptr = (int *)malloc(sarry.size * sizeof(int));    

    // construct sa_ptr
    divsufsort((unsigned char *)ref_string.c_str(),sarry.con_ptr,sarry.size);

    T0.Stop();
    printf("- Constructing Suffix Array Finished (%f seconds)\n",T0.GetTime());

    // save SA to disk
    sarry.tran2vec();
    {    
        std::ofstream safile("sa/sa_idx.bin");
        cereal::BinaryOutputArchive ar(safile);
        ar(CEREAL_NVP(sarry.con));
        std::cout << "- Save SA Finished" << std::endl;
    }

    // generate rpro according to the first kmer in reads. KMER_SIZE = 5
    std::vector<std::string> region_seq;
    bool start = false;
    std::string last_kmer = "";
    std::string kmer = "";
    std::ifstream ref_seq;
//    std::ofstream saidx_txt;
    ref_seq.open(refseq_file);

    T0.Reset();     T0.Start();
    for (size_t i = 0; i < sarry.size; ++i) 
    {
//        saidx_txt << sa_ptr[i] << std::endl;
        kmer.clear();
/*        ref_seq.clear();
        ref_seq.seekg(sa_ptr[i],std::ios_base::beg);
        saidx_txt << sa_ptr[i] << std::endl;
        char base;
        for (int j = 0; j < KMER_SIZE; ++j)
        {
            ref_seq.get(base);
            kmer.push_back(base);
        }*/

       kmer = ref_string.substr(sarry.con[i],KMER_SIZE);

        if (last_kmer != kmer && !start)
        {
            rpro.region_start_idx.push_back(i);
            region_seq.push_back(kmer);
            start = true;
        }else if (last_kmer != kmer && start)
        {
            rpro.region_end_idx.push_back(i);
            rpro.region_start_idx.push_back(i);
            region_seq.push_back(kmer);
        }
        last_kmer = kmer;
    }
//    saidx_txt.close();

    Stopwatch T1("");
    T1.Reset();     T1.Start();

    //-----compute the hashcode for ref_kmer using BBhash------
    // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
    // gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )
    double gammaFactor = 2.0; 
    int nthreads = 1;
    u_int64_t kmer_size = region_seq.size();
//    auto data_iterator = boomphf::range(static_cast<const REAL_TYPE*>(data), static_cast<const u_int64_t*>(data+nelem));

    //build the mphf
    bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(kmer_size,region_seq,nthreads,gammaFactor);
    //---------------------------------------------------------

    T1.Stop();
    printf("- Computing Hashcode Of Ref Kmer Using BBhash Finished (%f seconds)\n",T1.GetTime());

    std::ofstream region_file("sa/ref_suffix_region.txt");
    uint64_t idx;
    for (unsigned int i = 0; i < rpro.region_start_idx.size(); ++i)
    {
        idx = bphf->lookup(region_seq[i]);
        region_file << rpro.region_start_idx[i] << "\t" << rpro.region_end_idx[i] << "\t" << region_seq[i] << "\t" << idx 
                    << "\t" << "0" << std::endl;
        rpro.rkmer_idx.push_back(idx);
        rpro.region_visited.push_back(0);
    }
    region_file.close();

    //save sa region profile to disk using cereal
    {
        std::ofstream region_file("sa/ref_suffix_region.bin");
        cereal::BinaryOutputArchive ar(region_file);
        ar(rpro.region_start_idx,rpro.region_end_idx,rpro.rkmer_idx,rpro.region_visited);
    }
    T0.Stop();
    printf("- Save SA Regions According To The Kmers Finished (%f seconds)\n",T0.GetTime());
}

void REAL_TYPE_to_String(std::string &seq,REAL_TYPE *d)
{
    seq = "";
    for(int i=0;i < KMER_SIZE;i++)
    {
        seq.push_back((char)(d[i] + 48));
    }
}


/* TODO 
    1.how to search faster? perfect hashing?
        store region info in RAM,use bbhash to compute the kmer idx,and use it to look for the region
    2.how to get the rpro when the kmer of reads  mismatch?
*/
region_info Get_Read_Region(REAL_TYPE *read)
{
    region_info rinfo;
    std::fstream region_file;
//  region_file.open("sa_ptr/ref_suffix_region.txt",std::ios_base::in | std::ios_base::out);

    std::string kmer;
    REAL_TYPE_to_String(kmer,read);

    uint64_t rkmer_idx = bphf->lookup(kmer);
    uint64_t i;
    for (i = 0; i < rpro.rkmer_idx.size(); ++i)
    {
        if(rpro.rkmer_idx[i] == rkmer_idx)
        {
//            std::cout << i << "\t" << rpro.region_start_idx[i] << "\t" << rpro.region_end_idx[i] << std::endl;

            auto size = rpro.region_end_idx[i] - rpro.region_start_idx[i]; 
            rinfo.region_start_idx = rpro.region_start_idx[i];
            rinfo.region_visited = rpro.region_visited[i];
            for (size_t j = rpro.region_start_idx[i]; j < rpro.region_end_idx[i]; ++j)
            {
                rinfo.idx_set.push_back(sarry.con[j]);   
            }

            // change the region state to visited
            if (!rpro.region_visited[i])
            {
                rpro.region_visited[i] = 1;
            }
            
            break;
        }
    }
    return rinfo;
}

void Mapping_Process(mapped_res &mres,region_info &rinfo,bitset<BCODE_LEN> bCodeRead,SphericalHashing &src_sh)
{
    // mapping location in ref
    std::vector<size_t> ref_loc;
    bitset<BCODE_LEN> bCodeRef;
    if(rinfo.idx_set.size() != 0)
    {
        ref_loc.clear();
        int min_dis = 1000;

        std::string fname = "bin/";
        fname.append(std::to_string(rinfo.region_start_idx));
        fname.append(".bin");

        // if region hasn't been visited, computing their hashcode
        if (!rinfo.region_visited)
        {                
            std::ofstream recode_file; 
            recode_file.open(fname);
            if(!recode_file.is_open())
            {
                perror(fname.c_str());
            }
            for (unsigned int i = 0; i < rinfo.idx_set.size(); ++i)
            {
                REAL_TYPE *seq = new REAL_TYPE[DIM];
                for (int j = 0; j < DIM; ++j)
                {
                    seq[j] = (REAL_TYPE)(ref_string[rinfo.idx_set[i] + j] - 48);
                   // std::cout << seq[j];
                }
                //std::cout << std::endl;
                src_sh.Compute_BCode(seq,bCodeRef);
                recode_file << bCodeRef << std::endl;         
                
                // compute distance
                int dist = Compute_HD(bCodeRead,bCodeRef);
                if (min_dis > dist)
                {
                    min_dis = std::move(dist);
                    ref_loc.clear();
                    ref_loc.push_back(rinfo.idx_set[i]);
                }else if (min_dis == dist)
                {
                    ref_loc.push_back(rinfo.idx_set[i]);
                }
            
            }
            recode_file.close();
        }else if(rinfo.region_visited)
        {
            std::string code;
            std::ifstream recode_file;
            recode_file.open(fname);
            if(!recode_file.is_open())
            {
                perror(fname.c_str());
            }
            for (unsigned int i = 0; i < rinfo.idx_set.size(); ++i)
            {            
                getline(recode_file,code);
                bCodeRef = String2Bit<BCODE_LEN>(code);
                int dist = Compute_HD(bCodeRead,bCodeRef);
                if (min_dis > dist)
                {
                    min_dis = std::move(dist);
                    ref_loc.clear();
                    ref_loc.push_back(rinfo.idx_set[i]);
                }else if (min_dis == dist)
                {
                    ref_loc.push_back(rinfo.idx_set[i]);
                }
            } 
            recode_file.close();
        }    
        mres.mapped_ref_loc.push_back(ref_loc);
        mres.min_dis.push_back(min_dis);
    }
}

int Hash_Mapping_with_SA(SphericalHashing &src_sh,Points &read_buff,int pair)
{
    if (pair == 1)
    {
        system("rm bin/*");
        system("rm tmp/*");
    }
    
    int content = -1;
    std::string read;
/*    {
        cereal::BinaryInputArchive ar(ref_seq);
        ar(ref_string);
    }*/

    Stopwatch T0("");
    T0.Reset();     T0.Start();
    read_buff.srcfile.clear();
    read_buff.Initialize(READ_BUFFER_SIZE,DIM);

    std::vector<size_t> ref_loc;
    mapped_res mres;
    while(content < 0)
    {
        Stopwatch T1("");
        T1.Reset();     T1.Start();
        content = read_buff.Initialize_From_File(); 
        T1.Stop();
        int size;
        if (content == -1)
        {
            size = READ_BUFFER_SIZE;
        }else 
        {
            size = content;
        }
        printf("- Loading Reads Data In Buffer Finished (%f seconds)\n",T1.GetTime());
        
        T1.Reset();     T1.Start();
        
        region_info rinfo;
        for (int i = 0; i < size; ++i)
        {
            bitset<BCODE_LEN> bCodeRead;
            bitset<BCODE_LEN> bCodeRef;
           
            src_sh.Compute_BCode<REAL_TYPE>(read_buff.d[i], bCodeRead);

            rinfo.idx_set.clear();
            //std::cout << "---" << i << std::endl;
            rinfo = Get_Read_Region(read_buff.d[i]);
            Mapping_Process(mres,rinfo,bCodeRead,src_sh);
        }
        T1.Stop();
        printf("- Mapping Time: Mapping Reads Through Hashcodes and SA Finished (%f seconds)\n",T1.GetTime());
    }   
    Stopwatch T1("");
    //save to disk using cereal
    {
        std::string loc_fname;
        std::string dis_fname;
        if (pair == 1)
        {
            loc_fname = PAIR_1_LOC_FILE;
            dis_fname = PAIR_1_DIS_FILE;
        }else if(pair == 2)
        {
            loc_fname = PAIR_2_LOC_FILE;
            dis_fname = PAIR_2_DIS_FILE;
        }

        std::ofstream tmp_loc_file(loc_fname,std::ios::binary|std::ios::app);
        std::ofstream tmp_dis_file(dis_fname,std::ios::binary|std::ios::app);
        cereal::BinaryOutputArchive ar_loc(tmp_loc_file);
        cereal::BinaryOutputArchive ar_dis(tmp_dis_file);

        ar_loc(CEREAL_NVP(mres.mapped_ref_loc));
        ar_dis(CEREAL_NVP(mres.min_dis));
    }
    T0.Stop();
    printf("- Total Mapping Time: (%f seconds)\n",T0.GetTime());
 //   mres.mapped_ref_loc.~vector();
 //   mres.min_dis.~vector();
}