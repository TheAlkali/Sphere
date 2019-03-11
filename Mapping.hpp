#pragma once

#include <sstream>
#include <iostream>
#include <cassert>
#include <memory>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/bitset.hpp>

#include <divsufsort.h> 

#include "math.h"
#include "string.h"
#include "Points.h"
#include "BooPHF.h"
#include "Utils.hpp"
#include "fastxParser.hpp"
#include "SAMwriter.hpp"



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

class Mapping
{
private:
    // then tell BBhash to use this custom hash : (also appears below, line 104)
    typedef boomphf::mphf< std::string, Custom_string_Hasher > boophf_t;

    // BBhash perfect hash function
    boophf_t * bphf = NULL;

public:
    region_profile rpro;
    std::string ref_string;
    SArray sarry;
    Points read_buff;

    int read_size;
    int dim;

    Mapping(int size,int rdim)
    {
        read_size = size;
        dim = rdim; 
    }

    void Load_SA(int seg_len)
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

        // load SA region
        T0.Reset();     T0.Start();
        {
            std::ifstream region_file("sa/ref_suffix_region.bin");
            cereal::BinaryInputArchive ar(region_file);
            ar(rpro.region_start_idx,rpro.region_end_idx,rpro.rkmer_idx,rpro.region_visited);
        }
        rpro.region_loc_in_file.resize(rpro.region_start_idx.size());
        T0.Stop();
        printf("- Load Suffix Array Region Finished (%f seconds)\n",T0.GetTime());   
    //-----------------------------------
        std::string tmp = TRANSCRIPTS_FILE_NAME;
        ref_string = merge_ref_seq(strdup(tmp.c_str()),seg_len);
     
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
            //    std::cout << "region seq:" << region_seq.size() << std::endl;
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
    void Suffix_Array(int seg_len)
    {   
        system("rm sa/*");
        system("rm tmp/*");
        Stopwatch T0("");
        T0.Reset();     T0.Start();
        //------------ generate suffix array------------
        // get the whole transcripts from fasta file
        std::string tmp = TRANSCRIPTS_FILE_NAME;
        ref_string = merge_ref_seq(strdup(tmp.c_str()),seg_len);
     
        sarry.size = ref_string.size();
        std::cout << "- Total Ref Length:" << sarry.size << std::endl;

        // allocate sa_ptr
        sarry.con_ptr = (int *)malloc(sarry.size * sizeof(int));    

        // construct sa_ptr
        divsufsort((unsigned char *)ref_string.c_str(),sarry.con_ptr,sarry.size);

        T0.Stop();
        printf("- Constructing Suffix Array Finished (%f seconds)\n",T0.GetTime());

        // save SA to disk
        T0.Reset();     T0.Start();
        sarry.tran2vec();
        {    
            std::ofstream safile("sa/sa_idx.bin");
            cereal::BinaryOutputArchive ar(safile);
            ar(CEREAL_NVP(sarry.con));
        }
        T0.Stop();
        printf("- Save SA Finished (%f seconds)\n",T0.GetTime());

        // generate rpro according to the first kmer in reads. KMER_SIZE = 5
        std::vector<std::string> region_seq;
        bool start = false;
        std::string last_kmer = "";
        std::string kmer = "";
    //    std::ofstream saidx_txt;

        // save SA region
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

        std::cout << "- Start To Building SA Region Index ..." << std::endl;
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
        printf("- Building SA Region Index Using BBhash Finished (%f seconds)\n",T1.GetTime());

    //    std::ofstream region_file("sa/ref_suffix_region.txt");
        uint64_t idx;
        for (unsigned int i = 0; i < rpro.region_start_idx.size(); ++i)
        {
            idx = bphf->lookup(region_seq[i]);
            //region_file << rpro.region_start_idx[i] << "\t" << rpro.region_end_idx[i] << "\t" << region_seq[i] << "\t" << idx 
            //            << "\t" << "0" << std::endl;
            rpro.rkmer_idx.push_back(idx);
            rpro.region_visited.push_back(0);
        }
     //   region_file.close();

        //save sa region profile to disk using cereal
        {
            std::ofstream region_file("sa/ref_suffix_region.bin");
            cereal::BinaryOutputArchive ar(region_file);
            ar(rpro.region_start_idx,rpro.region_end_idx,rpro.rkmer_idx,rpro.region_visited);
        }
        T0.Stop();
        printf("- Save SA Regions According To The Kmers Finished (%f seconds)\n",T0.GetTime());
    }
/*
    void Generate_Trainging_Data()
    {
        Stopwatch T1("");
        T1.Reset();     T1.Start();
        std::ofstream training_i(TRAINING_DATA_I);
        std::ofstream training_j(TRAINING_DATA_J);
        std::ofstream similarity_osp(SIMILARITY_OSP);
        srand( (unsigned int)( time(NULL) ) );
        // generate similar seq
        for (int i = 0; i < NUM_TRAIN_SAMPLES / 2; ++i)
        {
            int rand_region = rand() % (rpro.region_start_idx.size());
            int region_size = rpro.region_end_idx[rand_region] - rpro.region_start_idx[rand_region];
            int rand_loc_1 = rpro.region_start_idx[rand_region];//rand() % (region_size) + rpro.region_start_idx[rand_region];
            int rand_loc_2 = rand_loc_1 + 1;

            int dis = 0;
            for (int j = 0; j < dim - KMER_SIZE; ++j)
            {
                if (j == dim - KMER_SIZE - 1)
                {
                    training_i << ref_string[sarry.con[rand_loc_1] + j];
                    training_j << ref_string[sarry.con[rand_loc_2] + j];
                }else
                {
                    training_i << ref_string[sarry.con[rand_loc_1] + j] << ',';
                    training_j << ref_string[sarry.con[rand_loc_2] + j] << ',';
                }
                dis += abs(ref_string[sarry.con[rand_loc_1] + j] - ref_string[sarry.con[rand_loc_2] + j]);
            }
            training_i << std::endl;
            training_j << std::endl;
            similarity_osp << dis / 10 << ',';    
        }

        // generate unsimilar seq
        for (int i = 0; i < NUM_TRAIN_SAMPLES / 2; ++i)
        {
            int rand_region_1 = rand() % (rpro.region_start_idx.size());
            int rand_region_2 = rand() % (rpro.region_start_idx.size());
            int rand_loc_1,rand_loc_2;
            int dis = 0;
            if (rand_region_2 != rand_region_1)
            {
                rand_loc_1 = rand() % (rpro.region_end_idx[rand_region_1] - rpro.region_start_idx[rand_region_1]) + rpro.region_start_idx[rand_region_1];
                rand_loc_2 = rand() % (rpro.region_end_idx[rand_region_2] - rpro.region_start_idx[rand_region_2]) + rpro.region_start_idx[rand_region_2];
            }

            for (int j = 0; j < dim - KMER_SIZE; ++j)
            {
               if (j == dim - 1)
                {
                    training_i << ref_string[sarry.con[rand_loc_1] + j];
                    training_j << ref_string[sarry.con[rand_loc_2] + j];
                }else
                {
                    training_i << ref_string[sarry.con[rand_loc_1] + j] << ',';
                    training_j << ref_string[sarry.con[rand_loc_2] + j] << ',';
                }
                dis += abs(ref_string[sarry.con[rand_loc_1] + j] - ref_string[sarry.con[rand_loc_2] + j]);
            }
            training_i << std::endl;
            training_j << std::endl;
            if (i == NUM_TRAIN_SAMPLES / 2 - 1)
            {
                similarity_osp << dis / 10;
            }else
            {
                similarity_osp << dis / 10 << ',';    
            }
        }
        training_j.close();
        training_i.close();
        similarity_osp.close();
        T1.Stop();
        std::cout << "- Generate Training Data Finished ("  << T1.GetTime() << " seconds)" << std::endl;
    }*/

    void REAL_TYPE_to_String(std::string &seq,REAL_TYPE *d)
    {
        seq = "";
        for(int i=0;i < KMER_SIZE;i++)
        {
            seq.push_back((char)(d[i] + 48));
        }
    }


    region_info Get_Read_Region(REAL_TYPE *read)
    {
        region_info rinfo;
        std::fstream region_file;

        std::string kmer;
        REAL_TYPE_to_String(kmer,read);

        uint64_t rkmer_idx = bphf->lookup(kmer);
        uint64_t i;
        for (i = 0; i < rpro.rkmer_idx.size(); ++i)
        {
            if(rpro.rkmer_idx[i] == rkmer_idx)
            {
                rinfo.region_size = rpro.region_end_idx[i] - rpro.region_start_idx[i]; 
                rinfo.region_start_idx = rpro.region_start_idx[i];
                rinfo.region_visited = rpro.region_visited[i];
                rinfo.region_loc_in_rpro = i;
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
    // TODO
    void Compute_Ref_Code(SphericalHashing &src_sh,int pair)
    {
        // save read region to disk
        std::ofstream read_region_file;
        std::ofstream bCodeRef_file;
        if (pair == PAIR_1)
        {
            system("rm bin/*");
            system("rm tmp/*");
            read_buff.Initialize (read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_1);

            read_region_file.open(PAIR_1_REGION_FILE,std::ios::binary|std::ios::app);
            if(!read_region_file.is_open())
            {
                perror(PAIR_1_REGION_FILE);
            }

            bCodeRef_file.open(PAIR_1_REF_HASH_FILE,std::ios::binary);
        }else if (pair == PAIR_2)
        {
            read_buff.Initialize (read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_2);

            read_region_file.open(PAIR_2_REGION_FILE,std::ios::binary|std::ios::app);
            if(!read_region_file.is_open())
            {
                perror(PAIR_2_REGION_FILE);
            }

            bCodeRef_file.open(PAIR_2_REF_HASH_FILE,std::ios::binary);
        }

        Stopwatch T0("");
        T0.Reset();     T0.Start();
        
        std::vector<bitset<BCODE_LEN>> bCodeRef_vec;
        std::vector<region_info> read_region;       
        region_info rinfo;
        size_t loc_in_file = 0;
        bitset<BCODE_LEN> bCodeRef;
        for (int i = 0; i < read_size; ++i)
        {
            rinfo = Get_Read_Region(read_buff.d[i]);
            if (rinfo.region_visited == 0)
            {
                rinfo.region_loc_in_file = loc_in_file;
                rpro.region_loc_in_file[rinfo.region_loc_in_rpro] = loc_in_file;
                
                size_t region_end_idx = rinfo.region_start_idx + rinfo.region_size;
                for (size_t i = rinfo.region_start_idx; i < region_end_idx; ++i)
                {
                    REAL_TYPE *seq = new REAL_TYPE[dim];
                    for (int j = 0; j < dim; ++j)
                    {
                        seq[j] = ictoi_table[ref_string[sarry.con[i] + j]];
                    }
                    
                    src_sh.Compute_BCode(seq,bCodeRef);
                    bCodeRef_vec.push_back(bCodeRef);            
                }
                loc_in_file += rinfo.region_size;
            }else
            {
                rinfo.region_loc_in_file = rpro.region_loc_in_file[rinfo.region_loc_in_rpro];
            }    
            read_region.push_back(rinfo);    
        }
      // save read region info 
        {
             cereal::BinaryOutputArchive ar(read_region_file);
             ar(read_region);
        }
        {
             cereal::BinaryOutputArchive ar(bCodeRef_file);
             ar(bCodeRef_vec);
        }
       
    }

    int Hash_Mapping_with_SA(SphericalHashing &src_sh,int pair)
    {
        std::ifstream read_region_file;
        std::vector<region_info> read_region;

        std::ifstream bCodeRef_file;
        std::vector<bitset<BCODE_LEN>> bCodeRef_vec;
        if (pair == PAIR_1)
        {
            read_buff.Initialize(read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_1);

            read_region_file.open(PAIR_1_REGION_FILE);

            bCodeRef_file.open(PAIR_1_REF_HASH_FILE);
        }else if(pair == PAIR_2)
        {
            read_buff.ReleaseMem();
            read_buff.Initialize(read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_2);

            read_region_file.open(PAIR_2_REGION_FILE);

            bCodeRef_file.open(PAIR_2_REF_HASH_FILE);
        }

        {
            cereal::BinaryInputArchive ar_read(read_region_file);
            ar_read(read_region);

            cereal::BinaryInputArchive ar_ref(bCodeRef_file);
            ar_ref(bCodeRef_vec);
        }

        Stopwatch T0("");
        T0.Reset();     T0.Start();

        int ref_size = 0;

        std::vector<size_t> ref_loc;
        mapped_res mres;

        region_info rinfo;
        bitset<BCODE_LEN> bCodeRead;
        std::vector<bitset<BCODE_LEN>> bCodeRead_vec;
        size_t loc_in_file = 0;
        for (int i = 0; i < read_size; ++i)
        {
            src_sh.Compute_BCode<REAL_TYPE>(read_buff.d[i], bCodeRead);
            bCodeRead_vec.push_back(bCodeRead);

         //   Mapping_Process
            
            ref_loc.clear();
            if(read_region[i].region_size != 0)
            {
                ref_loc.clear();
                int min_dis = 1000;

                for (unsigned int j = 0; j < read_region[i].region_size; ++j)
                {    
                    loc_in_file = read_region[i].region_loc_in_file + j;
                   
                    int dist = Compute_HD(bCodeRead, bCodeRef_vec[loc_in_file]);
                    if (min_dis > dist)
                    {
                        min_dis = std::move(dist);
                        ref_loc.clear();    
                        ref_loc.push_back(sarry.con[read_region[i].region_start_idx + j]);
                    }else if (min_dis == dist)
                    {
                        ref_loc.push_back(sarry.con[read_region[i].region_start_idx + j]);
                    }
                } 
                mres.mapped_ref_loc.push_back(ref_loc);
                mres.min_dis.push_back(min_dis);
            }

            ref_size += read_region[i].region_size;
        }
    
        std::cout << "per ref size:" << ref_size << std::endl;
        T0.Stop();
        printf("- Total Mapping Time: (%f seconds)\n",T0.GetTime());   
         //save to disk using cereal
        {
            std::string loc_fname;
            std::string dis_fname;
            std::string read_code_fname;
            if (pair == 1)
            {
                loc_fname = PAIR_1_LOC_FILE;
                dis_fname = PAIR_1_DIS_FILE;
                read_code_fname = "tmp/read_code_1.bin";
            }else if(pair == 2)
            {
                loc_fname = PAIR_2_LOC_FILE;
                dis_fname = PAIR_2_DIS_FILE;
                read_code_fname = "tmp/read_code_2.bin";
            }
            std::ofstream tmp_loc_file(loc_fname);
            std::ofstream tmp_dis_file(dis_fname);
            std::ofstream read_code_file(read_code_fname);
            cereal::BinaryOutputArchive ar_loc(tmp_loc_file);
            cereal::BinaryOutputArchive ar_dis(tmp_dis_file);
            cereal::BinaryOutputArchive ar_read(read_code_file);
            
            ar_read(CEREAL_NVP(bCodeRead_vec));

            ar_loc(CEREAL_NVP(mres.mapped_ref_loc));
            ar_dis(CEREAL_NVP(mres.min_dis));
        }
       // close memory mapping
        bCodeRef_file.close();   
    }
};