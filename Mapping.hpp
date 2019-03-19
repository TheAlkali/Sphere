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

    const int buffer_size = 10000;

    std::vector<std::string> region_kmer;

public:
    SphericalHashing src_sh;
    region_profile rpro;
    std::string ref_string;
    SArray sarry;

    int read_size;
    int dim;

    Mapping(int size,int rdim)
    {
        read_size = size;
        dim = rdim; 
    }

    void Learn_Spherical_Hashing(Points &buff,int code_len,int seg_len)
    {
        Stopwatch T0("");
        src_sh.Initialize(&buff,code_len,seg_len);
        T0.Reset();     T0.Start();
        src_sh.Set_Spheres();

        // TODO  save sh to file

        /*for(int i=0;i<BCODE_LEN;i++)
        {
            for (int j = 0; j < seg_len; ++j)
            {
                std::cout << src_sh.s[i].c[j] << "\t";
            }
            std::cout << std::endl;
            std::cout << src_sh.s[i].rSq << std::endl;
        }*/
        T0.Stop();
        printf("- Learning Spherical Hashing Finished (%f seconds)\n",T0.GetTime());
    //  sh.Save_Sphere_Info();
    }
    // TODO
    void Load_SH_Para()
    {

    }
    void Compute_All_Ref_Code()
    {
        system("rm bin/*");
        // save read region to disk
        std::ofstream bCodeRef_file;

        bCodeRef_file.open(REF_HASH_FILE,std::ios::binary);

        Stopwatch T0("");
        T0.Reset();     T0.Start();
        
        std::vector<unsigned long> bCodeRef_vec(buffer_size);

        size_t loc_in_file = 0;
        int ref_buffer_count = 0;
        bitset<BCODE_LEN> bCodeRef;
        REAL_TYPE *seq = new REAL_TYPE[dim];
        for (int i = 0; i < rpro.region_start_idx.size(); ++i)
        {
            if (region_kmer[i].find_first_of("9") == std::string::npos)
            {
                
                rpro.region_loc_in_file[i] = loc_in_file;
                
                for (size_t j = rpro.region_start_idx[i]; j < rpro.region_end_idx[i]; ++j)
                {
                    //if buffer is full, store ref codes to disk
                    if (ref_buffer_count == buffer_size)
                    {
                        bCodeRef_file.write(reinterpret_cast<const char*>(&bCodeRef_vec[0]),buffer_size * sizeof(unsigned long));
                        ref_buffer_count = 0;
                    }
                    
                    for (int k = 0; k < dim; ++k)
                    {
                        seq[k] = ictoi_table[ref_string[sarry.con[j] + k]];
                    }
                    src_sh.Compute_BCode(seq,bCodeRef);
                    bCodeRef_vec[ref_buffer_count] = bCodeRef.to_ulong();         
                    ref_buffer_count++;
                   
                }
                loc_in_file += rpro.region_end_idx[i] - rpro.region_start_idx[i];
            }
 
        }
        bCodeRef_file.close();    
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
            ar(rpro.rkmer_idx,rpro.region_start_idx,rpro.region_end_idx,rpro.region_loc_in_file);
        }

        {
            std::ifstream seq_file("sa/ref_suffix_region_kmer.bin");
            cereal::BinaryInputArchive ar(seq_file);
            ar(region_kmer);
        }
        T0.Stop();
        printf("- Load Suffix Array Region Finished (%f seconds)\n",T0.GetTime());   
    //-----------------------------------
        std::string tmp = TRANSCRIPTS_FILE_NAME;
        ref_string = merge_ref_seq(strdup(tmp.c_str()),seg_len); 

        T0.Reset();     T0.Start();

        //-----compute the hashcode for ref_kmer using BBhash------
        // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
        // gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )
        double gammaFactor = 2.0; 
        int nthreads = 1;
        u_int64_t kmer_size = region_kmer.size();
    //    auto data_iterator = boomphf::range(static_cast<const REAL_TYPE*>(data), static_cast<const u_int64_t*>(data+nelem));

        //build the mphf
        bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(kmer_size,region_kmer,nthreads,gammaFactor);
        //---------------------------------------------------------

        T0.Stop();
        printf("- Computing Hashcode Of Ref Kmer Using BBhash Finished (%f seconds)\n",T0.GetTime());

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

        // generate rpro according to the first kmer in reads. KMER_SIZE = 10
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
                region_kmer.push_back(kmer);
                start = true;
            }else if (last_kmer != kmer && start)
            {
                rpro.region_end_idx.push_back(i);
                rpro.region_start_idx.push_back(i);
                region_kmer.push_back(kmer);
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
        u_int64_t kmer_size = region_kmer.size();
    //    auto data_iterator = boomphf::range(static_cast<const REAL_TYPE*>(data), static_cast<const u_int64_t*>(data+nelem));

        //build the mphf
        bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(kmer_size,region_kmer,nthreads,gammaFactor);
        //---------------------------------------------------------

        T1.Stop();
        printf("- Building SA Region Index Using BBhash Finished (%f seconds)\n",T1.GetTime());

    //    std::ofstream region_file("sa/ref_suffix_region.txt");
        uint64_t idx;
        uint64_t max_idx = 0;
        for (unsigned int i = 0; i < rpro.region_start_idx.size(); ++i)
        {
            idx = bphf->lookup(region_kmer[i]);
            if (idx > max_idx)
            {
                max_idx = idx;
            }
        //    region_file << rpro.region_start_idx[i] << "\t" << rpro.region_end_idx[i] << "\t" << region_kmer[i] << "\t" << idx << std::endl;
            rpro.rkmer_idx.push_back(idx);
        }
        std::cout << "max idx:" << max_idx << std::endl;
        std::vector<size_t> region_bphf_idx(max_idx + 1);
        for (int i = 0; i < rpro.region_start_idx.size(); ++i)
        {
            idx = bphf->lookup(region_kmer[i]);
            region_bphf_idx[idx] = i;
        }

        rpro.region_loc_in_file.resize(rpro.region_start_idx.size());
    //    region_file.close();

        T0.Reset();     T0.Start();
        printf("- Starting To Compute Ref Hashcode ...\n");
        Compute_All_Ref_Code();
        T0.Stop();
        printf("- Computing Hashcode Of Ref Using SH (%f seconds)\n",T0.GetTime());

        //save sa region profile to disk using cereal
        {
            std::ofstream region_file("sa/ref_suffix_region.bin");
            cereal::BinaryOutputArchive ar(region_file);
            ar(rpro.rkmer_idx,rpro.region_start_idx,rpro.region_end_idx,rpro.region_loc_in_file);
        }

        {
            std::ofstream seq_file("sa/ref_suffix_region_kmer.bin");
            cereal::BinaryOutputArchive ar(seq_file);
            ar(region_kmer);
        }

        {
            std::ofstream region_bphf_idx_file("sa/region_bphf_idx.bin");
            cereal::BinaryOutputArchive ar(region_bphf_idx_file);
            ar(region_bphf_idx);
        }
        region_kmer.clear();

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
            seq.push_back(itoic_table[(int)d[i]]);
        }
    }
  
    // TODO
    void Get_Read_Region(int pair)
    {
        // save read region to disk
        std::ofstream read_region_file;
        Points read_buff;

        if (pair == PAIR_1)
        {
            system("rm tmp/*");
            read_buff.Initialize (read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_1,Points::point_type::mapping);

            read_region_file.open(PAIR_1_REGION_FILE,std::ios::binary | std::ios::app);
            if(!read_region_file.is_open())
            {
                perror(PAIR_1_REGION_FILE);
            }

        }else if (pair == PAIR_2)
        {
            read_buff.Initialize (read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_2,Points::point_type::mapping);

            read_region_file.open(PAIR_2_REGION_FILE,std::ios::binary | std::ios::app);
            if(!read_region_file.is_open())
            {
                perror(PAIR_2_REGION_FILE);
            }
        }

        std::vector<size_t> region_bphf_idx;
        {
            std::ifstream region_bphf_idx_file("sa/region_bphf_idx.bin");
            cereal::BinaryInputArchive ar(region_bphf_idx_file);
            ar(region_bphf_idx);
        }

        Stopwatch T0("");
        T0.Reset();     T0.Start();
      
        std::vector<uint64_t> read_region(read_size);       
        region_info rinfo;
        std::string kmer;
        uint64_t rkmer_idx;
        uint64_t loc_in_rpro = 0;
        std::fstream region_file;
        for (int i = 0; i < read_size; ++i)
        {     
            REAL_TYPE_to_String(kmer,read_buff.d[i]);

            rkmer_idx = bphf->lookup(kmer);
            read_region[i] = region_bphf_idx[rkmer_idx];
        /*    for (loc_in_rpro = 0; loc_in_rpro < rpro.rkmer_idx.size(); ++loc_in_rpro)
            {
                if(rpro.rkmer_idx[loc_in_rpro] == rkmer_idx)
                {
                    read_region[i] = loc_in_rpro;  
                    break;
                }
            }*/   
        }
        // save read region info 
        {
            cereal::BinaryOutputArchive ar(read_region_file);
            ar(read_region);
        }

        read_region.clear();
        read_region_file.close();
    }

    int Hash_Mapping_with_SA(int pair)
    {
        std::ifstream read_region_file;
        std::vector<uint64_t> read_region;

        
        Points read_buff;

        std::ofstream tmp_loc_file;
        std::ofstream tmp_dis_file;
        std::ofstream read_code_file;

        // load ref code using memory mapping
        //MemoryMapped bCodeRef_file(REF_HASH_FILE);
        std::ifstream bCodeRef_file(REF_HASH_FILE);

        if (pair == PAIR_1)
        {
            read_buff.Initialize(read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_1,Points::point_type::mapping);

            read_region_file.open(PAIR_1_REGION_FILE);

        }else if(pair == PAIR_2)
        {
            read_buff.Initialize(read_size,dim);
            read_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_2,Points::point_type::mapping);

            read_region_file.open(PAIR_2_REGION_FILE);

        }

        {
            cereal::BinaryInputArchive ar_read(read_region_file);
            ar_read(read_region);
        }

        Stopwatch T0("");
        T0.Reset();     T0.Start();

        int ref_size = 0;

        std::vector<size_t> ref_loc;
        mapped_res mres;

        region_info rinfo;
        bitset<BCODE_LEN> bCodeRead;
        bitset<BCODE_LEN> bCodeRef;
    //    std::vector<bitset<BCODE_LEN>> bCodeRead_vec;
        size_t loc_in_file = 0;
        int region_size = 0;
        int min_dis = 0;
        int dist = 0;
        size_t start = 0;
        std::vector<unsigned long> bCodeRef_vec;
        for (int i = 0; i < read_size; ++i)
        {
            src_sh.Compute_BCode<REAL_TYPE>(read_buff.d[i], bCodeRead);
        //    bCodeRead_vec.push_back(bCodeRead);

         //   Mapping_Process
            
            ref_loc.clear();
            region_size = rpro.region_end_idx[read_region[i]] - rpro.region_start_idx[read_region[i]];
            bCodeRef_vec.resize(region_size);
            if(region_size > 0)
            {
                bCodeRef_file.clear();
                bCodeRef_file.seekg(rpro.region_loc_in_file[read_region[i]] * sizeof(unsigned long),bCodeRef_file.beg);
                bCodeRef_file.read(reinterpret_cast<char*>(&bCodeRef_vec[0]),region_size * sizeof(unsigned long));

            //    start = rpro.region_loc_in_file[read_region[i]] * sizeof(unsigned long);
                ref_loc.clear();
                min_dis = 1000;
             //   std::cout << region_size << std::endl;
                for (unsigned int j = 0; j < region_size; ++j)
                {    
                    bCodeRef = bCodeRef_vec[j]; //bCodeRef_file.at(start + j * sizeof(unsigned long));
                    dist = Compute_HD(bCodeRead, bCodeRef);
                    if (min_dis > dist)
                    {
                        min_dis = std::move(dist);
                        ref_loc.clear();    
                        ref_loc.push_back(sarry.con[rpro.region_start_idx[read_region[i]] + j]);
                    }else if (min_dis == dist)
                    {
                        ref_loc.push_back(sarry.con[rpro.region_start_idx[read_region[i]] + j]);
                    }
                    /*else if (min_dis == 0 && min_dis < dist)
                    {
                        break;
                    }*/
                } 
                mres.mapped_ref_loc.push_back(ref_loc);
                mres.min_dis.push_back(min_dis);
            }
            ref_size += region_size;

        /*    read_buffer_size++;
            if (read_buffer_size == buffer_size)
            {
                read_buffer_size = 0;
                read_region.clear();
                read_region_file.read(reinterpret_cast<char *>(&read_region[0]),buffer_size * sizeof(region_info));
                tmp_loc_file.write(reinterpret_cast<char *>(&mres.mapped_ref_loc[0]),buffer_size * sizeof(size_t));
                tmp_dis_file.write(reinterpret_cast<char *>(&mres.min_dis[0]),buffer_size * sizeof(int8_t));
            }*/
        }
    
        std::cout << "per ref size:" << ref_size << std::endl;
        T0.Stop();
        printf("- Total Mapping Time: (%f seconds)\n",T0.GetTime());   
        
        if (pair == PAIR_1)
        {
            tmp_loc_file.open(PAIR_1_LOC_FILE);
            tmp_dis_file.open(PAIR_1_DIS_FILE);

        //    read_code_file.open("tmp/read_code_1.bin");
        }else if(pair == PAIR_2)
        {
            tmp_loc_file.open(PAIR_2_LOC_FILE);
            tmp_dis_file.open(PAIR_2_DIS_FILE);

        //    read_code_file.open("tmp/read_code_2.bin");
        }
         //save to disk using cereal
        {
        //    cereal::BinaryOutputArchive ar_read(read_code_file);
            cereal::BinaryOutputArchive ar_loc(tmp_loc_file);
            cereal::BinaryOutputArchive ar_dis(tmp_dis_file);
            
        //    ar_read(CEREAL_NVP(bCodeRead_vec));

            ar_loc(CEREAL_NVP(mres.mapped_ref_loc));
            ar_dis(CEREAL_NVP(mres.min_dis));
        }

    //    bCodeRead_vec.clear(); 
        mres.mapped_ref_loc.clear(); 
        mres.min_dis.clear();

        bCodeRef_file.close();

    }
};