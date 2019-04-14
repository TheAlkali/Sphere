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
    typedef boomphf::mphf< std::string, Custom_string_Hasher > sa_boophf_t;
    // BBhash perfect hash function
    sa_boophf_t * sa_bphf = NULL;

    const int buffer_size = 10000;

    std::vector<std::string> region_kmer;

    std::vector<std::vector<int>> ref_of_read_1;
    std::vector<std::vector<int>> ref_of_read_2;

    std::vector<std::string> rcode_of_read_1;
    std::vector<std::string> rcode_of_read_2;

    std::vector<std::string> ref_name;
    std::vector<int> loc_to_ref;

    std::string ref_string;
    SArray sarry;
    SphericalHashing src_sh;

public:
    //-------- preserve after mapping------
    Points read_1_buff;
    Points read_2_buff;

    std::vector<bool> is_read_1_rev;
    std::vector<bool> is_read_2_rev;
    region_profile rpro;

    int read_size;
    int dim;

    Mapping(int rdim)
    {
        dim = rdim; 
    }
    ~Mapping(){};

/*    void Release_Content()
    {
        region_kmer.~vector();
        ref_of_read_1.~vector();
        ref_of_read_2.~vector();
        rcode_of_read_1.~vector();
        rcode_of_read_2.~vector();
        ref_name.~vector();
        loc_to_ref.~vector();
        ref_string.~string();
        sarry.~SArray();
        src_sh.~SphericalHashing();
        std::cout << "- Release Mapping Content Finished" << std::endl;
    }*/

// --------friend func of SAMwriter--------------
    std::string Get_Read_1_Buff(int idx)
    {
        std::string read;
        read.resize(dim);
        for (int i = 0; i < dim; ++i)
        {
            read[i] = itos_table[(int8_t)read_1_buff.d[idx][i] ];
        }
        return read;
    }

    std::string Get_Read_2_Buff(int idx)
    {
        std::string read;
        read.resize(dim);
        for (int i = 0; i < dim; ++i)
        {
            read[i] = itos_table[(int8_t)read_2_buff.d[idx][i]];
        }
        return read;
    }

    bool Get_Is_Read_1_Rev(int idx)
    {
        return is_read_1_rev[idx];
    }

    bool Get_Is_Read_2_Rev(int idx)
    {
        return is_read_2_rev[idx];
    }
//------------------------------------------------

    void Learn_Spherical_Hashing(Points &buff,int code_len,int seg_len)
    {
        Stopwatch T0("");
        src_sh.Initialize(&buff,code_len,seg_len);
        T0.Reset();     T0.Start();
        src_sh.Set_Spheres();

        T0.Stop();
        src_sh.Save_Sphere_Info();
        printf("- Learning Spherical Hashing Finished (%f seconds)\n",T0.GetTime());
    }

    void Load_Spherical_Hashing(int rsize,int code_len,int seg_len)
    {
        read_size = rsize;
        Stopwatch T0("");
        T0.Reset();     T0.Start();
        src_sh.Load_Sphere_Info(code_len,seg_len);
        T0.Stop();
        printf("- Loading Spherical Hashing Finished (%f seconds)\n",T0.GetTime());
    }

    void Compute_All_Ref_Code()
    {
    //    system("rm bin/*");
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
                src_sh.Compute_BCode<REAL_TYPE>(seq,bCodeRef);
            //    std::cout << bCodeRef << std::endl;
                bCodeRef_vec[ref_buffer_count] = bCodeRef.to_ulong();         
                ref_buffer_count++;
               
            }
            loc_in_file += rpro.region_end_idx[i] - rpro.region_start_idx[i];
 
        }
        bCodeRef_file.close();    
    }

    void Merge_Duplicated_In_Region()
    {
        std::ofstream irinfo_file(REDUCED_REGION_INFO_FILE);
        std::ofstream bucket_file(REDUCED_REGION_CODE_BUCKET_FILE);
        std::ifstream bCodeRef_file(REF_HASH_FILE);

        std::vector<unsigned long> bCodeRef_vec;
        std::vector<unsigned long> reduced_region_code;

        std::vector<std::vector<size_t>> reduced_code_bucket;
        std::vector<size_t> each_bucket_idx;

        // record the first location that the unique hashcode appear in sa and store them in bucket
        std::vector<int> loc_in_bucket;
        int region_size = 0;
        int max_element = 0;
        int unique_size = 0;
        int bucket_idx = 0;

        std::vector<size_t> loc_tmp;
        size_t region_start = 0;


        for (int i = 0; i < rpro.region_start_idx.size() - 1; ++i)
        {
            region_size = rpro.region_end_idx[i] - rpro.region_start_idx[i];
        //    std::cout <<rpro.region_end_idx[i] << "\t"<< rpro.region_start_idx[i]<< "\t"<< region_size << std::endl;
            bCodeRef_vec.resize(region_size);
            bCodeRef_file.clear();
            bCodeRef_file.seekg(rpro.region_start_idx[i] * sizeof(unsigned long),bCodeRef_file.beg);
            bCodeRef_file.read(reinterpret_cast<char*>(&bCodeRef_vec[0]),region_size * sizeof(unsigned long));

            //find the max element in vec
            max_element = *(std::max_element(std::begin(bCodeRef_vec), std::end(bCodeRef_vec)));
            loc_in_bucket.clear();
            loc_in_bucket.resize(max_element + 1,-1);

            unique_size = 0;
            reduced_region_code.clear();
            reduced_code_bucket.clear();
            for (int j = 0; j < region_size; ++j)
            {
                if (loc_in_bucket[bCodeRef_vec[j]] < 0)
                {
                    loc_in_bucket[bCodeRef_vec[j]] = unique_size;
                    reduced_region_code.push_back(bCodeRef_vec[j]);
                    unique_size++;
                    reduced_code_bucket.resize(unique_size);
                }
                reduced_code_bucket[loc_in_bucket[bCodeRef_vec[j]]].push_back(sarry.con[rpro.region_start_idx[i] + j]);
            }
            rpro.region_start_idx[i] = region_start;
            region_start += unique_size;
            rpro.region_end_idx[i] = region_start;

            // save unique ref code to file
            irinfo_file.write(reinterpret_cast<const char*>(&reduced_region_code[0]),unique_size * sizeof(unsigned long));
        
            for (int j = 0; j < unique_size; ++j)
            {
                // save the bucket of unique ref code to bucket file
                bucket_file.write(reinterpret_cast<char*>(&reduced_code_bucket[loc_in_bucket[reduced_region_code[j]]][0]),
                    reduced_code_bucket[loc_in_bucket[reduced_region_code[j]]].size() * sizeof(size_t));

                // save the start of the bucket region of each unique ref code 
                each_bucket_idx.push_back(bucket_idx);
                bucket_idx += reduced_code_bucket[loc_in_bucket[reduced_region_code[j]]].size();
            }
            rpro.code_bucket_idx.push_back(each_bucket_idx);
            each_bucket_idx.clear();
        }
        irinfo_file.close();
        bucket_file.close();
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
            ar(rpro.rkmer_idx,rpro.region_start_idx,rpro.region_end_idx,rpro.code_bucket_idx);
        }

        {
            std::ifstream seq_file("sa/ref_suffix_region_kmer.bin");
            cereal::BinaryInputArchive ar(seq_file);
            ar(region_kmer);
        }

        {
            std::ifstream ref_string_file(MERGE_REF_SEQ_FILE);
            cereal::BinaryInputArchive ar(ref_string_file);
            ar(ref_string);
        }
        T0.Stop();
        printf("- Load Suffix Array Region Finished (%f seconds)\n",T0.GetTime());   
    //-----------------------------------
        std::string tmp = TRANSCRIPTS_FILE_NAME;
        //ref_string = merge_ref_seq(strdup(tmp.c_str()),seg_len); 

        T0.Reset();     T0.Start();

        //-----compute the hashcode for ref_kmer using BBhash------
        // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
        // gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )
        double gammaFactor = 2.0; 
        int nthreads = 1;
        u_int64_t kmer_size = region_kmer.size();
    //    auto data_iterator = boomphf::range(static_cast<const REAL_TYPE*>(data), static_cast<const u_int64_t*>(data+nelem));

        //build the mphf
        sa_bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(kmer_size,region_kmer,nthreads,gammaFactor);
        //---------------------------------------------------------
        region_kmer.clear();
        T0.Stop();
        printf("- Computing Hashcode Of Ref Kmer Using BBhash Finished (%f seconds)\n",T0.GetTime());


        //analyse the region size distribution when k=10
    /*    int max_size = 0;
        int tmp_size = 0;
        std::vector<int> bucket(82653,0);
        std::ofstream tmp_file("region_size_analysis.txt");
        for (int i = 0; i < rpro.region_start_idx.size() - 1; ++i)
        {   
            tmp_size = rpro.region_end_idx[i] - rpro.region_start_idx[i];
            bucket[tmp_size] = bucket[tmp_size] + 1;
        }
        for (int i = 0; i < 82653; ++i)
        {
            tmp_file << bucket[i] << std::endl;
        }
        tmp_file.close();

        std::cout << "- Max Region Size:" << max_size << std::endl;*/
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
           if (kmer.find_first_of("9") == std::string::npos)
           {
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
        sa_bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(kmer_size,region_kmer,nthreads,gammaFactor);
        //---------------------------------------------------------

        T1.Stop();
        printf("- Building SA Region Index Using BBhash Finished (%f seconds)\n",T1.GetTime());

        std::ofstream region_file("sa/ref_suffix_region.txt");
        T0.Reset();     T0.Start();
        printf("- Starting To Search Region ...\n");
        uint64_t idx;
        uint64_t max_idx = 0;
        for (unsigned int i = 0; i < rpro.region_start_idx.size(); ++i)
        {
            idx = sa_bphf->lookup(region_kmer[i]);
            if (idx > max_idx)
            {
                max_idx = idx;
            }
            region_file << rpro.region_start_idx[i] << "\t" << rpro.region_end_idx[i] << "\t" << region_kmer[i] << "\t" << idx << std::endl;
            rpro.rkmer_idx.push_back(idx);
        }
        std::vector<size_t> region_bphf_idx(max_idx + 1);
        for (int i = 0; i < rpro.region_start_idx.size(); ++i)
        {
            idx = sa_bphf->lookup(region_kmer[i]);
            region_bphf_idx[idx] = i;
        }

        region_file.close();
        T0.Stop();
        printf("- Searching Regions Finished (%f seconds)\n",T0.GetTime());

        T0.Reset();     T0.Start();
        printf("- Starting To Compute Ref Hashcode ...\n");
        Compute_All_Ref_Code();
        T0.Stop();
        printf("- Computing Hashcode Of Ref Using SH (%f seconds)\n",T0.GetTime()); 

        T0.Reset();     T0.Start();
        printf("- Starting To Merge Duplicated Hashcode In Regions ...\n");
        Merge_Duplicated_In_Region();
        T0.Stop();
        printf("- Merging Duplicated Hashcode In Regions Finished (%f seconds)\n",T0.GetTime());

        T0.Reset();     T0.Start();
        //save sa region profile to disk using cereal
        {
            std::ofstream region_file("sa/ref_suffix_region.bin");
            cereal::BinaryOutputArchive ar(region_file);
            ar(rpro.rkmer_idx,rpro.region_start_idx,rpro.region_end_idx,rpro.code_bucket_idx);
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

    void REAL_TYPE_to_String(std::string &seq,REAL_TYPE *d,bool is_rc)
    {
        seq = "";
        if (!is_rc)
        {
            for(int i=0;i < KMER_SIZE;i++)
            {
                seq.push_back(itoic_table[(int)d[i]]);
            }
        }else
        {
            for(int i=0;i < KMER_SIZE;i++)
            {
                seq.push_back(rc_itoic_table[(int)d[dim -1 -i]]);
            }
        }
    }
    
    void Get_Read_Region()
    {
        // save read region to disk
        std::ofstream read_1_region_file;
        std::ofstream rc_read_1_region_file;

        std::ofstream read_2_region_file;
        std::ofstream rc_read_2_region_file;

        system("rm tmp/*");
        read_1_buff.Initialize (read_size,dim);
        read_1_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_1,Points::point_type::mapping);

        read_1_region_file.open(PAIR_1_REGION_FILE,std::ios::binary | std::ios::app);
        rc_read_1_region_file.open(PAIR_1_RC_REGION_FILE,std::ios::binary | std::ios::app);
        if(!read_1_region_file.is_open())
        {
            perror(PAIR_1_REGION_FILE);
        }

        read_2_buff.Initialize (read_size,dim);
        read_2_buff.Initialize_MemoryMapped(INPUT_READ_FILE_NAME_2,Points::point_type::mapping);

        read_2_region_file.open(PAIR_2_REGION_FILE,std::ios::binary);
        rc_read_2_region_file.open(PAIR_2_RC_REGION_FILE,std::ios::binary);
        if(!read_2_region_file.is_open())
        {
            perror(PAIR_2_REGION_FILE);
        }


        std::vector<size_t> region_bphf_idx;
        {
            std::ifstream region_bphf_idx_file("sa/region_bphf_idx.bin");
            cereal::BinaryInputArchive ar(region_bphf_idx_file);
            ar(region_bphf_idx);
        }
      
        std::vector<uint64_t> forward_read_1_region(read_size);  
        std::vector<uint64_t> rc_read_1_region(read_size);    
        std::vector<uint64_t> forward_read_2_region(read_size);  
        std::vector<uint64_t> rc_read_2_region(read_size);  
        std::string kmer;
        std::string rc_kmer;
        uint64_t rkmer_idx;
        uint64_t loc_in_rpro = 0;

    #ifdef USE_PARALLELIZATION
        #pragma omp parallel for
    #endif
        for (int i = 0; i < read_size; ++i)
        {     
            REAL_TYPE_to_String(kmer,read_1_buff.d[i],false);

            rkmer_idx = sa_bphf->lookup(kmer);
            forward_read_1_region[i] = region_bphf_idx[rkmer_idx];

            REAL_TYPE_to_String(rc_kmer,read_1_buff.d[i],true);
            rkmer_idx = sa_bphf->lookup(rc_kmer);
            rc_read_1_region[i] = region_bphf_idx[rkmer_idx];


            REAL_TYPE_to_String(kmer,read_2_buff.d[i],false);

            rkmer_idx = sa_bphf->lookup(kmer);
            forward_read_2_region[i] = region_bphf_idx[rkmer_idx];

            REAL_TYPE_to_String(rc_kmer,read_2_buff.d[i],true);
            rkmer_idx = sa_bphf->lookup(rc_kmer);
            rc_read_2_region[i] = region_bphf_idx[rkmer_idx];
        }
        // save read region info 
        {
            cereal::BinaryOutputArchive ar(read_1_region_file);
            ar(forward_read_1_region);
        }

        {
            cereal::BinaryOutputArchive ar(read_2_region_file);
            ar(forward_read_2_region);
        }

        {
            cereal::BinaryOutputArchive ar_rc(rc_read_1_region_file);
            ar_rc(rc_read_1_region);
        }

        {
            cereal::BinaryOutputArchive ar_rc(rc_read_2_region_file);
            ar_rc(rc_read_2_region);
        }
        
        read_1_region_file.close();
        rc_read_1_region_file.close();

        read_2_region_file.close();
        rc_read_2_region_file.close();
    }
    std::pair<int,int> Mapping_Process(size_t read_region, std::vector<unsigned long> &reduced_region_code,bitset<BCODE_LEN> bCodeRead,int pair)
    {   
        bitset<BCODE_LEN> bCodeRef;
        // delete
        std::string res_code;

        int region_size = 0;
        int min_dis = 0;
        int min_dis_idx = 0;
        int dist = 0;

        region_size = rpro.region_end_idx[read_region] - rpro.region_start_idx[read_region];
        if(read_region < rpro.region_start_idx.size() && region_size > 0)
        {
            min_dis = 1000;
            for (unsigned int j = 0; j < region_size; ++j)
            {    
                bCodeRef = reduced_region_code[rpro.region_start_idx[read_region] + j];
                dist = Compute_HD(bCodeRead, bCodeRef);
                
                if (min_dis > dist)
                {
                    min_dis_idx = j;
                    min_dis = dist;
                //    res_code = bCodeRef.to_string();
                }
            } 
        }else
        {
            min_dis = 1000;
            min_dis_idx = -1;
        }
        if(pair == PAIR_1)
            rcode_of_read_1.push_back(res_code);
        else
            rcode_of_read_2.push_back(res_code);

        std::pair<int,int> para(min_dis,min_dis_idx);
        return para;
    }

    void Load_Ref_Info()
    {
        {
            std::ifstream ref_name_file(REF_NAME_FILE);
            cereal::BinaryInputArchive ar_ref_name(ref_name_file);
            ar_ref_name(ref_name);
        }

        {
            std::ifstream loc_to_ref_file(MERGE_REF_POS_FILE);
            cereal::BinaryInputArchive ar_ref_pos(loc_to_ref_file);
            ar_ref_pos(loc_to_ref);
        }
    }

    int Hash_Mapping_with_SA()
    {
        std::ifstream read_1_region_file;
        std::ifstream rc_read_1_region_file;

        std::ifstream read_2_region_file;
        std::ifstream rc_read_2_region_file;


        std::vector<uint64_t> read_1_region;
        std::vector<uint64_t> rc_read_1_region;
        std::vector<uint64_t> res_read_1_region;

        std::vector<uint64_t> read_2_region;
        std::vector<uint64_t> rc_read_2_region;
        std::vector<uint64_t> res_read_2_region;

        std::ofstream tmp_loc_1_file;
        std::ofstream tmp_dis_1_file;
        std::ofstream read_1_code_file;
        std::ofstream res_read_1_region_file;

        std::ofstream tmp_loc_2_file;
        std::ofstream tmp_dis_2_file;
        std::ofstream read_2_code_file;
        std::ofstream res_read_2_region_file;

        // load ref code 
        std::ifstream irinfo_file(REDUCED_REGION_INFO_FILE);
        irinfo_file.seekg(0,irinfo_file.end);
        int size = irinfo_file.tellg() / sizeof(unsigned long);
        std::vector<unsigned long> reduced_region_code;
        reduced_region_code.resize(size);
        irinfo_file.seekg(0,irinfo_file.beg);
        irinfo_file.read(reinterpret_cast<char*>(&reduced_region_code[0]),size * sizeof(unsigned long));

        is_read_1_rev.resize(read_size,false);
        is_read_2_rev.resize(read_size,true);

        Stopwatch T0("");
        T0.Reset();     T0.Start();

        read_1_region_file.open(PAIR_1_REGION_FILE);
        rc_read_1_region_file.open(PAIR_1_RC_REGION_FILE);

        read_2_region_file.open(PAIR_2_REGION_FILE);
        rc_read_2_region_file.open(PAIR_2_RC_REGION_FILE);

        {
            cereal::BinaryInputArchive ar_read(read_1_region_file);
            ar_read(read_1_region);
        }

        {
            cereal::BinaryInputArchive ar_read(read_2_region_file);
            ar_read(read_2_region);
        }

        {
            cereal::BinaryInputArchive ar_read(rc_read_1_region_file);
            ar_read(rc_read_1_region);
        }

        {
            cereal::BinaryInputArchive ar_read(rc_read_2_region_file);
            ar_read(rc_read_2_region);
        }

        Load_Ref_Info();

        T0.Stop();
        printf("- Loading Reads' Info: (%f seconds)\n",T0.GetTime());

        T0.Reset();     T0.Start();

        REAL_TYPE *rc_read_1 = new REAL_TYPE [dim];
        REAL_TYPE *rc_read_2 = new REAL_TYPE [dim];
        bitset<BCODE_LEN> bCodeRead_1,bCodeRead_2;
        std::vector<bitset<BCODE_LEN>> bCodeRead_1_vec,bCodeRead_2_vec;
        mapped_res mres_1,mres_2;
        std::pair<int,int> min_res_1,min_res_2;
        uint64_t region_1 = 0,region_2 = 0;

    //    std:cout << "sarry size:" << sarry.con.size() << std::endl;
        std::vector<unsigned long> bCodeRef_vec;
    #ifdef USE_PARALLELIZATION
        #pragma omp parallel for
    #endif
        for (int i = 0; i < read_size; ++i)
        {
            reverse_complete(read_2_buff.d[i],rc_read_2);
            src_sh.Compute_BCode<REAL_TYPE>(read_1_buff.d[i], bCodeRead_1);
            src_sh.Compute_BCode<REAL_TYPE>(rc_read_2, bCodeRead_2);

         //   Mapping_Process
            region_1 = read_1_region[i];
            region_2 = rc_read_2_region[i];
            min_res_1 = std::move(Mapping_Process(region_1,reduced_region_code,bCodeRead_1,PAIR_1));
            min_res_2 = std::move(Mapping_Process(region_2,reduced_region_code,bCodeRead_2,PAIR_2));
            if (min_res_1.first > 2 && min_res_2.first > 2)
            {
                reverse_complete(read_1_buff.d[i],rc_read_1);
                src_sh.Compute_BCode<REAL_TYPE>(rc_read_1, bCodeRead_1);
                src_sh.Compute_BCode<REAL_TYPE>(read_2_buff.d[i], bCodeRead_2);
                region_1 = rc_read_1_region[i];
                region_2 = read_2_region[i];
                min_res_1 = std::move(Mapping_Process(region_1,reduced_region_code,bCodeRead_1,PAIR_1));
                min_res_2 = std::move(Mapping_Process(region_2,reduced_region_code,bCodeRead_2,PAIR_2));
                is_read_1_rev[i] = true;
                is_read_2_rev[i] = false;
            }

            bCodeRead_1_vec.push_back(bCodeRead_1);
            bCodeRead_2_vec.push_back(bCodeRead_2);

            mres_1.min_code_idx.push_back(min_res_1.second);
            mres_1.min_dis.push_back(min_res_1.first);
            res_read_1_region.push_back(region_1);

            mres_2.min_code_idx.push_back(min_res_2.second);
            mres_2.min_dis.push_back(min_res_2.first);
            res_read_2_region.push_back(region_2);
        }
        T0.Stop();
        printf("- Total Mapping Time: (%f seconds)\n",T0.GetTime());   


        tmp_loc_1_file.open(PAIR_1_LOC_FILE);
        tmp_dis_1_file.open(PAIR_1_DIS_FILE);

        read_1_code_file.open("tmp/read_code_1.bin");
        res_read_1_region_file.open(PAIR_1_RES_REGION_FILE);

        tmp_loc_2_file.open(PAIR_2_LOC_FILE);
        tmp_dis_2_file.open(PAIR_2_DIS_FILE);

        read_2_code_file.open("tmp/read_code_2.bin");
        res_read_2_region_file.open(PAIR_2_RES_REGION_FILE);

         //save to disk using cereal
        {
            cereal::BinaryOutputArchive ar_read(read_1_code_file);
            cereal::BinaryOutputArchive ar_loc(tmp_loc_1_file);
            cereal::BinaryOutputArchive ar_dis(tmp_dis_1_file);

            cereal::BinaryOutputArchive ar_region(res_read_1_region_file);
            
            ar_read(CEREAL_NVP(bCodeRead_1_vec));
            ar_region(CEREAL_NVP(res_read_1_region));
            ar_loc(CEREAL_NVP(mres_1.min_code_idx));
            ar_dis(CEREAL_NVP(mres_1.min_dis));
        }

        {
            cereal::BinaryOutputArchive ar_read(read_2_code_file);
            cereal::BinaryOutputArchive ar_loc(tmp_loc_2_file);
            cereal::BinaryOutputArchive ar_dis(tmp_dis_2_file);

            cereal::BinaryOutputArchive ar_region(res_read_2_region_file);
            
            ar_read(CEREAL_NVP(bCodeRead_2_vec));
            ar_region(CEREAL_NVP(res_read_2_region));
            ar_loc(CEREAL_NVP(mres_2.min_code_idx));
            ar_dis(CEREAL_NVP(mres_2.min_dis));
        }
        irinfo_file.close();
    }


    void Output_Result(int pair)
    {

        std::ofstream output;
        std::ifstream tmp_loc;
        std::ifstream tmp_dis;
        std::ifstream tmp_code;

        std::ifstream bucket_file(REDUCED_REGION_CODE_BUCKET_FILE);
        std::ifstream read_region_file;
        REAL_TYPE **read;
        std::vector<bool> is_read_rev;
        std::vector<std::string> ref_code;
        if (pair == PAIR_1)
        {
            tmp_loc.open(PAIR_1_LOC_FILE);
            tmp_dis.open(PAIR_1_DIS_FILE);
            tmp_code.open("tmp/read_code_1.bin");
            output.open(PAIR_1_RES_FILE);
            read_region_file.open(PAIR_1_RES_REGION_FILE);
            is_read_rev = std::move(is_read_1_rev);
            read = &read_1_buff.d[0];
            ref_code = std::move(rcode_of_read_1);
        }else if (pair == PAIR_2)
        {
            tmp_loc.open(PAIR_2_LOC_FILE);
            tmp_dis.open(PAIR_2_DIS_FILE);
            tmp_code.open("tmp/read_code_2.bin");
            output.open(PAIR_2_RES_FILE);
            read_region_file.open(PAIR_2_RES_REGION_FILE);
            is_read_rev = std::move(is_read_2_rev);
            read = &read_2_buff.d[0];
            ref_code = std::move(rcode_of_read_2);
        }


        // load mapping result from disk using cereal
        mapped_res mres;
        {
            cereal::BinaryInputArchive ar_loc(tmp_loc);
            cereal::BinaryInputArchive ar_dis(tmp_dis);
            ar_loc(CEREAL_NVP(mres.min_code_idx));
            ar_dis(CEREAL_NVP(mres.min_dis));
        };

        std::vector<bitset<BCODE_LEN>> read_code;
        {
            cereal::BinaryInputArchive ar_code(tmp_code);
            ar_code(CEREAL_NVP(read_code));
        }

        std::vector<uint64_t> read_region;
        {
            cereal::BinaryInputArchive ar_read(read_region_file);
            ar_read(read_region);
        }

        std::vector<std::string> read_name;
        {
            std::ifstream read_name_file(PAIR_1_NAME_FILE);
            cereal::BinaryInputArchive ar(read_name_file);
            ar(read_name);
        }

        std::vector<size_t> ref_loc_vec;
        bucket_file.seekg(0,bucket_file.end);
        int bucket_size = bucket_file.tellg() / sizeof(size_t);
        ref_loc_vec.resize(bucket_size);
        bucket_file.seekg(0,bucket_file.beg);
        bucket_file.read(reinterpret_cast<char*>(&ref_loc_vec[0]),bucket_size * sizeof(size_t));

        Stopwatch T0("");
        T0.Reset();     T0.Start();

        size_t ref_loc = 0;
        int loc_size = 0;
        std::string tmp_read; 
        std::string rev_read;
        std::vector<std::vector<int>> ref_of_read(read_size);
        std::vector<int> ref_of_one_read;

    #ifdef USE_PARALLELIZATION
        #pragma omp parallel for
    #endif
        for (unsigned int qIndex = 0;qIndex < mres.min_code_idx.size();++qIndex)
        {      
            output << '>' << qIndex + 1  << ":" << mres.min_dis[qIndex] << ":" ;//<< std::endl;
            output << read_code[qIndex] << "\t" << "-" << is_read_rev[qIndex] << "\t" <<read_name[qIndex] << std::endl;
            output << "= " ;

            if (is_read_rev[qIndex])
            {
                reverse_complete_ictos(read[qIndex],rev_read);
                output << rev_read << std::endl;
            }else
            {
                for (int i = 0; i < dim; ++i)
                {
                    output << (char)itos_table[(int8_t)read[qIndex][i]];
                }
                output << std::endl;
            }
            
            if (mres.min_code_idx[qIndex] >= 0)
            {
                if (mres.min_code_idx[qIndex] == rpro.code_bucket_idx[read_region[qIndex]].size() - 1)
                {
                    loc_size = rpro.code_bucket_idx[read_region[qIndex] + 1][0] 
                        - rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex]];
                }else
                {
                    loc_size = rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex] + 1] 
                        - rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex]];
                }
    
                ref_of_one_read.clear();
                for (unsigned int i = 0; i < loc_size; ++i)
                {
                    ref_loc = ref_loc_vec[rpro.code_bucket_idx[read_region[qIndex]][mres.min_code_idx[qIndex]] + i];
                    if (ref_loc > 0)
                    {

                        //------ seek with sa------
                        output << "+ "; 
                        for (int j = 0; j < dim; ++j)
                        {
                            output << (char)ictos_table[(int8_t)ref_string[ref_loc + j]];
                        }
                        output << " " << ref_name[loc_to_ref[ref_loc]] << " " << loc_to_ref[ref_loc] << std::endl;
                        ref_of_one_read.push_back(loc_to_ref[ref_loc]);
                    }
                }
                ref_of_read[qIndex] = ref_of_one_read;
            }
        }

        if (pair == PAIR_1)
        {
            {
                std::ofstream ref_of_read_file("dataset/ref_of_read_1.bin");
                cereal::BinaryOutputArchive ar(ref_of_read_file);
                ar(ref_of_read);
            }
        }else if (pair == PAIR_2)
        {
            {
                std::ofstream ref_of_read_file("dataset/ref_of_read_2.bin");
                cereal::BinaryOutputArchive ar(ref_of_read_file);
                ar(ref_of_read);
            }
        }
        output.close();
        bucket_file.close();
        T0.Stop();
        printf("- Save Mapping Results To Disk Finished (%f seconds)\n",T0.GetTime() );
    }
};