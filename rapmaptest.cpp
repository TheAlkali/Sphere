#include<iostream>
#include "Stopwatch.hpp"
int main()
{
	Stopwatch T0("");
/*    	T0.Reset();     T0.Start();
	system("rapmap quasiindex -t ../reference/transcripts/Homo_sapiens.GRCh38.cdna.all.fa -i ../rapmap_res/ref_index");
	T0.Stop();
	std::cout << "rapmap index time:" << T0.GetTime() << std::endl;*/
	T0.Reset();     T0.Start();
	system("rapmap quasimap -i ../rapmap_res/ref_index -1 dataset/srrdata/small_SRR5337025_1.fastq -2 dataset/srrdata/small_SRR5337025_2.fastq -o ../rapmap_res/rapmap_mapped_reads.sam");
	T0.Stop();
	std::cout << "rapmap mapping time:" << T0.GetTime() << std::endl;
	return 0;
}
