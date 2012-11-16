//=====================================
// Name        : PairedreadFinder.cpp
// Author      : Matthias Dodt
// Copyright   : GPL V3
// Description : Extension to Flexbar
//=====================================


#include <iostream>
#include <seqan/misc/misc_cmdparser.h>
#include <tbb/pipeline.h>
#include <tbb/tick_count.h>
#include <tbb/task_scheduler_init.h>

#include "Enums.h"
#include "LookupFilter.h"
#include "HashmapFilter.h"
#include "SequenceInputFilter.h"


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace seqan;
	
	time_t start,end;
	double totalTime;

	time(&start);

	// Declare the supported options.
	string program = "PairedreadFinder, Version 1.01. This tool takes two fasta/q files and looks for matching readnames in both files.";
	stringstream ss;
	
	ss << program << " - Allowed options";

	CommandLineParser cmdParser(program.c_str());

	addOption(cmdParser,CommandLineOption("v","version","return program version",OptionType::Label));
	addOption(cmdParser,CommandLineOption("s1","source1", "input file 1", OptionType::String));
	addOption(cmdParser,CommandLineOption("s2","source2", "input file 2", OptionType::String));
	addOption(cmdParser,CommandLineOption("f","format", "input file format", OptionType::String));
	addOption(cmdParser,CommandLineOption("t1","target1", "target file 1", OptionType::String));
	addOption(cmdParser,CommandLineOption("t2","target2", "target file 2", OptionType::String));
	addOption(cmdParser,CommandLineOption("n","nr-threads", "nr of threads to use", (int)OptionType::Int,1));
	addOption(cmdParser,CommandLineOption("is","suffix-ignore", "nr of characters to ignore from the END of the readname (in case paired reads are named like /1 /2 it should be set to 2)", (int)OptionType::Int,0));
	addOption(cmdParser,CommandLineOption("ip","prefix-ignore", "nr of characters to ignore from the BEGINNING of the readname (in case paired reads are named like s_1.. s_2.. it should be set to 3)", (int)OptionType::Int,0));

	if(!parse(cmdParser,argc,argv,cerr)){
		cerr << "Error parsing command line. Check the mandatory parameters." << endl;
		help(cmdParser,cerr);
		return -1;
	}

	string source1,source2,target1,target2;

	getOptionValueLong(cmdParser,"source1",source1);
	if(length(source1)>0){
		cout << "source1 file was set to: " << source1 << ".\n";
	}
	else {
		cerr << "No source1 file specified!\n"<<endl;
		help(cmdParser,cerr);
		return -1;
	}

	getOptionValueLong(cmdParser,"source2",source2);
	if(length(source2)>0){
		cout << "source2 file was set to: " << source2 << ".\n";
	}
	else {
		cerr << "No source2 file specified!\n"<<endl;
		help(cmdParser,cerr);
		return -1;
	}

	getOptionValueLong(cmdParser,"target1",target1);
	if(length(target1)>0){
		cout << "target1 file prefix was set to: " << target1 << ".\n";
	}
	else {
		cerr << "No target1 file prefix specified!\n"<<endl;
		help(cmdParser,cerr);
		return -1;
	}

	getOptionValueLong(cmdParser,"target2",target2);
	if(length(target2)>0){
		cout << "target2 file was set to: " << target2 << ".\n";
	}
	else {
		cerr << "No target2 file specified!\n"<<endl;
		help(cmdParser,cerr);
		return -1;
	}

	int nrThreads = 1;
	getOptionValueLong(cmdParser,"nr-threads",nrThreads);
	cout << "Using " << nrThreads << " threads." << endl << endl;

	unsigned int suffixIgnore=0;
	getOptionValueLong(cmdParser,"suffix-ignore",suffixIgnore);
	if(suffixIgnore>0){
		cout << "Ignoring last " << suffixIgnore << " characters of readname for finding it's mate" << endl;
	}
	
	unsigned int prefixIgnore=0;
	getOptionValueLong(cmdParser,"prefix-ignore",prefixIgnore);
	if(prefixIgnore>0){
		cout << "Ignoring first " << prefixIgnore << " characters of readname for finding it's mate" << endl;
	}
	
	flexbar::FileFormat fformat = flexbar::FASTA;
	string format = "";
	
	if(getOptionValueLong(cmdParser,"format",format)){
		if (format=="fasta")
			{
				cout << "File format was set to FASTA" << endl;
			}
			else {
				if (format=="fastq"){
					fformat = flexbar::FASTQ;
					cout << "File format was set to FASTQ" << endl;
				}
				else{
					if (format=="csfasta")
					{
						fformat = flexbar::FASTA;
						cout << "File format was set to colorspace FASTA" << endl;

					}
					else {
						if (format=="csfastq"){
							fformat = flexbar::FASTQ;

							cout << "File format was set to colorspace FASTQ" << endl;
						}
						else {
							cout << "Please specify an input-file format!" << endl;
							help(cmdParser,cerr);
							return -1;
						}
					}
				}
			}
	}
	else {
		cerr << "Please specify an input-file format! "<< endl << endl;
		//help(cmdParser,cerr);
		return -1;
	}

	cout << "Single reads will be stored in files named <targetName>.single" << endl<< endl << "Starting caching reads... " << endl;
	try{
		// Start task scheduler
		tbb::task_scheduler_init init_serial( nrThreads );

		// Create the pipeline
		tbb::pipeline pipeline;

		// Create file-reading stage and add it to the pipeline
		SequenceInputFilter<seqan::String<char>,seqan::String<char> > input_filter( source1, fformat );
		pipeline.add_filter( input_filter );

		HashmapFilter<seqan::String<char>, seqan::String<char> > *hashmapFilter = new HashmapFilter<seqan::String<char>, seqan::String<char> >(target1);
		hashmapFilter->setSuffixIgnore(suffixIgnore);
		hashmapFilter->setPrefixIgnore(prefixIgnore);
		pipeline.add_filter( *hashmapFilter );

		// Run the pipeline
		//tbb::tick_count t0 = tbb::tick_count::now();
		pipeline.run( nrThreads );
		//tbb::tick_count t1 = tbb::tick_count::now();

		cout << "Reads cached. Processing 2nd file..." << endl;

		tbb::pipeline pipeline2;

		SequenceInputFilter<seqan::String<char>,seqan::String<char> > input_filter2( source2, fformat );
		pipeline2.add_filter( input_filter2 );

		LookupFilter<seqan::String<char>,seqan::String<char> > lookupFilter(hashmapFilter, target1, target2);
		lookupFilter.setSuffixIgnore(suffixIgnore);
		lookupFilter.setPrefixIgnore(prefixIgnore);
		pipeline2.add_filter(lookupFilter);

		pipeline2.run( nrThreads );

		cout << "Writing omitted reads of hashmap... " << endl;
		hashmapFilter->writeOmittedReads();
		cout << "Finished. " << endl<< endl;
		cout << "Nr. of pairs  : " << lookupFilter.getNrPairs() << endl;
		cout << "Nr. of singles: " << lookupFilter.getNrSingles() + hashmapFilter->getNrSingles() << endl;

		cout << "( Nr. single1 reads: " << hashmapFilter->getNrSingles() << ")" << endl;
		cout << "( Nr. single2 reads: " << lookupFilter.getNrSingles() << ")" << endl;
		cout << "=====================================" << endl;
		cout << "Nr. reads used: " << lookupFilter.getNrSingles() + hashmapFilter->getNrSingles() + lookupFilter.getNrPairs() * 2 << endl<< endl;
		time(&end);
		totalTime = difftime(end,start);
		cout << "Runtime: " << div(static_cast<int>(totalTime),60).quot << " minutes "<< div(static_cast<int>(totalTime), 60).rem << " seconds. " << endl;
	}
	catch(exception &e){
		cout << "Error running pipeline: " << e.what() << endl;
		exit(0);
	}
	
	return 0;
}

