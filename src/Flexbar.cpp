/*=======================================================================
Name:           Flexbar.cpp
Authors:        Matthias Dodt and Johannes Roehr

Description:    Flexbar - flexible barcode detection and adapter removal
Version:        2.21
Copyright:      GPL version 3

SeqAn library:  post release 1.3.1, revision 13058 on October 9, 2012
TBB   library:  version 4.0 update 5, stable release June 13, 2012
========================================================================*/


#include "Flexbar.h"

// #include <seqan/find.h>


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::ArgumentParser;
	using seqan::CharString;
	using seqan::UnicodeString;
	
	
	// CharString haystack = "ATGGATTGCG";
	// CharString needle   = "ATGCAT";
	// 
	// seqan::Finder<CharString> finder(haystack);
	// seqan::Pattern<CharString, seqan::DPSearch<seqan::SimpleScore, seqan::FindInfix> > pattern(needle, seqan::SimpleScore(0, -1, -7));
	// 
	// while (find(finder, pattern, -2)){
	//     while (findBegin(finder, pattern, getScore(pattern))){
	//         cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << endl;
	//         //cout << end(finder) << endl; //',' << position(pattern) << endl;
	// 	}
	// }
	// 
	// cout << "------" << endl;
	// 
	// clear(finder);
	// seqan::Pattern<CharString, seqan::AbndmAlgo > pattern2(needle, -2);
	// 
	// //seqan::Score<int> sc(0,-3,-2);  // = scoringScheme(pattern2);
	// //setScoringScheme(pattern2, sc);
	// 
	// while (find(finder, pattern2)){
	//     while (findBegin(finder, pattern2, getScore(pattern2))){
	//         cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << endl;
	// 	}
	// }
	
	
	time_t t_start;
	time(&t_start);
	
	Options o;
	ArgumentParser parser("flexbar");
	
	prepareHelpMessage(parser);
	parseCommandLine(o, parser, argc, argv);
	
	
	/// declaration and initialisation of variables ///
	
	string readsFile, readsFile2, targetName, shortReadFile, format;
	string adapterFile, barcodeFile, barReadsFile, adapterSeq, alignMode, log_level;
	
	int maxUncalled, min_readLen, a_min_overlap, b_min_overlap, nThreads;
	int match, mismatch, gapCost, b_match, b_mismatch, b_gapCost;
	
	float a_threshold, b_threshold;
	
	tbb::concurrent_vector<TAdapter> adapters, barcodes;
	
	bool isColorSpace   = false;
	bool useAdapterFile = false;
	bool useRemovalTag  = false;
	
	int cutLen_begin  = 0;
	int cutLen_end    = 0;
	int phred_preQual = 0;
	
	string a_trim_end;
	string b_trim_end;
	
	TrimEnd end, b_end;
	FileFormat fformat;
	QualityType qual;
	RunType runType;
	
	LogLevel logLevel       = NONE;
	BarcodeDetect barDetect = BOFF;
	AdapterRemoval adapRm   = AOFF;
	
	
	/// parsing of program options ///
	
	getOptionValue(targetName, parser, "target");
	cout << "Target name:           " << targetName << endl;
	
	getOptionValue(format, parser, "format");
	cout << "File format:           " << format;
	
	if(format == "fasta"){
		fformat = FASTA;
		qual = SANGER;
	}
	else if(format == "fastq-sanger"){
		fformat = FASTQ;
		qual = SANGER;
	}
	else if(format == "fastq-solexa"){
		fformat = FASTQ;
		qual = SOLEXA;
	}
	else if(format == "fastq-i1.3" || format == "fastq-illumina1.3"){
		fformat = FASTQ;
		qual = ILLUMINA13;
	}
	else if(format == "fastq-i1.5" || format == "fastq-illumina1.5"){
		fformat = FASTQ;
		qual = ILLUMINA13;
	}
	else if(format == "fastq-i1.8" || format == "fastq-illumina1.8"){
		fformat = FASTQ;
		qual = SANGER;
	}
	else if(format == "csfasta"){
		fformat = CSFASTA;
		qual = SANGER;
		isColorSpace = true;
	}
	else if(format == "csfastq"){
		fformat = CSFASTQ;
		qual = SANGER;
		cout << " (sanger quality scaling)";
		isColorSpace = true;
	}
	else{
		cerr << "Specified input file format is unknown!"<< endl << endl;
		return 1;
	}
	cout << endl;
	
	
	getOptionValue(readsFile, parser, "source");
	cout << "Source file:           " << readsFile << endl;
	runType = SINGLE;
	
	if(isSet(parser, "source2")){
		getOptionValue(readsFile2, parser, "source2");
		cout << "Source file2:          " << readsFile2 << "   (paired run)" << endl;
		runType = PAIRED;
	}
	
	
	if(isSet(parser, "barcodes")){
		
		getOptionValue(barcodeFile, parser, "barcodes");
		cout << "Barcode file:          " << barcodeFile << endl;
		
		if(isSet(parser, "barcode-reads")){
			getOptionValue(barReadsFile, parser, "barcode-reads");
			cout << "Barcode reads file:    " << barReadsFile << endl;
			
			barDetect = BARCODE_READ;
			
			if(runType == PAIRED)      runType = PAIRED_BARCODED;
			else if(runType == SINGLE) runType = SINGLE_BARCODED;
		}
		else {
			barDetect = WITHIN_READ;
			runType   = SINGLE_BARCODED_WITHIN_READ;
		}
	}
	
	
	if(isSet(parser, "adapters")){
		getOptionValue(adapterFile, parser, "adapters");
		cout << "Adapter file:          " << adapterFile << endl;
		adapRm = NORMAL;
		useAdapterFile = true;
	}
	else if(isSet(parser, "adapter-seq")){
		getOptionValue(adapterSeq, parser, "adapter-seq");
		adapRm = NORMAL;
	}
	else if(barDetect == BOFF){
		cerr << endl << "Specify either a barcodes file, adapters file or single adapter sequence!" << endl << endl;
		return 1;
	}
	
	
	if(isSet(parser, "short-read-file")){
		getOptionValue(shortReadFile, parser, "short-read-file");
		cout << "Short reads file:      " << shortReadFile << endl;
	}
	else{ shortReadFile = readsFile + ".omitted"; }
	
	cout << endl;
	
	
	if(barDetect == WITHIN_READ && isSet(parser, "barcode-remove")){
		barDetect = REMOVE_WITHIN_READ;
		cout << endl << "Removing barcodes in single alignment step based on barcoding options." << endl
		             << "Consider also to specify barcodes as adapters to use seperate parameters for removal." << endl << endl;
	}
	
	
	if(isSet(parser, "adapter-banded-algo")){
		alignMode = "BandedNW";
	} else { alignMode = "NW"; }
	
	if(isSet(parser, "removal-tag")) useRemovalTag = true;
	
	
	getOptionValue(nThreads, parser, "threads");
	cout << "threads:               " << nThreads << endl;
	
	getOptionValue(min_readLen, parser, "min-readlength");
	cout << "min-readlength:        " << min_readLen << endl;
	
	getOptionValue(maxUncalled, parser, "max-uncalled");
	cout << "max-uncalled:          " << maxUncalled << endl;
	
	if(isColorSpace) min_readLen++;
	
	
	if(isSet(parser, "pre-trim-phred")){
		getOptionValue(phred_preQual, parser, "pre-trim-phred");
		cout << "pre-trim-phred:        " << phred_preQual << endl;
	}
	
	if(isSet(parser, "pre-trim-front")){
		getOptionValue(cutLen_begin, parser, "pre-trim-front");
		cout << "pre-trim-front:        " << cutLen_begin << endl;
	}
	
	if(isSet(parser, "pre-trim-back")){
		getOptionValue(cutLen_end, parser, "pre-trim-back");
		cout << "pre-trim-back:       " << cutLen_end << endl;
	}
	
	
	/// barcode and adapter specific options ///
	
	if(barDetect != BOFF){
		
		getOptionValue(b_trim_end, parser, "barcode-trim-end");
		if     (b_trim_end == "LEFT")        b_end = LEFT;
		else if(b_trim_end == "RIGHT")       b_end = RIGHT;
		else if(b_trim_end == "ANY")         b_end = ANY;
		else if(b_trim_end == "LEFT_TAIL")   b_end = LEFT_TAIL;
		else if(b_trim_end == "RIGHT_TAIL")  b_end = RIGHT_TAIL;
		else{
			cerr << "Specified barcode trim-end is unknown!" << endl << endl;
			return 1;
		}
		cout << "barcode-trim-end:      " << b_trim_end << endl;
		
		if(isSet(parser, "barcode-min-overlap")){
			getOptionValue(b_min_overlap, parser, "barcode-min-overlap");
			cout << "barcode-min-overlap:   " << b_min_overlap << endl;
		}
		
		getOptionValue(b_threshold, parser, "barcode-threshold");
		cout << "barcode-threshold:     " << b_threshold << endl;
		
		
		getOptionValue(b_match, parser, "barcode-match");
		getOptionValue(b_mismatch, parser, "barcode-mismatch");
		getOptionValue(b_gapCost, parser, "barcode-gap-cost");
		
		cout << "barcode-match:        ";
		if(b_match >= 0) cout << " ";
		cout << b_match << endl;
		
		cout << "barcode-mismatch:     ";
		if(b_mismatch >= 0) cout << " ";
		cout << b_mismatch << endl;
		
		cout << "barcode-gap-cost:     ";
		if(b_gapCost >= 0) cout << " ";
		cout << b_gapCost << endl << endl;
	}
	
	
	if(adapRm != AOFF){
		
		getOptionValue(a_threshold, parser, "adapter-threshold");
		cout << "adapter-threshold:     " << a_threshold << endl;
		
		getOptionValue(a_min_overlap, parser, "adapter-min-overlap");
		cout << "adapter-min-overlap:   " << a_min_overlap << endl;
		
		getOptionValue(a_trim_end, parser, "adapter-trim-end");
		if     (a_trim_end == "LEFT")        end = LEFT;
		else if(a_trim_end == "RIGHT")       end = RIGHT;
		else if(a_trim_end == "ANY")         end = ANY;
		else if(a_trim_end == "LEFT_TAIL")   end = LEFT_TAIL;
		else if(a_trim_end == "RIGHT_TAIL")  end = RIGHT_TAIL;
		else {
			cerr << "Specified adapter trim-end is unknown!" << endl << endl;
			return 1;
		}
		cout << "adapter-trim-end:      " << a_trim_end;
		
		if(! isSet(parser, "adapter-no-adapt")){
			cout << "   (adaptive overlap)";
			adapRm = ADAPTIVE;
		}
		cout << endl;
		
		
		getOptionValue(match, parser, "adapter-match");
		getOptionValue(mismatch, parser, "adapter-mismatch");
		getOptionValue(gapCost, parser, "adapter-gap-cost");
		
		cout << "adapter-match:        ";
		if(match >= 0) cout << " ";
		cout << match << endl;
		
		cout << "adapter-mismatch:     ";
		if(mismatch >= 0) cout << " ";
		cout << mismatch << endl;
		
		cout << "adapter-gap-cost:     ";
		if(gapCost >= 0) cout << " ";
		cout << gapCost << endl << endl;
	}
	
	
	if(isSet(parser, "log-level")){
		getOptionValue(log_level, parser, "log-level");
		
		     if(log_level == "ALL") logLevel = ALL;
		else if(log_level == "TAB") logLevel = TAB;
		else if(log_level == "MOD") logLevel = MOD;
		else if(log_level != "NONE"){
			cerr << "Specified log-level is unknown!" << endl << endl;
			return 1;
		}
	}
	
	
	/// loading barcodes and adapters file ///
	
	if(barDetect != BOFF){
		
		tbb::task_scheduler_init init_serial(1);
		tbb::pipeline bpipeline;
		
		SequenceInputFilter<CharString,CharString > adapter_filter(barcodeFile, FASTA, maxUncalled);
		bpipeline.add_filter(adapter_filter);
		
		AdapterLoader<CharString, CharString> adapterLoader(fformat);
		bpipeline.add_filter(adapterLoader);
		
		bpipeline.run(1);
		barcodes = adapterLoader.getAdapters();
		
		adapterLoader.printAdapters("Barcode");
		
		if(barcodes.size() == 0){
			cerr << "No barcodes found in file!" << endl << endl;
			return 1;
		}
		
		if(! isSet(parser, "barcode-min-overlap")){
			b_min_overlap = length(barcodes.at(0).first->getSequence());
		}
	}
	
	
	if(adapRm != AOFF){
		
		AdapterLoader<CharString, CharString> adapterLoader(fformat);
		
		if(useAdapterFile){
			tbb::task_scheduler_init init_serial(1);
			tbb::pipeline prepipeline;
			
			SequenceInputFilter<CharString, CharString> adapter_filter(adapterFile, FASTA, maxUncalled);
			prepipeline.add_filter(adapter_filter);
			
			prepipeline.add_filter(adapterLoader);
			prepipeline.run(1);
			
			adapters = adapterLoader.getAdapters();
			
			if(adapters.size() == 0){
				cerr << "No adapters found in file!" << endl << endl;
				return 1;
			}
		}
		else {
			SequencingRead<CharString,CharString> *myRead;
			myRead = new SequencingRead<CharString, CharString>(adapterSeq, "cmdline");
			
			TAdapter adap;
			adap.first = myRead;
			adapters.push_back(adap);
			
			adapterLoader.setAdapters(adapters);
		}
		
		adapterLoader.printAdapters("Adapter");
	}
	
	
	/// main computation in parallel ///
	
	tbb::task_scheduler_init init_serial(nThreads);
	tbb::pipeline pipeline;
	
	// create file reading-writing stage
	MultiplexedInputFilter<CharString, CharString > input_filter(readsFile, fformat, maxUncalled, cutLen_begin, cutLen_end, min_readLen, phred_preQual, qual);
	
	if(barDetect == BARCODE_READ) input_filter.setBarcodeReadsFile(barReadsFile);
	if(runType == PAIRED || runType == PAIRED_BARCODED) input_filter.setPairedFile(readsFile2);
	
	pipeline.add_filter(input_filter);
	
	cout << endl << "Processing reads ...";
	
	if(logLevel != NONE)
	cout << endl << endl << "Generating " << log_level << " verbose output. This will slow down processing!" << endl << endl;
	
	
	MultiplexedAlignmentFilter<CharString, CharString, AlignmentAlgorithm<CharString> > filter(&adapters, &barcodes, a_threshold, b_threshold, a_min_overlap, b_min_overlap, b_end, match, mismatch, gapCost, b_match, b_mismatch, b_gapCost, end, barDetect, adapRm, useRemovalTag, logLevel, fformat);
	
	pipeline.add_filter(filter);
	
	// create file writing stage
	MultiplexedOutputFilter<CharString, CharString >
	output_filter(targetName, &adapters, &barcodes, fformat, min_readLen, shortReadFile, runType);
	
	pipeline.add_filter(output_filter);
	pipeline.run(nThreads);
	
	if(logLevel == TAB) cout << endl;
	cout << "done." << endl  << endl;
	printComputationTime(t_start);
	
	
	/// barcode identification and adapter removal statistics ///
	
	if(! isSet(parser, "no-length-dist")) output_filter.writeLengthDist();
	
	filter.printAdapterOverlapStats();
	
	unsigned long nReads     = input_filter.getNrProcessedReads();
	unsigned long nGoodReads = output_filter.getNrGoodReads();
	
	
	cout << "Filtering statistics " << endl;
	cout << "=========================="  << endl;
	cout << "Processed reads:                   " << nReads << endl;
	
	
	cout << "  skipped due to uncalled bases:   " << input_filter.getNrUncalledReads();
	
	if(runType == PAIRED || runType == PAIRED_BARCODED)
		cout << "   (discarded " << input_filter.getNrUncalledPairedReads() << " paired reads)" << endl;
	else cout << endl;
	
	
	if(isSet(parser, "pre-trim-phred"))
		cout << "  skipped due to low quality:      " << input_filter.getNrShortPhredReads()
	         << "   (trimmed " << input_filter.getNrLowPhredReads() << " reads)"  << endl;
	
	cout << "  leftover shorter min length:     " << output_filter.getNrShortReads() << endl;
	
	cout << "Discarded reads overall:           " << (nReads - nGoodReads) << endl;
	
	if(nReads > 0){
		cout << "Remaining reads:                   " << nGoodReads << "   (" << fixed <<
		        setprecision(2) << 100 * nGoodReads / nReads << "% of input reads)" << endl << endl;
	}
	
	output_filter.printFileSummary();
	
	if(adapRm != AOFF) output_filter.printAdapterRemovalStats();
	
	
	printCompletedMessage(barDetect, adapRm);
	
	return 0;
}

