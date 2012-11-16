/*
 *  Flexbar.h
 *
 *  Created on: Jul 31, 2012
 *      Author: jtr
 */


#ifndef FLEXBAR_FLEXBAR_H_
#define FLEXBAR_FLEXBAR_H_

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include <tbb/pipeline.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>
#include <seqan/arg_parse.h>

#include "Enums.h"
#include "SequencingRead.h"
#include "AdapterLoader.h"
#include "SequenceInputFilter.h"
#include "MultiplexedInputFilter.h"
#include "MultiplexedOutputFilter.h"
#include "MultiplexedAlignmentFilter.h"
#include "AlignmentAlgorithm.h"


std::string getFlexbarBanner(){
	
	std::string banner = "\n";
	
	banner += "                   ________          __              \n";
	banner += "                  / ____/ /__  _  __/ /_  ____ ______\n";
	banner += "                 / /_  / / _ \\| |/_/ __ \\/ __ `/ ___/\n";
	banner += "                / __/ / /  __/>  </ /_/ / /_/ / /    \n";
	banner += "               /_/   /_/\\___/_/|_/_.___/\\__,_/_/     \n\n";
	
	banner += "Flexbar - flexible barcode detection and adapter removal, version 2.21\n";
	banner += "Developed at the Berlin Institute for Medical Systems Biology, GPLv3\n\n";
	
	return banner;
}


void prepareHelpMessage(seqan::ArgumentParser &parser){
	
	using namespace seqan;
	
	setVersion(parser, "2.21");
	setDate(parser, "October 23, 2012");
	
	std::string banner = getFlexbarBanner();
	banner += "Available on: sf.net/projects/theflexibleadap\n";
	
	setShortDescription(parser, banner);
	
	addUsageLine(parser, "\\fB-t\\fP target \\fB-f\\fP format \\fB-s\\fP reads { \\fB-b\\fP barcodes | \\fB-a\\fP adapters } [options]");
	
	addOption(parser, ArgParseOption("n", "threads", "Number of threads", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("t", "target", "Prefix for output file names", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("s", "source", "Input file with reads, that may contain barcodes", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("p", "source2", "Second input file for paired read scenario", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("f", "format", "Input format of reads: csfasta, csfastq, fasta, fastq-sanger, fastq-solexa, fastq-i1.3, fastq-i1.5, fastq-i1.8 (illumina)", ArgParseArgument::STRING));
	
	addSection(parser, "Barcode detection");
	addOption(parser, ArgParseOption("b",  "barcodes", "Fasta file with barcodes, specify (br) to use seperate barcode reads", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("br", "barcode-reads", "Fasta or fastq file with barcode reads, if barcodes not within reads", ArgParseArgument::INPUTFILE));
	addOption(parser, ArgParseOption("be", "barcode-trim-end", "Type of barcoding within source reads, see section trim-end types", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("bt", "barcode-threshold", "Allowed mismatches and indels per 10 bases for barcode", ArgParseArgument::DOUBLE));
	addOption(parser, ArgParseOption("bo", "barcode-min-overlap", "Minimum overlap for barcodes (default is length of first barcode)", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bm", "barcode-match", "Match score", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bi", "barcode-mismatch", "Mismatch score", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bg", "barcode-gap-cost", "Gap score", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("bv", "barcode-remove", "Remove barcodes within reads based on barcoding parameters"));
	
	addSection(parser, "Adapter removal");
	addOption(parser, ArgParseOption("a",  "adapters", "Fasta file with adapters, or barcodes to remove within reads", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("as", "adapter-seq", "Single adapter sequence as alternative to adapters option", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("ae", "adapter-trim-end", "Type of alignment for removal, see section trim-end types", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("at", "adapter-threshold", "Allowed mismatches and indels per 10 bases for adapter", ArgParseArgument::DOUBLE));
	addOption(parser, ArgParseOption("ao", "adapter-min-overlap", "Minimum overlap of adapter and read in base pairs", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("am", "adapter-match", "Match score", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("ai", "adapter-mismatch", "Mismatch score", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("ag", "adapter-gap-cost", "Gap score", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("aa", "adapter-no-adapt", "Do not treat parameter min-overlap as adaptive measure, see doc"));
	addOption(parser, ArgParseOption("ab", "adapter-banded-algo", "Use banded version of alignment algorithm"));
	
	addSection(parser, "Filtering");
	addOption(parser, ArgParseOption("m", "min-readlength", "Minimum read length to remain after removal", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("u", "max-uncalled", "Allowed uncalled bases (N or .) in reads", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("x", "pre-trim-front", "Trim specified number of bases on 5' end of reads before removal", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("y", "pre-trim-back", "Trim specified number of bases on 3' end of reads before removal", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("q", "pre-trim-phred", "Trim reads from 3' end until specified or higher quality reached", ArgParseArgument::INTEGER));
	addOption(parser, ArgParseOption("d", "no-length-dist", "Prevent length distribution for each read output file"));
	addOption(parser, ArgParseOption("r", "removal-tag", "Tag reads for which adapter or barcode is removed"));
	addOption(parser, ArgParseOption("i", "short-read-file", "File for omitted reads being shorter min length", ArgParseArgument::STRING));
	addOption(parser, ArgParseOption("l", "log-level", "Print alignments for all or modified reads.", ArgParseArgument::STRING));
	
	addSection(parser, "Trim-end types");
	addText(parser._toolDoc, "ANY: longer part of read remains", false);
	addText(parser._toolDoc, "LEFT: adapters that align <= read end position", false);
	addText(parser._toolDoc, "RIGHT: adapters that align >= read start position", false);
	addText(parser._toolDoc, "LEFT_TAIL: considers first n (adapter length) bases", false);
	addText(parser._toolDoc, "RIGHT_TAIL: considers last n bases of reads", false);
	
	setRequired(parser, "t");
	setRequired(parser, "f");
	setRequired(parser, "s");
	
	hideOption(parser, "adapter-banded-algo");
	hideOption(parser, "short-read-file");
	
	setDefaultValue(parser, "threads", "1");
	
	setDefaultValue(parser, "barcode-trim-end", "ANY");
	setDefaultValue(parser, "barcode-threshold", "1.0");
	setDefaultValue(parser, "barcode-match", "1");
	setDefaultValue(parser, "barcode-mismatch", "-1");
	setDefaultValue(parser, "barcode-gap-cost", "-7");
	
	setDefaultValue(parser, "adapter-trim-end", "RIGHT");
	setDefaultValue(parser, "adapter-threshold", "3.0");
	setDefaultValue(parser, "adapter-min-overlap", "8");
	setDefaultValue(parser, "adapter-match", "1");
	setDefaultValue(parser, "adapter-mismatch", "-1");
	setDefaultValue(parser, "adapter-gap-cost", "-7");
	
	setDefaultValue(parser, "min-readlength", "18");
	setDefaultValue(parser, "max-uncalled", "0");
	
	// setMinValue(parser, "threads", "1");
	setValidValues(parser, "log-level", "ALL MOD TAB");
	
	addTextSection(parser, "EXAMPLES");
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-t\\fP target \\fB-f\\fP csfastq    \\fB-s\\fP reads.csfastq \\fB-a\\fP adapters.fasta", false);
	addText(parser._toolDoc, "\\fBflexbar\\fP \\fB-t\\fP target \\fB-f\\fP fastq-i1.3 \\fB-s\\fP reads.fastq   \\fB-b\\fP bar.fasta \\fB-a\\fP adap.fasta");
}


struct Options{
	
	bool g, l;
	int c;
    
    Options(){
		g = false;
		l = false;
        c = 40;
    }
};


void printLocalTime(){
	time_t t_current;
	time(&t_current);
	printf("Local time:            %s\n", asctime(localtime(&t_current)));
}


void parseCommandLine(Options &o, seqan::ArgumentParser &parser, int argc, char const ** argv){
	
	using seqan::ArgumentParser;
	
	ArgumentParser::ParseResult res = parse(parser, argc, argv);
	
	if(res != ArgumentParser::PARSE_OK)	exit(res == ArgumentParser::PARSE_ERROR);
	
	std::cout << getFlexbarBanner();
	printLocalTime();
}


void printComputationTime(time_t t_start){
	using namespace std;
	
	time_t t_end;
	time(&t_end);
	
	int totalTime = int(difftime(t_end, t_start));
	
	int hours   = div(totalTime, 3600).quot;
	int rest    = div(totalTime, 3600).rem;
	int minutes = div(rest, 60).quot;
	int seconds = div(rest, 60).rem;
	
	cout << "Computation time:  ";
	if(hours > 0)                               cout << hours     << " h ";
	if(hours > 0 || minutes > 0)                cout << minutes   << " min ";
	if(hours > 0 || minutes > 0 || seconds > 0) cout << seconds   << " sec\n\n\n";
	else                                        cout              << "< 1 sec\n\n\n";
}


void printCompletedMessage(flexbar::BarcodeDetect barDetect, flexbar::AdapterRemoval adapRm){
	
	using namespace std;
	using namespace flexbar;
	
	cout << "Flexbar completed ";
	
	if(barDetect != BOFF)                   cout << "barcode";
	if(barDetect == REMOVE_WITHIN_READ)     cout << " removal within reads";
	if(barDetect == WITHIN_READ)            cout << " detection within reads";
	if(barDetect == BARCODE_READ)           cout << " detection with seperate reads";
	if(barDetect != BOFF && adapRm != AOFF) cout << " and ";
	if(adapRm != AOFF)                      cout << "adapter removal";
	cout << "." << endl << endl;
}


#endif /* FLEXBAR_FLEXBAR_H_ */
