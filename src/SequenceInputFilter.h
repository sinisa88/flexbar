/*
 * SequenceInputFilter.h
 *
 *  Created on: Apr 12, 2010
 *      Author: mat
 *
 *  Revised on: Jul 16, 2012
 *      Author: jtr
 */

#ifndef FLEXBAR_SEQUENCEINPUTFILTER_H_
#define FLEXBAR_SEQUENCEINPUTFILTER_H_

#include <string>
#include <fstream>
#include <stdexcept>

#include <tbb/pipeline.h>

#include "Enums.h"
#include "SequencingRead.h"


/* This class reads a FASTA or FASTQ file and returns a SequencingRead for each valid sequence. */

template <class TString, class TIDString>
class SequenceInputFilter : public tbb::filter {

private:
    
	std::fstream fstrm;
	std::string m_filePath;
	flexbar::FileFormat m_format;
	TIDString m_nextTag;
	
	int m_prePhredTrim, m_maxUncalled;
	unsigned int m_preTrimBegin, m_preTrimEnd, m_minLength;
	unsigned long m_nrReads, m_count_lowPhred, m_count_shortPhred;
	
public:
	
	SequenceInputFilter(std::string filePath, flexbar::FileFormat format, int maxUncalled) : filter(serial_in_order){
		
		m_filePath    = filePath;
		m_format      = format;
		m_maxUncalled = maxUncalled;
		m_nextTag     = "";
		
		m_nrReads          = 0;
		m_count_lowPhred   = 0;
		m_count_shortPhred = 0;
		m_preTrimEnd       = 0;
		m_prePhredTrim     = 0;
		
		fstrm.open(m_filePath.c_str(), std::ios_base::in);
		
		if(!fstrm.is_open()){
			std::cout << "Error opening File: " << m_filePath << std::endl;
			exit(0);
		}
	};
	
	
	void setPreTrimBegin(int nr){
		if(nr > 0) this->m_preTrimBegin = nr;
	}
	
	
	void setPreTrimEnd(int nr){
		if(nr > 0) this->m_preTrimEnd = nr;
	}
	
	
	void setPrePhredTrim(int nr, flexbar::QualityType qtype, bool showText){
		
		using namespace std;
		using namespace flexbar;
		
		if(nr > 0) this->m_prePhredTrim = nr;
		else       return;
		
		switch(qtype){
			case SANGER:{
				this->m_prePhredTrim += 33;
				
				if(showText)
				cout << "Trimming reads from 3' end to phred quality " << nr << " (" << m_prePhredTrim << ")" << endl;
				break;
			}
			case SOLEXA:{
				this->m_prePhredTrim += 59;
				
				if(showText)
				cout << "Trimming reads from 3' end to phred quality " << nr << " (" << m_prePhredTrim << ")" << endl;
				break;
			}
			case ILLUMINA13:{
				this->m_prePhredTrim += 64;
				
				if(showText)
				cout << "Trimming reads from 3' end to phred quality " << nr << " (" << m_prePhredTrim << ")" << endl;
				break;
			}
		}
	}
	
	
	void setMinReadLength(unsigned int minReadLength){
		m_minLength = minReadLength;
	}
	
	
	unsigned long getNrLowPhredReads(){
		return m_count_lowPhred;
	}
	
	
	unsigned long getNrShortPhredReads(){
		return m_count_shortPhred;
	}
	
	
	unsigned long getNrProcessedReads(){
		return m_nrReads;
	}
	
	
	TString readLine(){
		char line[4096];
		TString text = "";
		
		if(fstrm.good()){
			fstrm.getline(line, 4096);
			text = line;
		}
		return text;
	}
	
	
	/* This is the core method for reading and parsing FASTA/FASTQ input.
	@return: single SequencingRead<TString, TIDString> or NULL if no more reads in file or error occured.
	@param valid is set to false if read is shorter than minimum readlength, e.g. due to phred trimming. */
	
	void* getRead(bool &isValid, bool &isUncalled){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		using seqan::length;
		
		SequencingRead<TString, TIDString> *myRead = NULL;
		
		TString source = "", quality = "", dummy = "";
		TIDString sequence_tag = "";
		
		if(!fstrm.eof()){
			
			isValid    = true;
			isUncalled = false;
			
			try{
				
				if(m_format == FASTA || m_format == CSFASTA){
					// FastA parsing
					
					// tag line will be read in previous iteration
					if(m_nextTag == "") sequence_tag = readLine();
					else                sequence_tag = m_nextTag;
					
					if(length(sequence_tag) > 0){
						if(seqan::isNotEqual(getValue(sequence_tag, 0), '>')){
							stringstream error;
							error << "Incorrect FASTA entry, missing new > line. Input: " << sequence_tag << endl;
							throw runtime_error(error.str());
						}
						else sequence_tag = suffix(sequence_tag, 1);
						
						if(length(sequence_tag) == 0){
							stringstream error;
							error << "Incorrect FASTA entry, missing readname after >. Input: " << sequence_tag << endl;
							throw runtime_error(error.str());
						}
					}
					else return NULL;
					
					
					source = readLine();
					
					if(length(source) < 1){
						stringstream error;
						error << "Warning, found sequence tag without read! Tag: " << sequence_tag << endl;
						throw runtime_error(error.str());
					}
					
					
					m_nextTag = readLine();
					
					// fasta files that have sequences splitted over several lines
					while(fstrm.good() && seqan::isNotEqual(getValue(m_nextTag, 0), '>')){
						append(source, m_nextTag);
						m_nextTag = readLine();
					}
					
					
					isUncalled = isUncalledSequence(source);
					
					if(this->m_preTrimBegin > 1 && length(source) > m_preTrimBegin){
						source = suffix(source, m_preTrimBegin);
						
						if(m_format == CSFASTA){
							TString colRead = "T";
							append(colRead, source);
							source = colRead;
						}
					}
					
					if(this->m_preTrimEnd > 1 && length(source) > m_preTrimEnd){
						source = prefix(source, length(source) - m_preTrimEnd);
					}
					
					
					myRead = new SequencingRead<TString, TIDString>(source, sequence_tag);
					
					++m_nrReads;
				}
				
				else{  // FastQ parsing
					
					source = readLine();
					
					if(length(source) > 0){
						if(seqan::isNotEqual(getValue(source, 0), '@')){
							stringstream error;
							error << "Incorrect FASTQ entry, missing new @ line. Input: " << source << endl;
							throw runtime_error(error.str());
						}
						else sequence_tag = suffix(source, 1);
						
						if(length(sequence_tag) == 0){
							stringstream error;
							error << "Incorrect FASTQ entry, missing readname after @. Input: " << source << endl;
							throw runtime_error(error.str());
						}
					}
					else return NULL;
					
					source = readLine();
					
					if(length(source) < 1){
						stringstream error;
						error << "Warning, found sequence tag without read! Tag: " << sequence_tag << endl;
						throw runtime_error(error.str());
					}
					
					
					dummy = readLine();
					if(length(dummy) == 0 || seqan::isNotEqual(getValue(dummy, 0), '+')){
							stringstream error;
							error << "Incorrect FASTQ entry, missing + line. Readname: " << sequence_tag << endl;
							throw runtime_error(error.str());
					}
					
					quality = readLine();
					
					// in case CSFASTQ format has same quality and read length it will be trimmed
					if(m_format == CSFASTQ){
						if(length(quality) == length(source)){
							quality = suffix(quality, 1);
						}
					}
					
					if(length(quality) < 1){
						stringstream error;
						error << "Warning, found sequence without quality values! Tag: " << sequence_tag << endl;
						throw runtime_error(error.str());
					}
					
					
					isUncalled = isUncalledSequence(source);
					
					if(this->m_preTrimBegin > 1 && length(source) > m_preTrimBegin){
						source = suffix(source, m_preTrimBegin);
						
						if(m_format == CSFASTQ){
							TString colRead = "T";
							append(colRead, source);
							source = colRead;
							
							quality = suffix(quality, m_preTrimBegin - 1);
						}
						else quality = suffix(quality, m_preTrimBegin);
					}
					
					if(this->m_preTrimEnd > 1 && length(source) > m_preTrimEnd){
						source = prefix(source, length(source) - m_preTrimEnd);
						
						if(m_format == CSFASTQ){
							quality = prefix(quality, (length(source) - m_preTrimEnd) - 1);
						}
						else quality = prefix(quality, length(source) - m_preTrimEnd);
					}
					
					
					// filtering based on phred quality
					if(this->m_prePhredTrim != 0){
						typename seqan::Iterator<TString >::Type it    = seqan::begin(quality);
						typename seqan::Iterator<TString >::Type itEnd = seqan::end(quality);
						
						--itEnd;
						
						unsigned int n = length(quality);
						bool isTooShort = n < this->m_minLength;
						
						bool nChanged = false;
						
						while(itEnd != it){
							if(static_cast<int>(*itEnd) >= this->m_prePhredTrim) break;
							--n;
							--itEnd;
							
							if(! nChanged){
								m_count_lowPhred++;
								nChanged = true;
							}
						}
						
						source = prefix(source, n);
						
						if(m_format == CSFASTQ) --n;
						quality = prefix(quality, n);
						
						if(! isTooShort && n < this->m_minLength){
							m_count_shortPhred++;
							isValid = false;
						}
					}
					
					
					myRead = new SequencingRead<TString, TIDString>(source, sequence_tag, quality);
					
					++m_nrReads;
				}
				
				return myRead;
			}
			catch(ios_base::failure &failure){
				cout << failure.what() << endl;
				fstrm.close();
				return NULL;
			}
		}
		else {  // end of stream reached
			return NULL;
		}
	}
	
	
	// returns TRUE if read contains too many uncalled bases
	bool isUncalledSequence(TString source){
		int n = 0;
		
		typename seqan::Iterator<TString >::Type it    = seqan::begin(source);
		typename seqan::Iterator<TString >::Type itEnd = seqan::end(source);
		
		while(it != itEnd){
			 if(*it == '.' || *it == 'N') n++;
			 ++it;
		}
		
		if(n <= m_maxUncalled) return false;
		else                   return true;
	}
 	
	
	// override
	void* operator()(void *){
		
		bool valid    = true;
		bool uncalled = false;
		
		return getRead(valid, uncalled);
	}
	
	
	virtual ~SequenceInputFilter(){};
	
};

#endif /* FLEXBAR_SEQUENCEINPUTFILTER_H_ */
