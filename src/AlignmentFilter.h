/*
 *  AlignmentFilter.h
 *
 *  Created on: Jun 29, 2010
 *      Author: mat
 *
 *  Revised on: Aug 3, 2012
 *      Author: jtr
 */

#ifndef FLEXBAR_ALIGNMENTFILTER_H_
#define FLEXBAR_ALIGNMENTFILTER_H_

#include <sstream>
#include <string>

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include "Enums.h"
#include "SequencingRead.h"
#include "SequenceConverter.h"
#include "AdapterLoader.h"


/* This class does the actual alignment via the passed Algorithm. It has
an internal vector of adapter sequences and will align each adapter to the
read. The one with the highest score will be used for the final alignment. */

template <class TString, class TIDString, class TAlgorithm >
class AlignmentFilter {

private:
	
	tbb::atomic<unsigned long> m_sumLength;
	tbb::atomic<flexbar::TrimEnd> m_trimEnd;
	
	tbb::atomic<bool> m_writeTag;
	tbb::atomic<int> m_minOverlapLength;
	tbb::atomic<int> m_maxOverlapLength;
	tbb::atomic<int> m_modified;
	
	tbb::concurrent_vector<TAdapter> *m_adapters;
	tbb::concurrent_vector<unsigned long> *m_rmOverlaps;
	
	flexbar::LogLevel m_verb;
	flexbar::FileFormat m_format;
	
	bool m_adaptiveOverlap;
	int m_minOverlap, m_match, m_mismatch, m_gapCost;
	float m_threshold;
	
	TAlgorithm *algo;
	
public:
	
	AlignmentFilter(tbb::concurrent_vector<TAdapter> *adapters, int match, int mismatch, int gapCost, bool writeTag, flexbar::TrimEnd end, flexbar::LogLevel logLevel, flexbar::FileFormat format){
		
		m_verb      = logLevel;
		m_trimEnd   = end;
		m_adapters  = adapters;
		m_match     = match;
		m_mismatch  = mismatch;
		m_gapCost   = gapCost;
		m_format    = format;
		
		m_writeTag   = writeTag;
		m_minOverlap = 10;
		m_threshold  = 1;
		m_modified   = 0;
		m_sumLength  = 0;
		
		m_maxOverlapLength = 0;
		m_minOverlapLength = 1000;
		
		algo = new TAlgorithm(m_match, m_mismatch, m_gapCost);
		
		m_rmOverlaps = new tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
	};
	
	
	void setThreshold(float threshold){
		m_threshold = threshold;
	}
	
	
	float getThreshold(){
		return m_threshold;
	}
	
	
	void setMinOverlap(int minOverlap, bool adaptiveOverlap){
		m_minOverlap      = minOverlap;
		m_adaptiveOverlap = adaptiveOverlap;
	}
	
	
	int getMinOverlap(){
		return m_minOverlap;
	}
	
	
	int getMinOverlapLength(){
		return m_minOverlapLength;
	}
	
	
	int getMaxOverlapLength(){
		return m_maxOverlapLength;
	}
	
	
	int getMedianOverlapLength(){
		unsigned long nValues = 0, halfValues = 0, cumValues = 0;
		
		for(int i = 0; i <= flexbar::MAX_READLENGTH; i++){ nValues += m_rmOverlaps->at(i); }
		
		if(nValues % 2 == 0) halfValues = nValues / 2;
		else                 halfValues = (nValues - 1) / 2;
		
		for(int i = 0; i <= flexbar::MAX_READLENGTH; i++){
			cumValues += m_rmOverlaps->at(i);
			
			if(cumValues >= halfValues) return i;
		}
        
        return 0;
	}
	
	
	int getMeanOverlapLength(){
		if(m_modified > 0) return m_sumLength / m_modified;
		else return -1;
	}
	
	
	int getNrModifiedReads(){
		return m_modified;
	}
	
	
	virtual ~AlignmentFilter(){};
	
	
	// function aligns all adapters to read, which will be cut if specified
	
	int align(void* item, bool postCutRead){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		
		
		SequencingRead<TString, TIDString> &myRead = *static_cast< SequencingRead<TString, TIDString>* >(item);
		
		int fmismatches, fgaps, foverlapLength;
		int fstartPos, fstartPosA, fstartPosS, fendPos, fendPosS, fendPosA;
		
		int fIndex   = -1;
		int scoreMax = -10000;
		
		float fallowedErrors;
		
		stringstream ss;
		TString read, quality, finalAliStr;
		
		
		TString readTag = myRead.getSequenceTag();
		
		// only align suffixes in colorspace
		switch(m_format){
			case CSFASTQ:
				read    = suffix<TString>(myRead.getSequence(), 2);
				quality = suffix<TString>(myRead.getQuality(), 1);
				break;
			case FASTQ:
				read    = myRead.getSequence();
				quality = myRead.getQuality();
				break;
			case CSFASTA:
				read    = suffix<TString>(myRead.getSequence(), 2);
				quality = "";
				break;
			case FASTA:
				read    = myRead.getSequence();
				quality = "";
				break;
		}
		
		int readLength = length(read);
		
		TString sequence = read;
		TString squality = quality;
		
		// iterate over all passed adapter sequences and remember the best one
		for(unsigned int i = 0; i < m_adapters->size(); ++i){
			
			TString currentAdapterTag = m_adapters->at(i).first->getSequenceTag();
			TString currentAdapter    = m_adapters->at(i).first->getSequence();
			
			int adapterLength = length(currentAdapter);
			
			
			// in tail mode trim read before alignment
			if(m_trimEnd == LEFT_TAIL || m_trimEnd == RIGHT_TAIL){
				
				// align only overlapping part of sequence and adapter - speedup
				if(adapterLength < readLength){
					
					if(m_trimEnd == LEFT_TAIL){
						sequence = prefix<TString>(read, adapterLength);
						
						if(m_format == FASTQ)   squality = prefix<TString>(quality, adapterLength);
						if(m_format == CSFASTQ) squality = prefix<TString>(quality, adapterLength - 1);
					}
					else {
						sequence = suffix<TString>(read, readLength - adapterLength);
						
						if(m_format == FASTQ)   squality = suffix<TString>(quality, readLength - adapterLength);
						if(m_format == CSFASTQ) squality = suffix<TString>(quality, readLength - adapterLength - 1);
					}
					
					if(m_verb == ALL || m_verb == MOD)
					ss << "Read trimmed to adapter length:  " << readTag << endl << endl;
				}
				else if(adapterLength > readLength){
					
					if(m_trimEnd == LEFT_TAIL){
						currentAdapter = suffix<TString>(currentAdapter, adapterLength - readLength);
					}
					else currentAdapter = prefix<TString>(currentAdapter, readLength);
					
					if(m_verb == ALL || m_verb == MOD)
					ss << "Adapter trimmed to read length:  " << currentAdapterTag << endl << endl;
				}
			}
			
			
			int startPos, endPos, startPosA, endPosA, startPosS, endPosS, aliScore, mismatches, gaps;
			stringstream aliString;
			
			// align currentAdapter via passed algorithm
			algo->align(currentAdapter, sequence, gaps, mismatches, startPos, endPos, startPosA, endPosA, startPosS, endPosS, aliScore, aliString, m_trimEnd);
			
			
			int overlapLength = endPos - startPos;
			
			float allowedErrors = m_threshold * overlapLength / 10.0f;
			float madeErrors    = static_cast<float>(mismatches + gaps);
			
			int minOverlapValue = m_minOverlap;
			
			
			bool validAli = true;
			
			if(((m_trimEnd == RIGHT_TAIL || m_trimEnd == RIGHT) && startPosA < startPosS) ||
			   ((m_trimEnd == LEFT_TAIL  || m_trimEnd == LEFT)  && endPosA > endPosS)     ||
			     overlapLength <= 0){
				
				validAli = false;
			}
			
			
			if(validAli && m_adaptiveOverlap && minOverlapValue <= adapterLength){
				
				if(((m_trimEnd == RIGHT_TAIL || m_trimEnd == RIGHT || m_trimEnd == ANY) &&
				    (startPosA >= startPosS && startPosA < endPosS && endPosA > endPosS))
				   ||
				   ((m_trimEnd == LEFT_TAIL || m_trimEnd == LEFT || m_trimEnd == ANY) &&
				    (endPosA <= endPosS && endPosA > startPosS && startPosA < startPosS))){
					
					minOverlapValue = overlapLength;
				}
			}
			
			
			// check if alignment is valid and score is max as well as if number of errors and overlap length are allowed
			if(validAli && aliScore > scoreMax && madeErrors <= allowedErrors && overlapLength >= minOverlapValue){
				
				fIndex      = i;
				scoreMax    = aliScore;
				fstartPos   = startPos;
				fstartPosA  = startPosA;
				fstartPosS  = startPosS;
				fendPos     = endPos;
				fendPosA    = endPosA;
				fendPosS    = endPosS;
				
				foverlapLength = overlapLength;
				
				if(m_verb == MOD || m_verb == ALL || m_verb == TAB){
					fgaps          = gaps;
					fmismatches    = mismatches;
					finalAliStr    = aliString.str();
					fallowedErrors = allowedErrors;
				}
			}
		}
		
		++fIndex;
		
		
		/// cut read depending on best adapter match ///
		
		bool rMod = false;
		
		if(postCutRead && fIndex > 0){
			
			TrimEnd trimEnd = m_trimEnd;
			
			if(trimEnd == ANY){
				
				if(fstartPosA <= fstartPosS && fendPosS <= fendPosA){
					rMod = true;
					myRead.setSequence("");
					myRead.setQuality("");
				}
				else if(fstartPosA - fstartPosS >= fendPosS - fendPosA){
					trimEnd = RIGHT;
					ss << "Trimming from right..." << endl;
				}
				else {
					trimEnd = LEFT;
					ss << "Trimming from left..." << endl;
				}
			}
			
			int rCutPos;
			int fstartPosMod = fstartPos;
			
			switch(trimEnd){
				case LEFT_TAIL:
					sequence = read;
					squality = quality;
				
				// now do the same as trimming from LEFT
				case LEFT:
					if(fendPosA <= fendPosS){
						
						rMod    = true;
						rCutPos = fendPos;
						
						// translate alignment end pos to read index
						if(fstartPosS > 0)       rCutPos -= fstartPosS;
						if(rCutPos > readLength) rCutPos  = readLength;
						
						erase(sequence, 0, rCutPos);
						
						if(m_format == FASTA || m_format == FASTQ){
							myRead.setSequence(sequence);
						}
						// prefix of read in colorspace
						else {
							TString colRead = "T";
							append(colRead, sequence);
							myRead.setSequence(colRead);
						}
						
						if(m_format == FASTQ || m_format == CSFASTQ){
							erase(squality, 0, rCutPos);
							myRead.setQuality(squality);
						}
					}
					break;
					
				case RIGHT_TAIL:
					sequence  = read;
					squality  = quality;
					// adjust cut pos to complete read length
					fstartPosMod = readLength - foverlapLength;
					
				case RIGHT:
					if(fstartPosA >= fstartPosS){
						
						rMod    = true;
						rCutPos = fstartPosMod;
						
						erase(sequence, rCutPos, readLength);
						
						if(m_format == FASTA || m_format == FASTQ){
							myRead.setSequence(sequence);
							
							if(m_format == FASTQ){
								erase(squality, rCutPos, readLength);
								myRead.setQuality(squality);
							}
						}
						// append original TX prefix in colorspace
						else {
							TString result = prefix(myRead.getSequence(), 2);
							append(result, sequence);
							myRead.setSequence(result);
							
							if(m_format == CSFASTQ){
								erase(squality, rCutPos, readLength);
								
								TString qresult = prefix(myRead.getQuality(), 1);
								append(qresult, squality);
								myRead.setQuality(qresult);
							}
						}
					}
					break;
                case ANY:;
			}
			
			
			if(rMod){
				
				if(m_writeTag){
					TString newTag = readTag;
					append(newTag, " - Flexbar removal");
					myRead.setSequenceTag(newTag);
				}
				
				++m_modified;
				
				// count for each adapter how often it was removed
				m_adapters->at(fIndex - 1).second++;
				
				// output alignment
				if(m_verb == MOD || m_verb == ALL){
					ss << "Sequence Removal"             << endl;
					ss << "read-tag:        " << readTag << endl;
					ss << "adapter-tag:     " << m_adapters->at(fIndex - 1).first->getSequenceTag() << endl;
					ss << "read:            " << read    << endl;
					
					if(m_format == FASTQ || m_format == CSFASTQ)
					ss << "quality:         " << quality    << endl;
					ss << "read pos:        " << fstartPosS << "-" << fendPosS << endl;
					ss << "adapter pos:     " << fstartPosA << "-" << fendPosA << endl;
					ss << "overlap:         " << foverlapLength       << endl;
					ss << "score:           " << scoreMax             << endl;
					ss << "indels:          " << fgaps                << endl;
					ss << "mismatches:      " << fmismatches          << endl;
					ss << "allowed errors:  " << fallowedErrors       << endl;
					ss << "remaining seq:   " << myRead.getSequence() << endl;
					
					if(m_format == FASTQ || m_format == CSFASTQ)
					ss << "remaining qual:  " << myRead.getQuality()  << endl;
					
					ss << "Alignment:" << endl << endl << finalAliStr;
				}
				
				if(m_verb == TAB)
				ss << readTag << "\t" << m_adapters->at(fIndex-1).first->getSequenceTag() << "\t" << fstartPosA << "\t"
				   << fendPosA << "\t" << foverlapLength << "\t" << fmismatches << "\t" << fgaps << "\t" << fallowedErrors << endl;
				
				// calculate max and min overlap length
				if(foverlapLength > m_maxOverlapLength) m_maxOverlapLength = foverlapLength;
				if(foverlapLength < m_minOverlapLength) m_minOverlapLength = foverlapLength;
				
				// calculate sum for mean, and occurrences for median
				m_sumLength += foverlapLength;
				m_rmOverlaps->at(foverlapLength)++;
			}
		}
		
		// read not modified
		if(! rMod && m_verb == ALL){
			ss << "No valid alignment" << endl;
			ss << "read-tag:  " << readTag << endl;
			ss << "read:      " << read    << endl << endl << endl;
		}
		
		// bundeled output for multithreading
		if(m_verb != NONE && fIndex != -1) cout << ss.str();
		
		return fIndex;
	};
	
};

#endif /* FLEXBAR_ALIGNMENTFILTER_H_ */
