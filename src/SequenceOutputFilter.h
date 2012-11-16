/*
 * SequenceOutputFilter.h
 *
 *  Created on: Apr 12, 2010
 *      Author: mat
 */

#ifndef FLEXBAR_SEQUENCEOUTPUTFILTER_H_
#define FLEXBAR_SEQUENCEOUTPUTFILTER_H_

#include <map>
#include <fstream>

#include "Enums.h"
#include "SequencingRead.h"


/* This class writes sequencing reads in specified format to file. */

template <typename TString, typename TIDString>
class SequenceOutputFilter {

private:
	
	std::map<unsigned int, unsigned int> m_lengthDist;
	
	std::fstream m_targetStream;
	std::string m_filePath;
	
	unsigned long m_countGood;
	
	flexbar::FileFormat m_format;
	
public:
	
	SequenceOutputFilter(std::string filePath, flexbar::FileFormat format){
		using namespace std;
		
		m_format    = format;
		m_filePath  = filePath;
		m_countGood = 0;
		
		m_targetStream.open(m_filePath.c_str(), ios_base::out | ios_base::binary);
		
		if(!m_targetStream.is_open()){
			cerr << "Error opening file: " << m_filePath << endl;
		}
	};
	
	
	std::string getFileName(){
		return m_filePath;
	}
	
	
	virtual ~SequenceOutputFilter(){
		m_targetStream.close();
	};
	
	
	void writeLengthDist(){
		using namespace std;
		
		string fname = m_filePath + ".lengthdist";
		fstream lstream;
		
		lstream.open(fname.c_str(), ios_base::out | ios_base::binary);
		
		if(!lstream.is_open()){
			cerr << "Error opening File: " << fname << endl;
		}
		else {
			map<unsigned int, unsigned int>::iterator iter;
			lstream << "Readlength\tCount" << endl;
			for( iter = m_lengthDist.begin(); iter != m_lengthDist.end(); iter++ ) {
				lstream << iter->first << "\t" << iter->second << endl;
			}
			lstream.close();
		}
	}
	
	
	std::string getFastString(SequencingRead<TString, TIDString> *myRead){
		
		using namespace std;
		using namespace flexbar;
		
		stringstream ss;
		
		switch(m_format){
			case FASTQ:   ss << "@" << myRead->getSequenceTag() << endl << myRead->getSequence() << endl << "+" << endl << myRead->getQuality() << endl; break;
			case FASTA:   ss << ">" << myRead->getSequenceTag() << endl << myRead->getSequence() << endl; break;
			case CSFASTQ: ss << "@" << myRead->getSequenceTag() << endl << myRead->getSequence() << endl << "+" << endl << myRead->getQuality() << endl; break;
			case CSFASTA: ss << ">" << myRead->getSequenceTag() << endl << myRead->getSequence() << endl; break;
		}
		return ss.str();
	}
	
	
	unsigned long getNrGoodReads(){
		return m_countGood;
	}
	
	
	void *writeRead(void* item){
		if(item){
			SequencingRead<TString, TIDString> *myRead = static_cast< SequencingRead<TString, TIDString>* >(item);
			
			unsigned int readLength = length(myRead->getSequence());
			
			if(m_targetStream.is_open() && m_targetStream.good()){
				
				++m_countGood;
				
				// store read length distribution
				std::map<unsigned int, unsigned int>::iterator it;
				it = m_lengthDist.find(readLength);
				
				if(it != m_lengthDist.end()) it->second++;
				else                         m_lengthDist[readLength] = 1;
				
				m_targetStream << getFastString(myRead);
				m_targetStream.flush();
				
				return NULL;
			}
			else {
				std::cerr << "Error writing target file!" << m_filePath << std::endl;
				exit(1);
			}
		}
		return NULL;
	}
	
};

#endif /* FLEXBAR_SEQUENCEOUTPUTFILTER_H_ */
