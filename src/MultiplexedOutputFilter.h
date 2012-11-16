/*
 * MultiplexedInputFilter.h
 *
 *  Created on: May 31, 2011
 *      Author: mat
 *
 *  Revised on: Jul 16, 2012
 *      Author: jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_
#define FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>

#include "Enums.h"
#include "MultiplexedRead.h"
#include "SequenceOutputFilter.h"
#include "OutputFileStruct.h"
#include "AdapterLoader.h"


/* This class will process a MultiplexedRead and write it to a file depending on the runtype (single-end, paired-end and/or barcoded). */

template <typename TString, typename TIDString>
class MultiplexedOutputFilter : public tbb::filter {

private:
	
	unsigned int m_mapsize, m_minLength;
	long m_cnt_total;
	
	std::string m_filePath;
	flexbar::FileFormat m_format;
	flexbar::RunType m_runType;
	
	typedef SequenceOutputFilter<TString, TIDString> TOutputFilter;
	typedef OutputFileStruct<TString, TIDString> filters;
	
	filters *m_outputMap;
	
	tbb::concurrent_vector<TAdapter> *m_adapters, *m_barcodes;
	
public:
	
	MultiplexedOutputFilter(std::string filePath, tbb::concurrent_vector<TAdapter> *adapters, tbb::concurrent_vector<TAdapter> *barcodes, flexbar::FileFormat format, int min_read_length, std::string shortFilename, flexbar::RunType runType) : filter(serial_in_order){
		
		using namespace std;
		using namespace flexbar;
		
		m_filePath  = filePath;
		m_format    = format;
		m_adapters  = adapters;
		m_barcodes  = barcodes;
		m_runType   = runType;
		m_minLength = min_read_length;
		m_mapsize   = 0;
		
		
		switch(m_runType){
			
			case PAIRED_BARCODED:{
				
				m_mapsize = barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				stringstream ss;
				
				for(unsigned int i = 0; i < barcodes->size(); ++i){
					
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_1" << toFormatString(m_format);
					TOutputFilter *filter1 = new TOutputFilter(ss.str(), m_format);
					
					ss.str("");
					ss.clear();

					// single
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_1_single" << toFormatString(m_format);
					TOutputFilter *single1 = new TOutputFilter(ss.str(), m_format);

					ss.str("");
					ss.clear();

					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_2"<< toFormatString(m_format);
					TOutputFilter *filter2 = new TOutputFilter(ss.str(), m_format);
					
					ss.str("");
					ss.clear();

					// single reads
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << "_2_single"<< toFormatString(m_format);
					TOutputFilter *single2 = new TOutputFilter(ss.str(), m_format);
					
					ss.str("");
					ss.clear();

					filters f;
					f.f1 = filter1;
					f.f2 = filter2;
					f.single1 = single1;
					f.single2 = single2;

					m_outputMap[i + 1] = f;
				}
				
				ss << filePath << "_barcode_unassigned_1"<< toFormatString(m_format);
				TOutputFilter *filter3 = new SequenceOutputFilter<TString, TIDString>(ss.str(), m_format);

				ss.str("");
				ss.clear();
				
				ss << shortFilename << "_barcode_unassigned_1"<< toFormatString(m_format);
				
				ss.str("");
				ss.clear();
				
				
				// single reads
				ss << filePath << "_barcode_unassigned_1_single"<< toFormatString(m_format);
				TOutputFilter *single3 = new SequenceOutputFilter<TString, TIDString>(ss.str(), m_format);

				ss.str("");
				ss.clear();

				ss << shortFilename << "_barcode_unassigned_1_single"<< toFormatString(m_format);
				
				ss.str("");
				ss.clear();
				
				ss << filePath << "_barcode_unassigned_2"<< toFormatString(m_format);
				TOutputFilter *filter4 = new SequenceOutputFilter<TString, TIDString>(ss.str(), m_format);
				
				ss.str("");
				ss.clear();
				
				
				// single reads
				ss << filePath << "_barcode_unassigned_2_single"<< toFormatString(m_format);
				TOutputFilter *single4 = new SequenceOutputFilter<TString, TIDString>(ss.str(), m_format);
				
				ss.str("");
				ss.clear();
				
				
				filters f2;
				f2.f1 = filter3;
				f2.f2 = filter4;
				f2.single1 = single3;
				f2.single2 = single4;
				
				m_outputMap[0] = f2;
				break;
			}
			
			case PAIRED:{  //if no barcodes have been specified write paired output
				
				m_mapsize = 1;
				m_outputMap = new filters[m_mapsize];
				stringstream ss;
				
				ss << filePath << "_1"<< toFormatString(m_format);
				TOutputFilter *filter5 = new TOutputFilter(ss.str(),m_format);
				
				ss.str("");
				ss.clear();
				
				// single
				ss << filePath << "_1_single" << toFormatString(m_format);
				TOutputFilter *single5 = new TOutputFilter(ss.str(), m_format);
				
				ss.str("");
				ss.clear();
				
				ss << shortFilename << "_1"<< toFormatString(m_format);
				
				ss.str("");
				ss.clear();

				ss << filePath << "_2"<< toFormatString(m_format);
				TOutputFilter *filter6 = new TOutputFilter(ss.str(), m_format);
				
				ss.str("");
				ss.clear();
				
				
				// single
				ss << filePath << "_2_single"<< toFormatString(m_format);
				TOutputFilter *single6 = new TOutputFilter(ss.str(), m_format);
				
				ss.str("");
				ss.clear();
				
				filters f3;
				f3.f1 = filter5;
				f3.f2 = filter6;
				f3.single1 = single5;
				f3.single2 = single6;
				
				m_outputMap[0] = f3;
				break;
			}
			
			case SINGLE:{
				
				m_mapsize = 1;
				m_outputMap = new filters[m_mapsize];
				stringstream ss;
				
				ss << filePath << toFormatString(m_format);
				TOutputFilter *filter7 = new SequenceOutputFilter<TString,TIDString>(ss.str(), m_format);
				
				// ss.str("");
				// ss.clear();
				// ss << shortFilename;
				
				filters f4;
				f4.f1 = filter7;
				f4.f2 = NULL;
				
				m_outputMap[0] = f4;
				break;
			}
			
			case SINGLE_BARCODED:
			case SINGLE_BARCODED_WITHIN_READ:{
				
				m_mapsize = m_mapsize = barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				stringstream ss;
				
				for(unsigned int i=0; i < barcodes->size(); ++i){
					
					ss << filePath << "_barcode_" << barcodes->at(i).first->getSequenceTag() << toFormatString(m_format);
					TOutputFilter *filter8 = new TOutputFilter(ss.str(), m_format);
					
					ss.str("");
					ss.clear();
					
					filters f5;
					f5.f1 = filter8;
					f5.f2 = NULL;
					f5.single1 = NULL;
					f5.single2 = NULL;

					m_outputMap[i + 1] = f5;
				}
				
				ss << filePath << "_barcode_unassigned"<< toFormatString(m_format);
				TOutputFilter *filter9 = new SequenceOutputFilter<TString, TIDString>(ss.str(), m_format);
				
				filters f6;
				f6.f1 = filter9;
				f6.f2 = NULL;
				
				m_outputMap[0] = f6;
				break;
			}
		}
	}
	
	
	virtual ~MultiplexedOutputFilter(){};
	
	
	void setMinReadlength(int minLength){
		m_minLength = minLength;
	}
	
	
	unsigned long getNrGoodReads(){
		unsigned long good = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			good += m_outputMap[i].f1->getNrGoodReads();
			if(m_outputMap[i].f2 != NULL) good += m_outputMap[i].f2->getNrGoodReads();
		}
		return good;
	}
	
	
	void* operator()(void* item) {
		
		using namespace flexbar;
		
		MultiplexedRead<TString, TIDString> *myRead = static_cast< MultiplexedRead<TString, TIDString>* >(item);
		bool l1ok = false, l2ok = false;
		
		switch(m_runType){
			case PAIRED_BARCODED:
			case PAIRED:{
				
				if(myRead->m_r1 != NULL && myRead->m_r2 != NULL){
					
					// now check if both reads have minLength
					if(length(myRead->m_r1->getSequence()) >= m_minLength) l1ok = true;
					if(length(myRead->m_r2->getSequence()) >= m_minLength) l2ok = true;
					
					if(l1ok && l2ok){
						m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
						m_outputMap[myRead->m_barcode_id].f2->writeRead(myRead->m_r2);
					}
					else if(l1ok && ! l2ok){
						m_outputMap[myRead->m_barcode_id].single1->writeRead(myRead->m_r1);
					}
					else if(! l1ok && l2ok){
						m_outputMap[myRead->m_barcode_id].single2->writeRead(myRead->m_r2);
					}
					
					if(! l1ok) m_outputMap[myRead->m_barcode_id].m_cnt_short_1 += 1;
					if(! l2ok) m_outputMap[myRead->m_barcode_id].m_cnt_short_2 += 1;
					
				} break;
			}
			
			case SINGLE_BARCODED:
			case SINGLE_BARCODED_WITHIN_READ:
			case SINGLE:{
				if(myRead->m_r1 != NULL){
					if(length(myRead->m_r1->getSequence()) >= m_minLength){
						
						m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
					}
					else m_outputMap[myRead->m_barcode_id].m_cnt_short_1 += 1;
				}
			}
		}
		
		delete myRead;
		return NULL;
	}
	
	
	void writeLengthDist(){
		for(unsigned int i = 0; i < m_mapsize; i++){
			m_outputMap[i].f1->writeLengthDist();
			if(m_outputMap[i].f2 != NULL) m_outputMap[i].f2->writeLengthDist();
		}
	}
	
	
	unsigned long getNrShortReads(){
		unsigned long nShortReads = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			nShortReads += m_outputMap[i].m_cnt_short_1;
			
			if((m_runType == flexbar::PAIRED_BARCODED) || (m_runType == flexbar::PAIRED)){
				nShortReads += m_outputMap[i].m_cnt_short_2;
			}
		}
		return nShortReads;
	}
	
	
	void printAdapterRemovalStats(){
		using namespace std;
		
		cout << endl << "Adapter removal statistics" << endl;
		cout <<         "==========================" << endl;
		
		const unsigned int maxSpaceLen = 18;
		
		cout << "Adapter:" << string(maxSpaceLen - 8, ' ') << "Removed:" << endl;
		
		for(unsigned int i = 0; i < m_adapters->size(); i++){
			seqan::CharString seqTag = m_adapters->at(i).first->getSequenceTag();
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			cout << seqTag << whiteSpace << m_adapters->at(i).second << endl;
		}
		cout << endl;
	}
	
	
	void printFileSummary(){
		using namespace std;
		
		cout << "Output file statistics" << endl;
		cout << "==========================" << endl;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			cout << "Output read file:            " << m_outputMap[i].f1->getFileName()    << endl;
			cout << "Discarded short reads:       " << m_outputMap[i].m_cnt_short_1        << endl;
			cout << "Reads written to file:       " << m_outputMap[i].f1->getNrGoodReads() << endl;
			
			if(m_runType == flexbar::PAIRED_BARCODED || m_runType == flexbar::PAIRED){
				cout << "Output read file2:           " << m_outputMap[i].f2->getFileName()    << endl;
				cout << "Discarded short reads:       " << m_outputMap[i].m_cnt_short_2        << endl;
				cout << "Reads written to file:       " << m_outputMap[i].f2->getNrGoodReads() << endl;
			}
		}
	}
	
};

#endif /* FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_ */
