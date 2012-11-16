/*
 *  MultiplexedInputFilter.h
 *
 *  Created on: May 31, 2011
 *      Author: mat
 *
 *  Revised on: Jul 16, 2012
 *      Author: jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDINPUTFILTER_H_
#define FLEXBAR_MULTIPLEXEDINPUTFILTER_H_

#include <tbb/pipeline.h>

#include "MultiplexedRead.h"
#include "SequenceInputFilter.h"


/* This class will handle up to 3 file sources at the same time (paired read input plus barcode reads)
   and create a MultiplexedRead depending on the run type (single-end, paired-end and/or barcoded). */

template <typename TString, typename TIDString>
class MultiplexedInputFilter : public tbb::filter {

private:
	
	bool m_isPaired, m_useBarcodeRead;
	int m_maxUncalled, m_cutLen_begin, m_cutLen_end, m_phredPreQual, m_minReadLen;
	long m_uncalled, m_uncalledPairs;
	
	std::string m_filePath;
	flexbar::FileFormat m_format;
	flexbar::QualityType m_qual;
	
	SequenceInputFilter<TString,TString > *m_f1;
	SequenceInputFilter<TString,TString > *m_f2;
	SequenceInputFilter<TString,TString > *m_b;
	
public:
	
	MultiplexedInputFilter(std::string filePath, flexbar::FileFormat format, int maxUncalled, int cutLen_begin, int cutLen_end, int minReadLen, int phredPreQual, flexbar::QualityType qual) : filter(serial_in_order){
		
		m_format   = format;
		m_qual     = qual;
		m_filePath = filePath;
		
		m_cutLen_begin = cutLen_begin;
		m_cutLen_end   = cutLen_end;
		m_phredPreQual = phredPreQual;
		m_minReadLen   = minReadLen;
		m_maxUncalled  = maxUncalled;
		
		m_uncalled       = 0;
		m_uncalledPairs  = 0;
		m_isPaired       = false;
		m_useBarcodeRead = false;
		
		m_f1 = new SequenceInputFilter<TString, TString>(m_filePath, m_format, m_maxUncalled);
		
		m_f1->setPreTrimBegin(cutLen_begin);
		m_f1->setPreTrimEnd(cutLen_end);
		m_f1->setMinReadLength(minReadLen);
		m_f1->setPrePhredTrim(phredPreQual, qual, true);
		
		m_f2 = NULL;
		m_b  = NULL;
	}
	
	
	virtual ~MultiplexedInputFilter(){}
	
	
	void setPairedFile(std::string filePath2){
		
		m_f2 = new SequenceInputFilter<TString, TString>(filePath2, m_format, m_maxUncalled);
		m_f2->setPreTrimBegin(m_cutLen_begin);
		m_f2->setPreTrimEnd(m_cutLen_end);
		m_f2->setMinReadLength(m_minReadLen);
		m_f2->setPrePhredTrim(m_phredPreQual, m_qual, false);
		
		m_isPaired = true;
	}
	
	
	void setBarcodeReadsFile(std::string bfileName){
		
		m_b = new SequenceInputFilter<TString, TString>(bfileName, m_format, m_maxUncalled);
		m_b->setPreTrimBegin(0);
		m_b->setPreTrimEnd(0);
		
		m_useBarcodeRead = true;
	}
	
	
	void* operator()(void*){
		
		using namespace std;
		
		SequencingRead<TString, TIDString> *myRead1 = NULL, *myRead2 = NULL, *myBarcode = NULL;
		MultiplexedRead<TString, TIDString> *mRead = NULL;
		
		while(true){
			
			bool valid = false, valid2 = false, vBR = false, uncalled = true, uncalled2 = true, uBR = true;
			
			// single read input
			if(! m_isPaired){
				
				while(uncalled || ! valid){
					myRead1 = static_cast< SequencingRead<TString, TIDString>* >(m_f1->getRead(valid, uncalled));
					
					if(m_useBarcodeRead) myBarcode = static_cast< SequencingRead<TString, TIDString>* >(m_b->getRead(vBR, uBR));
					
					if(myRead1 == NULL) return NULL;
					
					if(m_useBarcodeRead && myBarcode == NULL){
						cerr << "Error: read without barcode read, or file reading error!" << endl << endl;
						exit(1);
					}
					
					if(uncalled) m_uncalled++;
				}
			}
			
			// paired read input
			else{
				
				// find at least one valid read
				while(uncalled || uncalled2 || (! valid && ! valid2)){
					
					// (!v1 || myRead1 == NULL) || (!v2 || myRead2 == NULL)  to remove outer while(true)
					
					myRead1 = static_cast< SequencingRead<TString, TIDString>* >(m_f1->getRead(valid, uncalled));
					myRead2 = static_cast< SequencingRead<TString, TIDString>* >(m_f2->getRead(valid2, uncalled2));
					
					if(m_useBarcodeRead) myBarcode = static_cast< SequencingRead<TString, TIDString>* >(m_b->getRead(vBR, uBR));
					
					// end of files reached
					if(myRead1 == NULL && myRead2 == NULL) return NULL;
					
					else if(myRead1 == NULL || myRead2 == NULL){
						cerr << "Error: single read in paired mode, or file reading error!" << endl << endl;
						exit(1);
					}
					
					if(m_useBarcodeRead && myBarcode == NULL){
						cerr << "Error: reads without barcode read or file reading error!" << endl << endl;
						exit(1);
					}
					
					if(uncalled)  m_uncalled++;
					if(uncalled2) m_uncalled++;
					
					if(uncalled || uncalled2) m_uncalledPairs++;
				}
				
				if(! valid2) myRead2 = NULL;
				
				else if(! valid){
					myRead1 = myRead2;
					myRead2 = NULL;
				}
			}
			
			mRead = new MultiplexedRead<TString, TIDString>(myRead1, myRead2, myBarcode);
			return mRead;
		}
	}
	
	
	long getNrUncalledReads(){
		return m_uncalled;
	}
	
	
	long getNrUncalledPairedReads(){
		return m_uncalledPairs;
	}
	
	
	long getNrProcessedReads(){
		if(m_isPaired) return m_f1->getNrProcessedReads() + m_f2->getNrProcessedReads();
		else return m_f1->getNrProcessedReads();
	}
	
	
	long getNrLowPhredReads(){
		if(m_isPaired) return m_f1->getNrLowPhredReads() + m_f2->getNrLowPhredReads();
		else           return m_f1->getNrLowPhredReads();
	}
	
	
	long getNrShortPhredReads(){
		if(m_isPaired) return m_f1->getNrShortPhredReads() + m_f2->getNrShortPhredReads();
		else           return m_f1->getNrShortPhredReads();
	}
	
};

#endif /* FLEXBAR_MULTIPLEXEDINPUTFILTER_H_ */
