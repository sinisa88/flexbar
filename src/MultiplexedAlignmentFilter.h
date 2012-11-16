/*
 * MultiplexedAlignmentFilter.h
 *
 *  Created on: Jun 1, 2011
 *      Author: mat
 *
 *  Revised on: Aug 3, 2012
 *      Author: jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDALIGNMENTFILTER_H_
#define FLEXBAR_MULTIPLEXEDALIGNMENTFILTER_H_

#include <sstream>
#include <string>

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include "Enums.h"
#include "MultiplexedRead.h"
#include "AlignmentFilter.h"
#include "AlignmentAlgorithm.h"
#include "AdapterLoader.h"


/* This class processes a MultiplexedRead. It assigns a barcode to the read (if processing barcoded run)
and removes adapter sequences. */

template <class TString, class TIDString, class TAlgorithm >
class MultiplexedAlignmentFilter : public tbb::filter {

private:
	
	flexbar::LogLevel m_verb;
	flexbar::BarcodeDetect m_barType;
	flexbar::AdapterRemoval m_adapRem;
	
	tbb::atomic<flexbar::TrimEnd> m_aEnd;
	tbb::atomic<flexbar::TrimEnd> m_bEnd;
	
	tbb::concurrent_vector<TAdapter> *m_adapters;
	tbb::concurrent_vector<TAdapter> *m_barcodes;
	
	typedef AlignmentFilter<TString, TIDString, AlignmentAlgorithm<TString> > AliFilter;
	AliFilter *m_afilter, *m_bfilter;
	
public:
	
	MultiplexedAlignmentFilter(tbb::concurrent_vector<TAdapter> *adapters, tbb::concurrent_vector<TAdapter> *barcodes, float aThresh, float bThresh, int a_minOverlap, int b_minOverlap, flexbar::TrimEnd bEnd, int match, int mismatch, int gapCost, int b_match, int b_mismatch, int b_gapCost, flexbar::TrimEnd aEnd, flexbar::BarcodeDetect barType, flexbar::AdapterRemoval adapRem, bool remTag, flexbar::LogLevel logLevel, flexbar::FileFormat format) : tbb::filter(parallel){
		
		m_aEnd     = aEnd;
		m_bEnd     = bEnd;
		m_barType  = barType;
		m_adapRem  = adapRem;
		m_verb     = logLevel;
		m_adapters = adapters;
		m_barcodes = barcodes;
		
		m_afilter = new AliFilter(m_adapters, match, mismatch, gapCost, remTag, m_aEnd, m_verb, format);
		m_bfilter = new AliFilter(m_barcodes, b_match, b_mismatch, b_gapCost, remTag, m_bEnd, m_verb, format);
		
		m_afilter->setThreshold(aThresh);
		m_afilter->setMinOverlap(a_minOverlap, (m_adapRem == flexbar::ADAPTIVE));
		
		m_bfilter->setThreshold(bThresh);
		m_bfilter->setMinOverlap(b_minOverlap, false);
		
		if(m_verb == flexbar::TAB)
		std::cout << "read-ID\tadapter-ID\tadapter-start\tadapter-end\toverlap-length\tmismatches\tindels\tallowed-errors" << std::endl;
	}
	
	// BarcodeFilter:  &myRead = *static_cast< SequencingRead<TString,TIDString>* >(myMultiRead->m_b);
	
	
	virtual ~MultiplexedAlignmentFilter(){};
	
	
	void* operator()(void* item){
		
		using namespace flexbar;
		
		if(item != NULL){
			MultiplexedRead<TString, TIDString> *myRead = static_cast< MultiplexedRead<TString, TIDString>* >(item);
			
			// barcode detection
			if(m_barType != BOFF){
				switch(m_barType){
					case BARCODE_READ:        myRead->m_barcode_id = m_bfilter->align(myRead->m_b,  false); break;
					case WITHIN_READ:         myRead->m_barcode_id = m_bfilter->align(myRead->m_r1, false); break;
					case REMOVE_WITHIN_READ:  myRead->m_barcode_id = m_bfilter->align(myRead->m_r1, true);  break;
					case BOFF:                                                                              break;
				}
			}
			
			// adapter removal
			if(m_adapRem != AOFF){
				m_afilter->align(myRead->m_r1, true);
				
				if(myRead->m_r2 != NULL)
				m_afilter->align(myRead->m_r2, true);
			}
			return myRead;
		}
		else return NULL;
	}
	
	
	int getNrModifiedReads(){
		return m_afilter->getNrModifiedReads();
	}
	
	
	void printAdapterOverlapStats(){
		if(m_afilter->getNrModifiedReads() > 0)
		std::cout << "Min, max, mean and median adapter overlap length: "<< m_afilter->getMinOverlapLength()  << " / "<< m_afilter->getMaxOverlapLength() << " / " << m_afilter->getMeanOverlapLength() << " / " << m_afilter->getMedianOverlapLength() << std::endl << std::endl;
	}
	
};

#endif /* FLEXBAR_MULTIPLEXEDALIGNMENTFILTER_H_ */
