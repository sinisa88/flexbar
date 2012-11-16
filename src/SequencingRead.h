/*
 * Read.h
 *
 *  Created on: Apr 12, 2010
 *      Author: mat
 *
 *  Revised on: Jul 16, 2012
 *      Author: jtr
 */

#ifndef FLEXBAR_SEQUENCINGREAD_H_
#define FLEXBAR_SEQUENCINGREAD_H_

#include <seqan/basic.h>


/* A Sequencing read consists of a nucleotide sequence (color or basepair space),
a sequence name and optionally a quality string plus the quality scaling. */

template <class TString, class TIDString>
class SequencingRead {

private:
	
	TString m_seq;
	TIDString m_tag, m_qual;
	
public:
	
	SequencingRead(){
		m_tag = "";
		m_seq = "";
	}
	
	
	SequencingRead(TString source, TIDString sequence_tag){
		
		m_tag = sequence_tag;
		m_seq = source;
	}
	
	
	SequencingRead(TString source, TIDString sequence_tag, TIDString qual){
		
		m_tag  = sequence_tag;
		m_seq  = source;
		m_qual = qual;
	}
	
	
	void setSequenceTag(TString tag){
		m_tag = tag;
	}
	
	void setSequence(TString seq){
		m_seq = seq;
	}
	
	void setQuality(TString qual){
		m_qual = qual;
	}
	
	
	TIDString getSequenceTag(){
		return m_tag;
	}
	
	TString getSequence(){
		return m_seq;
	}
	
	TIDString getQuality(){
		return m_qual;
	}
	
	
	virtual ~SequencingRead(){};
	
};

#endif /* FLEXBAR_SEQUENCINGREAD_H_ */
