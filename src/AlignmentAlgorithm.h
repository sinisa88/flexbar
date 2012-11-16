/*
 * AlignmentAlgorithm.h
 *
 *  Created on: Jul 15, 2010
 *      Author: mat
 */

#ifndef FLEXBAR_ALIGNMENTALGORITHM_H_
#define FLEXBAR_ALIGNMENTALGORITHM_H_

#include <iostream>
#include <string>

#include <seqan/basic.h>
#include <seqan/align.h>


template <typename TString>
class AlignmentAlgorithm {

private:
	
	seqan::Score<int> *m_score;
	
public:
	
	AlignmentAlgorithm(int match, int mismatch, int gapCost){
		m_score     = new seqan::Score<int>(match, mismatch, gapCost);
	};
	
	
	virtual ~AlignmentAlgorithm(){
		delete m_score;
	};
	
	
	void align(TString &querySeq, TString &read, int &gaps, int &mismatches, int &startPos, int &endPos, int &startPosA, int &endPosA, int &startPosS, int &endPosS, int &aliScore, std::stringstream &aliString, flexbar::TrimEnd trimEnd){
		
		using namespace std;
		using namespace seqan;
		
		typedef Align<CharString, ArrayGaps> TAlign;
		typedef Row<TAlign>::Type TRow;
		typedef Iterator<TRow>::Type TRowIterator;
		
		TAlign align;
		resize(rows(align), 2);
	    assignSource(row(align, 0), read);
	    assignSource(row(align, 1), querySeq);
		
		
		if(trimEnd == flexbar::RIGHT || trimEnd == flexbar::RIGHT_TAIL){
			AlignConfig<true, false, true, true> ac;
			
			// aliScore = globalAlignment(align, *m_score, ac, 2, 2, NeedlemanWunsch());
			
			aliScore = globalAlignment(align, *m_score, ac, NeedlemanWunsch());  // MyersHirschberg();
		}
		else if(trimEnd == flexbar::LEFT  || trimEnd == flexbar::LEFT_TAIL){
			AlignConfig<true, true, false, true> ac;
			
			aliScore = globalAlignment(align, *m_score, ac, NeedlemanWunsch());
		}
		else{
			AlignConfig<true, true, true, true> ac;
			
			aliScore = globalAlignment(align, *m_score, ac, NeedlemanWunsch());
		}
		
		
		TRow &row1 = row(align,0);
		TRow &row2 = row(align,1);
		
		startPosS = toViewPosition(row1, 0);
		startPosA = toViewPosition(row2, 0);
		endPosS   = toViewPosition(row1, length(source(row1)));
		endPosA   = toViewPosition(row2, length(source(row2)));
		
		// calculate overlap start and end
		if(startPosA > startPosS) startPos = startPosA;
		else                      startPos = startPosS;
		
		if(endPosA > endPosS) endPos = endPosS;
		else                  endPos = endPosA;
		
		aliString << align;
		
		
		// compute number of mismatches and gaps
		TRowIterator it1 = begin(row1);
		TRowIterator it2 = begin(row2);
		
		int aliPos = 0;
		gaps       = 0;
		mismatches = 0;
		
		for(; it1 != end(row1); ++it1){
			if(startPos <= aliPos && aliPos < endPos){
				
				if(isGap(it1) || isGap(it2)) gaps++;
				else if(*it1 != *it2)        mismatches++;
			}
			aliPos++;
			it2++;
		}
		
		// cout << endl << endl << gaps << endl<< mismatches << endl << align;
	}
	
};

#endif /* FLEXBAR_ALIGNMENTALGORITHM_H_ */
