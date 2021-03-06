// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_QGRAM_FIND_H
#define SEQAN_HEADER_INDEX_QGRAM_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// QGram finders

	struct FinderQGramLookup_; //Finder that simply looks up the q-gram in the hash table

/**
.Tag.Index Find Algorithm
..tag.QGramFindLookup:q-gram search.
Finds q-grams in a @Spec.IndexQGram@ index using the hash table.
..include:seqan/index.h
*/

	typedef Tag<FinderQGramLookup_> const QGramFindLookup;

//____________________________________________________________________________


	template < typename TText, typename TShapeSpec, typename TSpec >
	struct DefaultFinder<Index<TText, IndexQGram<TShapeSpec, TSpec> > > {
        typedef QGramFindLookup Type;
    };


//////////////////////////////////////////////////////////////////////////////
// _findFirstIndex implementation

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void _findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		QGramFindLookup const)
	{
		typedef Index<TText, TSpec>									TIndex;
		typedef typename Fibre<TIndex, QGramSA>::Type				TSA;
		typedef typename Fibre<TIndex, QGramShape>::Type			TShape;
		typedef typename Fibre<TIndex, QGramDir>::Type				TDir;
		typedef typename Iterator<TSA const, Standard>::Type		TSAIterator;
		typedef typename Iterator<TPattern const, Standard>::Type	TPatternIterator;

		TIndex &index = haystack(finder);
		indexRequire(index, QGramSADir());

		TSAIterator saIt = begin(indexSA(index), Standard());
		TPatternIterator pIt = begin(pattern, Standard());
		TDir const &dir = indexDir(index);
		TShape &shape = indexShape(index);

		finder.range.i1 = saIt + dir[getBucket(index.bucketMap, hash(shape, pIt, length(pattern)))];
		finder.range.i2 = saIt + dir[getBucket(index.bucketMap, hashUpper(shape, pIt, length(pattern)))];
	}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_

