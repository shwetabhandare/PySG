// $Id: aln.h 296 2011-05-11 07:24:40Z satoken $

#ifndef __INC_ALN_H__
#define __INC_ALN_H__

#include <boost/version.hpp>
#if BOOST_VERSION >= 103800
#include <boost/spirit/include/classic.hpp>
#else
#include <boost/spirit.hpp>
#endif
#include <list>
#include "../common/rna.h"

#ifndef BOOST_SPIRIT_CLASSIC_NS
#define BOOST_SPIRIT_CLASSIC_NS boost::spirit
#endif

template < class Seq >
bool
load_aln(MASequence<Seq>& ma, BOOST_SPIRIT_CLASSIC_NS::file_iterator<>& fi);

template < class Seq >
bool
load_aln(std::list<Seq>& ma, BOOST_SPIRIT_CLASSIC_NS::file_iterator<>& fi);

template < class Seq >
bool
load_aln(std::list< MASequence<Seq> >& ma, const char* filename);

#endif	// __INC_ALN_H__

// Local Variables:
// mode: C++
// End: