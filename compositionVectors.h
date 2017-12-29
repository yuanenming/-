#ifndef COMPOSITIONVECTORS_H_
#define COMPOSITIONVECTORS_H_

std::map<std::string,double> compositionVector_getScoreVectors(std::map<int, std::map<std::string,int> > KMERS_MAP, int SeqLength, int K);
double compositionVector_getDistance( std::map<std::string,double> seq1, std::map<std::string,double> seq2 );
double ** compositionVector_getDistanceMatrix( std::vector<std::map<std::string,double> >Vectors, int size );
double ** compositionVector( std::map<int, std::map<int, std::map<std::string,int> > > Map, int seqL[], int K);

#endif 


