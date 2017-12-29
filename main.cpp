#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <map>
#include <math.h> 
#include <iomanip>
#include <algorithm>
#include <stdlib.h>

#include "compositionVectors.h"
#include "wordContext.h"
#include "spacedWordFrequency.h"
#include "neighborJoining.h"
#include "baseSub.h"

/* Classes */
#include "GraphClass.h"


using namespace std;


//Define Iterators for each Level in the Map
typedef map<int, map<int, map<string,int> > >::iterator sequence_it;
typedef map<int, map<string,int> >::iterator kmers_it;
typedef map<string,int>::iterator freq_it;
typedef map<string,double>::iterator score_it;

/*
 * All small operational functions after main function.
 * Declare their names.
 */
bool endsWith (std::string const &fullString, std::string const &ending);

/**
 * Function for reading a FASTA file.
 * @arguments: 
 *        string filename: the name of the file to be read
 * @return
 *         A set of sets of sequences as a vector
 */
vector< vector<string> > getSequencesFastaFile( string filename ){
    
    ifstream FASTA_FILE;
    
    FASTA_FILE.open(filename.c_str());
    
    //File error
    if (FASTA_FILE.fail()){
        cerr << "Problem with FASTA file." << endl;
        exit(1);
    }
    
    vector <string>  FASTA_SEQUENCE_SET;
    string  FASTA_SEQUENCE = "";
    string line;
    
    string names="";
    
    //Read File
    while ( getline (FASTA_FILE,line) ){
        
        if (line.empty()) continue;
        
        //Check whether the line contains sequence's chars or sequence's identification
        if ( line[0] == '>' ){//line starts with '>' , save the current in vector and set new sequence
            names.append(line);
            if ( !FASTA_SEQUENCE.empty() ){
                FASTA_SEQUENCE_SET.reserve(1);
                FASTA_SEQUENCE_SET.push_back(FASTA_SEQUENCE);
                FASTA_SEQUENCE = "";
            }
            
            continue;    
        }
        FASTA_SEQUENCE.append( line );
    }
    
    //save the last sequence
    if ( !FASTA_SEQUENCE.empty() ){
        FASTA_SEQUENCE_SET.reserve(1);
        FASTA_SEQUENCE_SET.push_back(FASTA_SEQUENCE);
        FASTA_SEQUENCE = "";
    }
    FASTA_SEQUENCE_SET.reserve(1);
    FASTA_SEQUENCE_SET.push_back(names);
    
    //Close file
    FASTA_FILE.close();
    
    vector< vector<string> > sets;
    sets.push_back(FASTA_SEQUENCE_SET);
    
    return sets;
}

/**
 * Function for reading a PAML file.
 * @arguments: 
 *        string filename: the name of the file to be read
 * @return
 *         a set of sequences as a vector
 *
 * TODO: ITS IMPLEMENTATION... What we have here is the copy of fasta file reading.
 */
vector< vector<string> > getSequencesPamlFile( string filename ){
    
    ifstream PAML_FILE;
    
    PAML_FILE.open(filename.c_str());
    
    //File error
    if (PAML_FILE.fail()){
        cerr << "Problem with PAML file." << endl;
        exit(1);
    }
    
    vector< vector<string> > ALLSets;
    vector <string> PAML_SEQUENCE_SET;
    string line;
    
    string names="";
    
    //Read File
    while ( getline (PAML_FILE, line) ){
        
        if (line.empty()) continue;
        
        std::size_t found = line.find_first_of(" ");
        
        if ( found != 0 ){    
            string sequence_name = line.substr(0,found);
            string sequence = line.substr(found+1, ( line.size()-sequence_name.size() ) );
            //Remove all spaces from the sequence
            char space = ' ';
            sequence.erase ( remove(sequence.begin(), sequence.end(), space), sequence.end() );
            
            PAML_SEQUENCE_SET.reserve(1);
            PAML_SEQUENCE_SET.push_back(sequence);
            names.append(">");
            names.append(sequence_name);
            
        }else if ( !names.empty() ){
            //Add the sequence names at the last position
            PAML_SEQUENCE_SET.reserve(1);
            PAML_SEQUENCE_SET.push_back(names);
            
            //Add the sequences set at the ALL Sets
            ALLSets.reserve(1);
            ALLSets.push_back(PAML_SEQUENCE_SET);
            names = "";
            PAML_SEQUENCE_SET.clear();
        }        
    }
    //Close file
    PAML_FILE.close();
    
    //Add the last sequence set!
    PAML_SEQUENCE_SET.reserve(1);
    PAML_SEQUENCE_SET.push_back(names);
    
    //Add the sequences set at the ALL Sets
    ALLSets.reserve(1);
    ALLSets.push_back(PAML_SEQUENCE_SET);
       
    return ALLSets;
}

/**
 * Generic function for reading a file.
 * @arguments: 
 *        string filename: the name of the file to be read
 * @return
 *         a set of sequences as a vector
 */
vector< vector<string> > getSequencesFromFile( string filename ){
    
    if ( endsWith(filename, ".fasta") ) //check whether the file has extension .fasta 
        return getSequencesFastaFile(filename);
    /**
     * TODO: Read Paml file. Bakker hasn't reply yet.
     * How to read it!
     */
    if ( endsWith(filename, ".paml") )
        return getSequencesPamlFile(filename);
    
    cout << "No acceptable file given!";
    exit(1);
}

vector<string> getSequencesNames(string SequencesNames){
    vector<string> NameSet;
    std::string delimiter = ">";
    size_t pos = 0;
    std::string token;
    while ((pos = SequencesNames.find(delimiter)) != std::string::npos) {
        token = SequencesNames.substr(0, pos);
        SequencesNames.erase(0, pos + delimiter.length());
        if (token.empty()){continue;} 
        NameSet.reserve(1);
        NameSet.push_back(token);
    }
    NameSet.reserve(1);
    NameSet.push_back(SequencesNames);//push back the last name
    return NameSet;
}

std::map<int,std::map<std::string,int> > getKMers( string Sequence ){
    
    map<int, map<string,int> > Maps;
    
    int k = 1;
    int kmer;
    while ( k < 11 ){ // up to 10-mer
        kmer = k;
        map<string,int> kmer_map;
        for ( int i = k; i < Sequence.size() - k + 1; i++ ){
            if( kmer_map[Sequence.substr(i-k,k)]   ){//non zero -- increase
                kmer_map[Sequence.substr(i-k,k)] += 1;
            }else{//zero -- add
                kmer_map[Sequence.substr(i-k,k)] = 1;
            }
        }
        Maps[k] = kmer_map;
        k++;
    }
    return Maps;    
}

/*
 * Function for correction NaN values
 *
 */
double** MatrixCorrection(double** DistanceMatrix, int size){
        
    double** Matrix = new double*[size];
    for ( int i = 0; i < size; i++ ){
        Matrix[i] = new double[size];
        for ( int j = 0; j < size; j++ ){
            Matrix[i][j] = DistanceMatrix[i+1][j+1];
        }
    }
    return Matrix;        
}

int* getSequencesLength(vector<string> Set){
    int * seqL = new int[Set.size()];
    int j = 0;
    for( vector<string>::iterator i = Set.begin(); i != Set.end(); ++i ){
        seqL[j] = (*i).length();
        j++;
    }
    return seqL;
}

double** normalization(double** Matrix, int size){
    
    //MIN VALUE IS ALWAYS ZERO BECAUSE OF DIAGONAL
    double MAX = 0;
    
    for ( int i = 0; i < size; i++ ){
        for ( int j = i+1; j < size; j++ ){
            if ( MAX < Matrix[i][j]  ){
                MAX = Matrix[i][j];
            }
        }
    }
    
    for ( int i = 0; i < size; i++ ){
        for ( int j = i+1; j < size; j++ ){
            Matrix[i][j] = Matrix[i][j] / MAX ;
            Matrix[j][i] = Matrix[i][j] ;
        }
    }
    
    return Matrix;
}

void writeToFile(string filename, double ** Matrix, int SEQUENCES_SIZE){
    ofstream of(filename.c_str());
    for ( int i =0; i< SEQUENCES_SIZE; i++){
        for ( int j =0; j< SEQUENCES_SIZE; j++){
            of << setw(12) <<Matrix[i][j];
        }
        of << endl;
    }
    of.close();
    cout << "Results on: "<< filename << " file." << endl;    
    
}

void calculatePhyloTrees( string filename, int option, int kmer, bool writeMatrices ){
    
    double** Matrix;
    string f_out_name;
    
    cout << "Choose the name for the file of Netwick Tree Format: ";
    string Netwick = "";
    cin >> Netwick;
    ofstream ov(Netwick.c_str());

    /**
     * Read a (fasta or paml) file as an input.
     * Get all sequences.
     * Return k-mers for each sequence.
     * Each k-mer is an hash map which contains as key the k-mer (e.g. 'AT') and as value its frequency.
     */
    
    vector< vector< string > > Sets = getSequencesFromFile(filename);
    
    int iter = 1;
    for( std::vector< std::vector <std::string> >::iterator Set = Sets.begin(); Set != Sets.end(); ++Set ){
        system("cls");
        cout << "Sequence Set " << iter <<" /" << Sets.size() << endl << endl;
        string SequencesNames  = (*Set).back();//get the last element of the Sequence Set
        (*Set).pop_back();//Remove it! It's the Sequences' names!
        vector<string> SequencesNamesV = getSequencesNames(SequencesNames);
        
        int SEQUENCES_SIZE = (*Set).size();
                
        /*
         * BEFORE READING ALL SEQUENCES WE SHOULD 
         * PREPROCESS THEM FOR 'N' OCCURENCES.
         */
        (*Set) = cleanData((*Set));
        
        map<int, map<int, map<string,int> > > kmers_freq;
        
        /*
         * READ SEQUENCES AND PARSE THEM
         */
        int j = 1;
        cout << "Read sequences.."<<endl;
        for( vector<string>::iterator i = (*Set).begin(); i != (*Set).end(); ++i ){
                cout << "Sequence " << j << " / " << SEQUENCES_SIZE<<endl;
                kmers_freq[j] = getKMers( (*i) );
                j++;
//                system("cls");
        }
        int* seqL = getSequencesLength( (*Set) );
        
        switch( option ){
            
            case 1: {
                Matrix = compositionVector( kmers_freq , seqL, kmer);
                string s = static_cast<ostringstream*>( &(ostringstream() << kmer) )->str();
                f_out_name = "CompositionVectors.";
                f_out_name.append(s);
                f_out_name.append("mer.Results.txt");
                break;
            }
            case 2:{
                Matrix = wordContext(kmers_freq);    
                f_out_name = ("WordContext.txt");
                Matrix = MatrixCorrection(Matrix, SEQUENCES_SIZE);
                break;
            }
            case 3:{
                Matrix = spacedWordFrequency(kmers_freq);    
                f_out_name = ("SpacedWordFrequency.txt");
                Matrix = MatrixCorrection(Matrix, SEQUENCES_SIZE);
                break;
            }
            default:
                break;
        }
        
        /* If true write the DistanceMatrix on file */
        if ( writeMatrices ){
            ofstream of(f_out_name.c_str());
            for ( int i =0; i< SEQUENCES_SIZE; i++){
                for ( int j =0; j< SEQUENCES_SIZE; j++){
                    of << setw(12) <<Matrix[i][j];
                }
                of << endl;
            }
            of.close();
            cout << "Results on: "<< f_out_name << " file." << endl;    
        } 
        /* NORMALIZE MATRIX TO [0, 1] RANGE */
        Matrix = normalization(Matrix, SEQUENCES_SIZE);
        
        /*GET THE NEIGHBOR JOINING FORM*/
        Graph tree = neighborJoining(Matrix, SEQUENCES_SIZE, SequencesNamesV);
//        tree.print();
        ov << tree.getNetwickTreeFormat() << ";" << endl;
        iter++;
    }
    ov.close();
    cout << "Netwick Tree Format on: "<< Netwick << endl;
}


/*
 * TODO: Add Argument list for giving the user the option to add a filename
 */

int main( int argc, char* argv[] ){
    string filename;    
    if ( argc < 2){
        cout << "Filename is missing" << endl;
        exit(0);
    }
    filename = argv[1];
    
    if ( !endsWith(filename, ".fasta") && !endsWith(filename, ".paml") ){
        cout << "File extension is wrong" << endl;
        exit(0);
        
    }
    
    
    cout << "Greetings my lord!" << endl;
    cout << "Choose you destiny.." << endl << endl;
    
    cout << "\t1 - Composition Vectors" << endl;
    cout << "\t2 - Word Context" << endl;
    cout << "\t3 - Spaced Word Frequency" << endl;
    
    int option, k;
    cout << "Option: ";
    cin >> option;    
    
    
    while( option != 1 && option != 2 && option != 3 ){
        cout << "You fool.\nChoose wiser.." << endl;
        cout << "Option: ";
        cin >> option;
    }
        
    /* K-mers choice */
    cout << "Now, choose the number of k-mers you want. But beware! It's only 3 to 10\nk-mers: ";
    cin >> k;
    while( k < 3 || k >10){
        cout << "I feel like you like games. Please, choose wiser: ";
        cin >> k;
    }        
    
    /* Start calculating. . . */
    calculatePhyloTrees(filename, option, k, false);
    system("Rscript R_draw_tree.R");
    
    cout << "\n\nGoodbye my lord.." << endl;
    
    return 0;
}


/*
 * Function endsWith
 * @attr: 
 *        fullString: the string we want to check
 *        ending: the string with which we check the fullString
 * @return: 
 *         true if fullString ends with ending
 *         false otherwise
 */
bool endsWith (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}


