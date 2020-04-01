#include <iostream>
using std::boolalpha;
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <sstream>
using std::ostringstream;
//yan cong
/*
This project concerns bioinformatics (specifically genetics).
We need to write functions involved with turning a DNA sequence into 
an amino acid sequence (protein).
*/

/*
Citation:

I watached this video https://www.youtube.com/watch?v=oefAI2x2CQM;

I visited these links:
 https://en.cppreference.com/w/cpp;
and http://www.cplusplus.com/;

I checked piazza posts frequently.
*/

// This is the header file for "bio.cpp"
// it has declarations and explanations for the implementation file.
#include "bio.h"

/*
This function returns true if and only if
every character in the input is one of ATCG.
*/
bool IsValidDNASequence(const std::string &input)
{
    // assign input to a string variable, so that it can be modified
    string temp_ipt = input;
    // define the string variable "ATCG"
    string atcg = "ATCG";
    // create a flage, initially true
    bool flag = true;
    // input size
    int iptsize = temp_ipt.size();
    // loop over each character in temp_ipt
    for (int i = 0; i < iptsize; ++i)
    {
        string::size_type pos = atcg.find(temp_ipt[i]);
        // the flag becomes flase if the char is not found in atcg
        if (pos == string::npos)
        {
            flag = false;
        }
    }
    return flag; // return a bool value
}

/*
This function calculates the reverse complement DNA sequence.

Takes: 
the first argument is the sequence, 
the second argument is a pointer to an empty string, 
which is to be modified to store the result.

This is obtained by reversing the input sequence and swaping each
nucleotide/letter with it's complement:
A <-> T
C <-> G

Example:
input = AAATTCGGGG
reverse = GGGGCTTAAA
reverse complement = CCCCGAATTT
*/

void GetReverseComplementSequence(const std::string &input, std::string *const output)
{
    // assign input to a string variable, so that it can be modified
    string temp_ipt = input;
    // define a string variable (for now, an empty string)
    string temp_opt;
    // the size of temp_ipt
    int iptsize = temp_ipt.size();
    // loop over each char in temp_ipt
    for (int i = 0; i < iptsize; ++i)
    {
        // append the last element in temp_ipt to temp_opt first
        // then apend the second to last, etc.
        // this will give the reverse
        temp_opt.push_back(temp_ipt[iptsize - i - 1]);
    }
    // the size of the temp_opt
    int optsize = temp_opt.size();
    // loop over each char in temp_opt
    for (int j = 0; j < (optsize); ++j)
    {
        char C = temp_opt[j];
        switch (C)
        { // Classify each character
        case 'A':
            temp_opt[j] = 'T';
            break;
        case 'T':
            temp_opt[j] = 'A';
            break;
        case 'C':
            temp_opt[j] = 'G';
            break;
        case 'G':
            temp_opt[j] = 'C';
            break;
        default:   // Other character
            break; // it's always safe to have a break here
        }
    }
    // this will give the complement
    // deferencing output
    *output = temp_opt;
}

/*
This function returns the RNA transcript from a DNA sequence.
A RNA transcript is the reverse complement of the DNA sequence, 
but RNA has U (uracil) instead of T (thiamine).

Example:
ASSERT_EQ(GetRNATranscript("AAATTCGGGG"), "CCCCGAAUUU");
*/
std::string GetRNATranscript(const std::string &input)
{
    // assign input to a string variable, so that it can be modified
    string temp_ipt = input;
    // define a string variable (for now, an empty string)
    string temp_opt;
    // call GetReverseComplementSequence function
    GetReverseComplementSequence(temp_ipt, &temp_opt);
    // everytime there is a 'T', replace it with 'U'
    int iptsize = temp_ipt.size();
    for (int i = 0; i < (iptsize); ++i)
    {
        if (temp_opt[i] == 'T')
        {
            temp_opt[i] = 'U';
        }
    }
    return temp_opt;
}

/*
This function returns a vector of vector of strings
with each possible RNA reading frame from the given DNA sequence.
*/
std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string &input)
{
    // assign input to a string variable, so that it can be modified
    string temp_ipt = input;
    // define original RNA transcript
    string original_trs = temp_ipt;
    // define antiparallel RNA transcript (for now, an empty stirng)
    string antip_trs;
    // the size of the temp_ipt
    int iptsize = temp_ipt.size();
    // make sure 'T' is changed to 'U' in the original RNA transcript
    for (int i = 0; i < iptsize; ++i)
    {
        if (original_trs[i] == 'T')
        {
            original_trs[i] = 'U';
        }
    }
    // call function GetRNATranscript
    // assign the result to antiparallel RNA transcript
    antip_trs = GetRNATranscript(temp_ipt);
    // The offsets (starting at pos 0) of the RNA transcripts
    int original_s1 = original_trs.size() / 3;
    // create a vector of string
    vector<string> vec1(original_s1, "abc");
    // fill in the vector
    for (int i = 0; i < original_s1; ++i)
    {
        vec1[i][0] = original_trs[3 * i];
        vec1[i][1] = original_trs[3 * i + 1];
        vec1[i][2] = original_trs[3 * i + 2];
    }
    // The offsets (starting at pos 1) of the RNA transcripts
    int original_s2 = (original_trs.size() - 1) / 3;
    // create a vector of string
    vector<string> vec2(original_s2, "abc");
    // fill in the vector
    for (int i = 0; i < original_s2; ++i)
    {
        vec2[i][0] = original_trs[3 * i + 1];
        vec2[i][1] = original_trs[3 * i + 2];
        vec2[i][2] = original_trs[3 * i + 3];
    }
    // The offsets (starting at pos 2) of the RNA transcripts
    int original_s3 = (original_trs.size() - 2) / 3;
    // create a vector of string
    vector<string> vec3(original_s3, "abc");
    // fill in the vector
    for (int i = 0; i < original_s3; ++i)
    {
        vec3[i][0] = original_trs[3 * i + 2];
        vec3[i][1] = original_trs[3 * i + 3];
        vec3[i][2] = original_trs[3 * i + 4];
    }
    // The offsets (starting at pos 0) of the RNA transcripts
    int antip_s1 = antip_trs.size() / 3;
    // create a vector of string
    vector<string> vec4(antip_s1, "abc");
    // fill in the vector
    for (int i = 0; i < antip_s1; ++i)
    {
        vec4[i][0] = antip_trs[3 * i];
        vec4[i][1] = antip_trs[3 * i + 1];
        vec4[i][2] = antip_trs[3 * i + 2];
    }
    // The offsets (starting at pos 1) of the RNA transcripts
    int antip_s2 = (antip_trs.size() - 1) / 3;
    // create a vector of string
    vector<string> vec5(antip_s2, "abc");
    // fill in the vector
    for (int i = 0; i < antip_s2; ++i)
    {
        vec5[i][0] = antip_trs[3 * i + 1];
        vec5[i][1] = antip_trs[3 * i + 2];
        vec5[i][2] = antip_trs[3 * i + 3];
    }
    // The offsets (starting at pos 2) of the RNA transcripts
    int antip_s3 = (antip_trs.size() - 2) / 3;
    // create a vector of string
    vector<string> vec6(antip_s3, "abc");
    // fill in the vector
    for (int i = 0; i < antip_s3; ++i)
    {
        vec6[i][0] = antip_trs[3 * i + 2];
        vec6[i][1] = antip_trs[3 * i + 3];
        vec6[i][2] = antip_trs[3 * i + 4];
    }
    // create a vector of vector of string
    // fill in the board with the six vectors above
    vector<vector<string>> board = {vec4, vec5, vec6, vec1, vec2, vec3};
    /*
    this is to test the result - print out the table and double check
    for(int i = 0; i < 6; ++i)
    {
        int boardsize = board[i].size();
        for (int j = 0; j < boardsize; ++j)
        {
           cout << board[i][j]<<", ";
        }
        cout << endl;
    }
*/
    return board;
}

/*
This function translates/converts a vector<string> (vector of codons)
into a string of amino acids using the genetic code

For example, the codons:
{"UUU", "GCC", "CAA"}
translates to:
F (Phenylalanine), A (Alanine), Q (Glutamine)
abreviated:
FAQ
*/
std::string Translate(const std::vector<std::string> &codon_sequence)
{
    // define a string variable (for now, empty)
    string result = "";
    // assign condon_sequence to a vector of string,
    // so that it can be modified
    vector<string> codon_seq = codon_sequence;
    // here's a list of the possible codons:
    // provided in bio.h
    vector<string> psb_codon = {
    "GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
    "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
    "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
    "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
    "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
    "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
    "GUG", "UAG", "UGA", "UAA"};
    // there corresponding amino acids ("*" represents STOP codons):
    // given in bio.h
    vector<string> amino_acids = {
    "A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
    "C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
    "I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
    "P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
    "Y", "V", "V", "V", "V", "*", "*", "*"};
    // the size of codon sequence
    int codon_seq_size = codon_seq.size();
    // the size of possible codons
    int psb_codon_size = psb_codon.size();
    // loop over codon sequence and possible codons;
    // when there is a match,
    // use the index to find the corresponding str in amino acids
    for (int i = 0; i < codon_seq_size; ++i)
    {
        for (int j = 0; j < psb_codon_size; ++j)
        {
            if (codon_seq[i] == psb_codon[j])
            {
                result += amino_acids[j];
            }
        }
    }
    return result;
}

/*
This function takes a DNA sequence and returns the longest possible
amino acid sequence / protein that is encoded by that sequence
(open reading frame). 

A valid open reading frame begins with the
codon AUG (the amino acid, Methionine (M)) and runs until 
a stop codon (*) is encountered. 

There may be multiple open reading frames in a sequence, 
and I need to check all six reading frames in order given by
get_reading_frames_as_codons.

If there are ties for longest, favor the first one found.

Return the longest open reading frame as an amino acid sequence. 
It must start with an 'M' and end with a '*' with no other '*'s within.
*/

std::string GetLongestOpenReadingFrame(const std::string &DNA_sequence)
{
    // define a string variable (for now, empty)
    string finalresult;
    // assign DNA_sequence to a string variable, so that it's editable
    string DNA_seq = DNA_sequence;
    // define a string variable raw DNA sequence
    string DNA_seq_raw;
    // call function GetReadingFramesAsCodons
    // assign the result to a vector of vector of string
    vector<vector<string>> result1 = GetReadingFramesAsCodons(DNA_seq);
    // the size of result1
    int result1s = result1.size();
    // create a vector of string
    vector<string> vec;
    // loop each RNA reading frame
    // translate them into string of amino acids
    // then push them back to a vector of string
    for (int i = 0; i < result1s; ++i)
    {
        vector<string> result2 = result1[i];
        DNA_seq_raw = Translate(result2);
        vec.push_back(DNA_seq_raw);
    }
    // this gives the longest possible amino acid sequence
    int longer = 0; // initialize an int as 0
    for (int i = 0; i < 6; ++i)
    {                                // loop over the six vectors of string
        int vecsize = vec[i].size(); // the size of each string
        for (int k = 0; k < vecsize; ++k)
        { // loop over each char
            if (vec[i][k] == 'M')
            { // starts comparison once hit 'M'
                for (int j = k; j < vecsize; ++j)
                { // stop when hit '*'
                    if (vec[i][j] == '*')
                    {
                        if (j - k <= longer)
                        {          // if it's not longer
                            break; // continue to next increment
                        }
                        else
                        {
                            longer = j - k; // else, store it to finalresult
                            finalresult = vec[i].substr(k, (j - k + 1));
                            break; // break the nearest for loop
                        }
                    }
                }
            }
        }
    }
    return finalresult;
}
