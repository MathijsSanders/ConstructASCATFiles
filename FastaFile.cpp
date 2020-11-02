//----------------------------------------------------------------
// Name        : FastaFile.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : 
//----------------------------------------------------------------
#include <fstream>
#include <iostream>
#include "FastaFile.h"
//----------------------------------------------------------------
CFastaFile::CFastaFile(void)
{
    pSequence=NULL;
}
//----------------------------------------------------------------
CFastaFile::~CFastaFile(void)
{
    Clear();
}
//----------------------------------------------------------------
bool CFastaFile::Open(const string & sFileName)
{
    char *              pSequence;
    char *              pTarget;
    const char *        pLine;
    
    TARGET              target;

    ifstream            fastaFile;
    streampos           fileSize;       
    string              sLine;
    string              sTargetName;
    
    map<string,TARGET>  targets;
    
    //----------------------------------------------------------------
    //Remove previous sequences
    //----------------------------------------------------------------
    
    Clear();
    
    //----------------------------------------------------------------
    //Open the fasta file
    //----------------------------------------------------------------
    
    fastaFile.open(sFileName.c_str(),ifstream::in);
    
    if(fastaFile.good()==false)
    {
        cerr << "Error: Could not open fasta file!" << endl;
        return false;
    }
    
    //----------------------------------------------------------------
    //Compute file size
    //----------------------------------------------------------------
    
    fastaFile.seekg(0,ios_base::end);
    fileSize=fastaFile.tellg();
    fastaFile.seekg(0,ios_base::beg);
    
    //----------------------------------------------------------------
    //Allocate memory 
    //----------------------------------------------------------------
    
    pSequence=pTarget=new char[fileSize];
    
    //----------------------------------------------------------------    
    //Read fasta file 
    //----------------------------------------------------------------
    
    for(getline(fastaFile,sLine);fastaFile.good();)
    {
        if(sLine[0]=='>')
        {
            sTargetName=sLine.substr(1);
            target.pTarget=pTarget;
            
            for(getline(fastaFile,sLine);fastaFile.good();getline(fastaFile,sLine))
            {
                if(sLine[0]=='\0')  //Skip empty lines
                {
                    continue;
                }

                if(sLine[0]=='>')   //New target?
                {
                    break;
                }
                
                for(pLine=sLine.c_str();*pLine;)
                {
                    *pTarget++=*pLine++;
                }
            }
            
            target.iTargetLength=pTarget-target.pTarget;
            targets[sTargetName]=target;
        }
    }
    
    //----------------------------------------------------------------
    //End of file?
    //----------------------------------------------------------------
    
    if(fastaFile.eof()==false)
    {
        cerr << "Error: Could not read fasta file!" << endl;
        delete [] pSequence;
        return false;
    }
    
    //----------------------------------------------------------------
    //Save results
    //----------------------------------------------------------------
    
    this->pSequence=pSequence;
    this->targets=targets;
    
    //----------------------------------------------------------------
    //Done
    //----------------------------------------------------------------
    
    return true;
}
//----------------------------------------------------------------
void CFastaFile::Clear(void)
{
    if(pSequence)
    {
        delete [] pSequence;
        pSequence=NULL;
    }
    
    targets.clear();
}
//----------------------------------------------------------------