// Copyright (C) 2006-2011 by Cecile Ane

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License (the file gpl.txt included with this
// distribution or http://www.gnu.org/licenses/gpl.txt) for more
// details.

// File:     mdl.C

// Input:    [Options]

// Options:  type mdl -h 
// Output:   a file (default name: "mdl.out") listing the best partitions,
//           their DL criterion, their number of blocks, and the location
//           of breakpoints between blocks.
//           a file (default name: "mdl.mb") with a MrBayes block that
//           defines the character sets of the best partition.

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>

#include "mdl.h"

using namespace std;

string VERSION = "1.1";
string DATE = "February 2013";
string COPYRIGHT = "Copyright (C) 2007-2013 by Cecile Ane";

// readFile
//
// In version 1.0, 
// input file is expected to be of the following form:
// Tree    Length
// 1       2930
// Tree    Length
// 1       4385
// Tree    Length   etc., as made by paup (through mdl.pl)
//
// In version 1.1, new input file. Blocks no longer come in a specific order.
// Input file expected to show "block number : tree score" like this:
// 2 : 470
// 9 : 315
// 5 : 1138
// 1 : 284
// 3 : 644
// 6 : 1299
// 19 : 146

struct TreePair {
  TreePair(int t, double s) : treenb(t), score(s) {}
  // treenb: number (ID) for the block of segments,
  // score:  parsimony score of the best tree for that block.

  int treenb; 
  double score;
};

bool treeCmp( const TreePair a, const TreePair b ) {
	    return a.treenb < b.treenb ;
}

void readFile(vector<double> &groupscore, RunParameters &runpar)
{
  // temporary vector used for sorting
  vector<TreePair> tmpVec;
  
  string filename = runpar.getScoreFileName(); 
  ifstream f(filename.c_str());
  char tempc;
  if(f.fail()) {
    cerr <<"Error: Cannot open file " << filename << "." << endl;
    exit(1);
  }
  cout<<"Reading file "<<filename<<" ... "<<flush;
  bool line_even = false;
  while(f.good() && !f.eof()) {
    int treenb;
    string colon;
    double score;
    bool done;
    //if (line_even){
    f >> treenb >> colon >> score;
    if (!f.eof()) {
      tmpVec.push_back(TreePair(treenb, score));
    }
    //cout.precision(8);
    //cout<<fixed<<" read group score = "<<score<<endl;
    //} /* then discard the end of the line */
    do {f.get(tempc);} while (tempc !='\n' && !f.eof() && f.good());
    if(f.fail()) break;
    else done = false;
    //line_even = !line_even;
  }
  f.close();
  unsigned long Ngroups = tmpVec.size();
  runpar.setTotalNgroups(Ngroups);

  // sort the vector of TreePairs
  sort(tmpVec.begin(), tmpVec.end(), treeCmp);
  // print out sorted pairs and put them into groupscore vector
  for (int i = 0; i < tmpVec.size(); i++) {
    //cout << tmpVec[i].treenb << " " << tmpVec[i].score << endl;
    groupscore.push_back(tmpVec[i].score);
  }

  cout<< "done." << endl;
}

void usage()
{
  cout << "Usage: mdl " << endl;
  cout << "[-ntax      number-of-taxa]" << endl;
  cout << "[-nchar     number-of-characters]           (total alignment length)" << endl;
  cout << "[-ncharbase number-of-characters-per-block] (except for the last block)" << endl;
  cout << "[-nletters  number-of-letters]              (default to 4 for DNA)" << endl;
  cout << "[-ngroupmax maximum-number-of-groups-in-partitions] (default to max value)" << endl;
  cout << "[-nbestpart maximum-number-of-partitions-retained] " << endl;
  cout << "[-s         score-bound]                    (partitions with higher" << endl;
  cout << "                                             scores will not be kept)" << endl;
  cout << "[-scorefile file-name-for-group-scores]     (default to 'scores'. File produced by mdl.pl)" << endl;
  cout << "[-o         output-file]                    (default to 'mdl.out' and 'mdl.mb')" << endl;
  cout << "[-overhead  one-group-overhead]             (default is determined by Ntax)" <<endl;
  cout << "[-h] or [--help]                            (help message then quits)" << endl;
  exit(1);
}

void intro(ostream& f) {
  f << "\nMinimum Description Length for genome partitioning" << endl;
  f << "version " << VERSION << ", " << DATE << endl;
  f << COPYRIGHT << endl << endl;
}

int readArguments(int argc, char *argv[], DataParameters &datapar, RunParameters &runpar)
{
  int k=1;
  bool done = (argc>1 ? false : true);
  while(!done && k<argc) {
    string flag=argv[k];
    if(flag=="-h" || flag=="--help")
      usage();
    else if(flag=="-ntax") {
      double ntax;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> ntax) )
	usage();
      if(ntax <= 0.0)
	cerr << "Warning: parameter ntax must be positive." 
	     << " Ignoring argument -ntax " << ntax << "." << endl;
      else
	datapar.setNtax(ntax);
      k++;
    }
    else if(flag=="-nchar") {
      double nchar;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> nchar) )
	usage();
      if(nchar <= 0.0)
	cerr << "Warning: parameter nchar must be positive." 
	     << " Ignoring argument -nchar " << nchar << "." << endl;
      else
	datapar.setNchar(nchar);
      k++;
    }
    else if(flag=="-ncharbase") {
      unsigned long ncharbase;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> ncharbase) )
	usage();
      runpar.setNcharBase(ncharbase);
      k++;
    }
    else if(flag=="-nletters") {
      double nletters;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> nletters) )
	usage();
      if(nletters <= 0.0)
	cerr << "Warning: parameter nletters must be positive." 
	     << " Ignoring argument -nletters " << nletters << "." << endl;
      else
	datapar.setNletters(nletters);
      k++;
    }
    else if(flag=="-ngroupmax") {
      unsigned long ngroupmax;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> ngroupmax) )
	usage();
      runpar.setNgroupMax(ngroupmax);
      k++;
    }
    else if(flag=="-nbestpart") {
      int nbestpart;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> nbestpart) )
	usage();
      if(nbestpart <= 0.0)
	cerr << "Warning: parameter nbestpart must be positive." 
	     << " Ignoring argument -nbestpart " << nbestpart << "." << endl;
      else
      runpar.setNbestPart(nbestpart);
      k++;
    }
    else if(flag=="-s") {
      double score;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> score) )
	usage();
      if(score <= 0.0)
	cerr << "Warning: parameter s must be positive." 
	     << " Ignoring argument -s " << score << "." << endl;
      else
	runpar.setScoreBound(score);
      k++;
    }
    else if(flag=="-scorefile") {
      string scorefile = argv[++k];
      runpar.setScoreFileName(scorefile);
      k++;
    }
    else if(flag=="-o") {
      string outputfile = argv[++k];
      runpar.setOutputFile(outputfile);
      runpar.setMrbayesFile(outputfile + ".mb");
      k++;
    }
    else if(flag=="-overhead") {
      double overhead;
      string num = argv[++k];
      istringstream f(num);
      if( !(f >> overhead) )
	usage();
      if(overhead <= 0.0){
	cerr << "Warning: parameter overhead must be positive. It is recommended not "
	     << "to tune this parameter by hand. However, this option can be used to "
	     << "test the sensitivity of the analysis to this parameter. In this case, "
	     << "it is recommended to choose values of the order of the number of taxa. "
	     << endl;
	exit(1);
      }
      datapar.setonegroup_overhead(overhead);
      k++;
    }
    else if(flag=="-datafile") { // this is not used actually. 
      string datafile = argv[++k];
      datapar.setDataname(datafile);
      k++;
    }
    else
      done = true;
  }

  if( argc-k>0 )
    usage();
  return k;
}

void update_smallest(double &Ycandidate, vector<double> &Ybest,
		     unsigned long &j, vector<unsigned long> &jbest,
		     unsigned long &Ngcandidate, vector<unsigned long> &Ngbest,
		     int &kcandidate, vector<int> &kbest, bool add2vectors){
  // Ybest ought to be a vector already sorted, smallest (=best) values with largest 
  // indices. Ycandidate will replace a value in Ybest if it is < than at least one of
  // the Ybest values. If Ycandidate replaces a value in Ybest, then j Ngcandidate and
  // kcandidate should replace a value of jbest, Ngbest and kbest 
  // following the same permutation of indices.
  int k, k0, s;
  vector<double> temp;

  s = Ybest.size();
  if (add2vectors){
    if (s){
      for (k0=0; k0<s; k0++){
	if (Ybest[k0]<=Ycandidate) break;
      } // now k0 is the first index such that Ybest[k0]<=Ycandidate, and
      // if there is no such index, then k0 is the size of Ybest
      // and Ycandidate will just be appended to Ybest.
      Ybest.resize(s+1);
      jbest.resize(s+1);
      Ngbest.resize(s+1);
      kbest.resize(s+1);
      for (k=s; k>k0; k--){
	Ybest[k]=Ybest[k-1];
	jbest[k]=jbest[k-1];
	Ngbest[k]=Ngbest[k-1];
	kbest[k]=kbest[k-1];
      }
      Ybest[k0] = Ycandidate;
      jbest[k0]=j;
      Ngbest[k0]=Ngcandidate;
      kbest[k0]=kcandidate;
    } 
    else { // Ybest was empty
      Ybest.push_back(Ycandidate);
      jbest.push_back(j);
      Ngbest.push_back(Ngcandidate);
      kbest.push_back(kcandidate);
    }
  }
  else{
    for (k0=0; k0<s; k0++){
      if (Ybest[k0]<=Ycandidate) break;
    } // now k0 is the first index such that Ybest[k0]<=Ycandidate, and
      // if there is no such index, then k0 is the size of Ybest.
    for (k=0; k<k0-1; k++){
      Ybest[k] = Ybest[k+1];
      jbest[k]=jbest[k+1];
      Ngbest[k]=Ngbest[k+1];
      kbest[k]=kbest[k+1];
    }
    if (k0){
      Ybest[k0-1] = Ycandidate; 
      jbest[k0-1] = j;
      Ngbest[k0-1]= Ngcandidate;
      kbest[k0-1] = kcandidate;
    }
  }
  //cout<<"Ybest  is: "; for (k=0; k<Ybest.size(); k++){cout<<Ybest[k] <<" ";} cout<<endl;
  //cout<<"jbest  is: "; for (k=0; k<jbest.size(); k++){cout<<jbest[k] <<" ";} cout<<endl;
  //cout<<"Ngbest is: "; for (k=0; k<Ngbest.size();k++){cout<<Ngbest[k]<<" ";} cout<<endl;
}

void BestPartitions::fillup(vector<double> &gscore, DataParameters &datapar, RunParameters &runpar){
  cout <<"Filling in the best partitions..."<<flush;

  Group temp;
  temp.setNblocks(runpar.getTotalNblocks());
  unsigned long i,j, sb, eb;
  double DLcandidate, scoreMin = runpar.getScoreBound();
  vector<double> DLbest;    // vector of maximum size NbestPart
  vector<unsigned long> jbest, Ngbest;
  vector<int> kbest;
  int k,  ksmall;
  unsigned long Ngcandidate, Ngroupmax = runpar.getNgroupMax();
  bool add2vectors;

  for (i=0; i<numBlock; i++){
    eb=startBlock+i;
    sb=startBlock;

    DLbest.resize(0);
    jbest.resize(0);
    Ngbest.resize(0);
    kbest.resize(0);

    for (j=0; j<=i; j++){
      ksmall = (j==0 ? 1 : mdl[j-1].size()) ;
      for (k=0; k<ksmall; k++){
	Ngcandidate = 1+ (j==0 ? 0 : ngroups[j-1][k]);
	if (Ngcandidate <= Ngroupmax){
	  DLcandidate = (j==0 ? 0 : mdl[j-1][k])
	    + gscore[temp.getIndex(sb,eb)]+datapar.getonegroup_overhead();

	  add2vectors = ( DLbest.size() < runpar.getNbestPart() );
	  if (scoreMin<0 || DLcandidate<= scoreMin){ 
	    update_smallest(DLcandidate, DLbest, j, jbest, Ngcandidate,Ngbest,k,kbest, add2vectors);
	  }
	}
      }
      sb++; // should equal startblock+j
    }
    mdl[i].resize(DLbest.size());
    lastgroupstart[i].resize(jbest.size());
    ngroups[i].resize(Ngbest.size());
    k2goback[i].resize(kbest.size());
    for (k=0; k<mdl[i].size(); k++){
      mdl[i][k] = DLbest[k]; // copies the vector. Resizing included.
      lastgroupstart[i][k] = jbest[k]; // vector also.
      ngroups[i][k] = Ngbest[k];
      k2goback[i][k] = kbest[k];
    }
  }
  //this->print();
  cout<<" done."<<endl;
}

void BestPartitions::reconstruct(DataParameters &datapar, RunParameters &runpar,
				 ostream &fout, ostream &fmbout) const{
  cout << "Reconstructing the best partitions..."<< flush;
  fout<<"MDLscore Ngroups startingChar_list"<<endl<<fixed;
  fout.precision(8);
  int Npart = mdl[numBlock-1].size(); /* number of best partitions to reconstruct */
  Partition part;
  unsigned long j, i, ng, g, Ng;
  int p, k;

  for (p=Npart-1; p>=0; --p){
    /* best partitions are on the right, with large indices */
    part.reset();
    Ng=ngroups[numBlock-1][p];
    part.setNgroups(Ng);
    part.setMDLscore(mdl[numBlock-1][p]);

    j=numBlock-1;
    k=p;
    for (i=0; i<Ng; ++i){
      g  = lastgroupstart[j][k];
      ng = Ng-1-i;
      part.setgroup(ng, g);
      if(ng){
	k = k2goback[j][k];
	j = g-1; /* the order in which we update j and k is highly important */
      }
    } 
    /* at the end we should have that j=-1, k=0 and ngroups[j][k]=1
       before last update of j and k. */

    part.print(fout, runpar.getNcharBase(), datapar.getstepcost());
    if (p == (Npart-1)) // comment out this line and the next to suppress the .mb file.
      part.printMrBayesBlock(fmbout, runpar.getNcharBase(), datapar.getNchar());
  }
  cout<<" done."<<endl;
}


int main(int argc, char *argv[]){
  time_t beginTime, endTime;
  time(&beginTime);
  intro(cout);
  vector<double> groupscore;

  /* set default parameters, using 4 letters for DNA without gaps */
  //DataParameters datapar(5.0,69039.0,4.0);
  DataParameters datapar(4.0); // default to nletters=4
  RunParameters runpar(0);     // dedault to totalNblocks=0


  /* Read arguments from the command line */
  /* return value k is the position of the first input file, in case we want to 
     use input file after options.  Then argv[k] is the first input file */
  int k=readArguments(argc,argv,datapar,runpar);
  readFile(groupscore, runpar);
  runpar.checkNcharBase(datapar);

  ofstream fout((runpar.getOutputFile()).c_str());
  if(fout.fail()) {
    cerr <<"Error: Cannot open file " << runpar.getOutputFile() << "." << endl;
    exit(1);
  }
  fout<<"[";
  intro(fout);
  datapar.print(fout);
  runpar.print(fout);
  fout << "Program initiated at " << ctime(&beginTime) << flush;
  fout<<"]"<<endl;

  
  ofstream fmbout((runpar.getMrbayesFile()).c_str());
  if(fmbout.fail()) {
    cerr <<"Error: Cannot open file " << runpar.getMrbayesFile() << "." << endl;
    exit(1);
  }
  

  BestPartitions bestpart(1UL,runpar.getTotalNblocks());
  bestpart.fillup(groupscore,datapar,runpar);
  bestpart.reconstruct(datapar,runpar,fout,fmbout);
  fmbout.close();

  time(&endTime);
  fout << "[\nProgram ended at " << ctime(&endTime) << flush;
  int diff=endTime-beginTime;
  int days=diff/(24*60*60), hours=diff%(24*60*60)/(60*60), minutes=diff%(60*60)/60,seconds=diff%60;
  fout << "Elapsed time: ";
  if(days>0)                         fout<<days<< (days==1 ? " day, " : " days, ");
  if(days>0 || hours>0)              fout<<hours<< (hours==1 ? " hour, " : " hours, ");
  if(days>0 || hours>0 || minutes>0) fout<<minutes<< (minutes==1 ? " minute, " : " minutes, ");
  fout << seconds << (seconds==1 ? " second." : " seconds.") << "\n]"<<endl;
  fout.close();

  return 0;
}
