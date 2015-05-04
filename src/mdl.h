#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;


double lg(double n){
  // returns the smallest integer >= the log of n in base 2
  if (n<=0){
    cout<<"Error: no log of a negative number\n";
    return 0;
  } else{
    return ceil(log2(n));
  }
}

double ln(double n){
  // this one is the largest integer <= log of n in base 2.
  if (n<=0){
    cout<<"Error: no log of a negative number\n";
    return 0;
  } else{
    return floor(log2(n));
  }
}

double des(double n){
  // returns the number of bits required to describe n,
  // when there is no information about how large n is.
  // Uses logarithmic ramp. See Li & Vitanyi 1997. 
 double res = 1;
 while (n>=2){
   n=ln(n);
   res += n+1;
 }
 return res;
}

double log2UB(double n){
  // log2 of the # of Unrooted Binary trees with n taxa.
  double res = 0;
  while (n>3){
    res += log2(2*n-5);
    n -=1;
  }
  return(res);
}

class DataParameters {
 public:
  DataParameters(double nl) : Nletters(nl) {
    //b = lg(Nletters);
    b = log2(Nletters);
    Nchar = 0;
  }
  double getNchar() const{ return Nchar;}
  double getNletters() const{ return Nletters;}
  double getb() const{ return b;}
  double getlgNedge() const{ return lgNedge;}
  double getonegroup_overhead() const{ return onegroup_overhead;}
  double getstepcost() const{ return stepcost;}
  string getDataname() const{ return Dataname;}
  void setNchar(double n) { Nchar=n;}
  void setNtax(double n) {
    // b needs to be correct at this point, but we don't care about Nchar
    Ntax=n;
    lgNedge = log2(2*Ntax-3);
    onetreedes = log2UB(Ntax);
    stepcost = b + lgNedge;
    onegroup_overhead = log2UB(Ntax+1)/stepcost;
    //lgNedge = lg(2*Ntax-3);  
    //onetreedes = 2*Ntax-4+Ntax*lg(Ntax);
    //onegroup_overhead = (onetreedes + lgNedge)/stepcost;
    // that one is aproximately Ntax when Ntax is large, 
    // so that we separate 2 groups when Delta L > Ntax approximately.
  }
  void setNletters(double n) { 
    Nletters=n; 
    b=log2(Nletters); //b=lg(Nletters);
  }
  void setonegroup_overhead(double x){ onegroup_overhead=x;}
  void setDataname(string s) { Dataname=s;}
  void print(ostream &f) const{
    f<<"Data: "<<Dataname<<"\nNtax="<<Ntax<<" Nchar="<<Nchar
     <<" Nletters="<<Nletters<<" b="<<b<<endl;
    f<<"lgNedge="<<lgNedge<<", onetreedes="<<onetreedes<<", stepcost="<<
      stepcost<<", onegroup_overhead="<<onegroup_overhead<<endl;
  }
 private:
  double Ntax, Nchar, Nletters;
  double b;          // number of bits required to describe a single letter.
  double lgNedge;    // number of bits required to describe one edge of the tree
  double onetreedes; // number of bits required to describe the first tree.
  double onegroup_overhead; // (onetreedes + lgNedge) / (b+lgNedge)
  double stepcost; // b+lgNedge
  string Dataname;
};

class RunParameters {
 public:
  RunParameters(unsigned long n): TotalNblocks(n) {
    describeTreesSeparately = true;
    scorefilename = "scores";  /* default input  file name  */
    outputfile = "mdl.out";    /* default output file names */
    mrbayesfile = "mdl.mb";
    NgroupMax=0; NbestPart=0; NcharBase=1; scoreBound=-1;
  }
  void setdescribeTreesSeparately(bool x){ describeTreesSeparately = x; }
  void setTotalNblocks(unsigned long n) {TotalNblocks = n;}
  void setTotalNgroups(unsigned long n) {
    TotalNgroups = n;
    double Ng = static_cast<double>(TotalNgroups);
    TotalNblocks = static_cast<unsigned long>(sqrt(1/4 + 2*Ng)-1/2);
    if (NgroupMax==0 || TotalNblocks < NgroupMax) NgroupMax = TotalNblocks; 
    // by default, partitions may have as many groups as there are blocks
  }
  void setNgroupMax(unsigned long n){NgroupMax = n;}
  void setNbestPart(int n){NbestPart = n;}
  void setNcharBase(unsigned long n){NcharBase = n;}
  void checkNcharBase(DataParameters &dp){
    double p = dp.getNchar() / NcharBase - TotalNblocks+1;
    if (p>1 || p<=0){
      cerr << "Bad NcharBase="<<NcharBase<<" while I found "<<TotalNblocks<<" unit blocks"
	   <<" and Nchar="<<dp.getNchar()<<" (there should be about "<<p<<" blocks)"<<endl;
      exit(1);
    }
  }
  void setScoreFileName(string s){scorefilename = s;}
  void setOutputFile(string s){ outputfile  = s;}
  void setMrbayesFile(string s){mrbayesfile = s;}
  void setScoreBound(double x){scoreBound = x;}
  unsigned long getTotalNblocks() const{ return TotalNblocks;}
  unsigned long getTotalNgroups() const{ return TotalNgroups;}
  unsigned long getNgroupMax() const{ return NgroupMax;}
  unsigned long getNcharBase() const{ return NcharBase;}
  double getScoreBound() const{ return scoreBound;}
  int getNbestPart() const{ return NbestPart;}
  string getScoreFileName() const{ return scorefilename;}
  string getOutputFile()  const{ return outputfile; }
  string getMrbayesFile() const{ return mrbayesfile;}
  bool getdescribeTreesSeparately() const{ return describeTreesSeparately; }
  void print(ostream &f) const{
    f<<"NcharBase="<<NcharBase<<" TotalNblocks="<<TotalNblocks
     <<" TotalNgroups="<<TotalNgroups
     <<"\nNgroupMax="<<NgroupMax<<" NbestPart="<<NbestPart;
    if(scoreBound>=0)
      f<<" ScoreBound="<<scoreBound;
    f<<"\nscorefilename="<<scorefilename<<endl;
  }
 private:
  unsigned long TotalNblocks, TotalNgroups, NgroupMax, NcharBase;
  /* only partitions with a maximum of NgroupMax groups will be considered */
  int NbestPart; /* Number of best partitions to keep (first best, second best, etc.) */
  bool describeTreesSeparately;
  double scoreBound; /* Partitions with higher scores will not be kept */
  /* Name of input file containing scores of all groups */
  string scorefilename; 
  string outputfile;
  string mrbayesfile;
};



double blockDL(unsigned long gi, vector<double> &s, DataParameters &datapar){
  // Returns the contribution of the group with index gi
  // to the Description Length.
  return  datapar.getonegroup_overhead() + s[gi];
}


class Group {
 public:
  Group(unsigned long sb = 1, unsigned long eb = 1, unsigned long n = 1):
    startblock(sb), endblock(eb), nblocks(n){
    index = ((sb-1)*(2*n - sb))/2 + eb -1;
  }
  void setStartblock(unsigned long sb){startblock = sb;}
  void setEndpoint(unsigned long eb){endblock = eb;}
  void setNblocks(unsigned long n){nblocks = n;}
  void setIndex(unsigned long i){
    index = i;
    unsigned long s = 1;
    // while ((s*(2*nblocks-s+1))/2 <= i){ s++;}
    double N = static_cast<double>(nblocks);
    double sapprox = (1-sqrt(1-8*i/((2*N+1)*(2*N+1))))*(2*N+1)/2;
    s = static_cast<unsigned long>(sapprox);
    if ((s*(2*nblocks-s+1))/2 <= i){ s++;}
    startblock = s;
    endblock = i - ((s-1)*(2*nblocks - s))/2 + 1;
    // endblock has to be >= startblock and <= nblocks.
  }
  void updateIndex(){
     index = ((startblock-1)*(2*nblocks - startblock))/2 +endblock -1;
  }
  unsigned long getStartblock() const { return startblock; }
  unsigned long getEndblock() const { return endblock; }
  unsigned long getIndex() const { return index; }
  unsigned long getIndex(unsigned long sb, unsigned long eb) const {
    return (((sb-1)*(2*nblocks - sb))/2 +eb -1);
  }
  void print() const{
    cout<<"start and end blocks are: "<<startblock<<" and "<< endblock<<endl;
    cout<<"nblocks (i.e. TotalNblocks) is: "<<nblocks<<endl;
    cout<<"index is: "<<index<<endl;
  }
 private:
  unsigned long index, startblock, endblock;
  unsigned long nblocks;
  // index start at 0, and go up to nblocks*(nblocks+1)/2
  // nblocks is the total number of possible blocks, not the 
  // number of blocks actually in the group.
};

class Partition{
 public:
  Partition(unsigned long ng=0):ngroups(ng){
    startingBlock.resize(ng);
  }
  void setNgroups(unsigned long &n){
    ngroups=n;
    startingBlock.resize(n);
  }
  void reset(){ ngroups=0; startingBlock.resize(0);}
  void setgroup(unsigned long &ng, unsigned long &g){ startingBlock[ng] = g;}
  void setMDLscore(const double &s){ score=s;}

  double getScore() const{ return score;}
  unsigned long getNgroups() const{ return ngroups;}
  void print()const{
    cout << "score="<<score<<", ngroups=" <<ngroups<< ", starting blocks are"<<endl;
    for (unsigned long i=0; i<ngroups; i++){ cout<<" "<<startingBlock[i]; }
    cout<<endl;
  }
  void print(ostream &f, unsigned long nc, double stepcost) const{
    /* in case one wants to add the cost of describing the number of groups:
       double truescore = score + log2(ngroups)/stepcost;
       f<<score<<" "<<truescore<<" "<<ngroups;
    */
    f<<score<<" "<<ngroups;
    //for (unsigned long i=0; i<ngroups; i++){ f<<" "<<startingBlock[i]; }
    for (unsigned long i=0; i<ngroups; i++){ f<<" "<<1+startingBlock[i]*nc; }
    f<<endl;
  }
  void printMrBayesBlock(ostream &f, unsigned long nc, double Nchar) const{
    f<<"begin mrbayes;\n\n";
    for (unsigned long i=0; i<ngroups; i++)
      f<<"charset group"<<i+1<<" = "<<1+startingBlock[i]*nc<<"-"
       <<(i==(ngroups-1)? Nchar : startingBlock[i+1]*nc)<<";\n"; 
    f<<"\nend;\n"<<endl;
  }

 private:
  unsigned long ngroups;
  vector<unsigned long> startingBlock;
  /* StartingBlock = vector of size ngroups consisting of groups' starting blocks */
  double score;
};

class BestPartitions{
 public:
  BestPartitions(unsigned long sB=1, unsigned long nB=1):
    /* Default initialization: from block 1 to block nB which =1 by default.
       Top-left element of each matrix initialized. */
    startBlock(sB), numBlock(nB){
    endBlock = sB + nB-1;

    lastgroupstart.resize(nB);
    mdl.resize(nB);
    ngroups.resize(nB);
    k2goback.resize(nB);
  }
  unsigned long getnumBlock() const{ return numBlock;}
  void print() const{
    cout<<"startBlock="<<startBlock<<", endBlock="<<endBlock<<", numBlock="
	<<numBlock<<endl;
    for (unsigned int i=0; i<mdl.size(); i++){
      cout<<"MDLs for best partitions up to block "<<startBlock+i<<": ";
      for (int j=0; j<mdl[i].size(); j++){
	cout<<mdl[i][j]<<" ";
      }
      cout<<endl;
      cout<<"lastgroupsstart's for best partitions up to block "<<startBlock+i<<": ";
      for (int j=0; j<lastgroupstart[i].size(); j++){
	cout<<startBlock+lastgroupstart[i][j]<<" ";
      }
      cout<<endl;
      cout<<"ngroup's for best partitions up to block "<<startBlock+i<<": ";
      for (int j=0; j<ngroups[i].size(); j++){
	cout<<ngroups[i][j]<<" ";
      }
      cout<<endl;
      cout<<"k2goback's for best partitions up to block "<<startBlock+i<<": ";
      for (int j=0; j<k2goback[i].size(); j++){
	cout<<k2goback[i][j]<<" ";
      }
      cout<<endl;
    }
  }
  void fillup(vector<double> &, DataParameters &, RunParameters &);
  void reconstruct(DataParameters &, RunParameters &, ostream &, ostream &) const;

 private:
  unsigned long startBlock, endBlock, numBlock;
  vector<vector<unsigned long> > lastgroupstart;
  vector<vector<double> > mdl;
  vector<vector<unsigned long> >  ngroups;
  vector<vector<int> > k2goback;
  /* so mdl[i][kmax] is the Min Description Length from 
     the first best partition of blocks startBlock through startBlock+i.
     We need i<numBlocks.
     mdl[i][kmax-1] is for the second best, etc.
  */
};

double DL(Partition * partition, DataParameters * param){
  // Returns the Description Length of a data with "Ntax" taxa and
  // "Nchar" characters, using the partition.
  // The code assumes a separate description of trees,
  // not using nni distances between trees.
  unsigned int Ng = partition->getNgroups();
  double overhead, res;

  overhead = 2 * param->getb() * param->getNchar()
             + (Ng==1?0:des(static_cast<double>(Ng)));
  res = overhead + Ng * param->getonegroup_overhead()
                 + partition->getScore() * param->getstepcost();
  return res;
}
