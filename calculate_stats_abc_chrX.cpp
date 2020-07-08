#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstring>
#include <list>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <functional>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <time.h>

using namespace std;

std::vector<double> statsCal (vector<string> haplotypes, unordered_map<string, int> hapsmap1, unordered_map<string, int> hapsmap2) {
  std::vector<double> stats(18);
  unsigned hapnum1 = hapsmap1.size();
  unsigned hapnum2 = hapsmap2.size();
  double hapfreq1 = 0;
  double hapfreq2 = 0;
  double hapheter1 = 1;
  double hapheter2 = 1;
  double hapfreqmost1 = 0;
  double hapfreqmost2 = 0;
  
  for (auto it = hapsmap1.begin(); it != hapsmap1.end(); ++it) {
    hapfreq1 = (double)it -> second / (double)8;
    if (hapfreq1 > hapfreqmost1) {
      hapfreqmost1 = hapfreq1;
    }
    hapheter1 = hapheter1 - hapfreq1 * hapfreq1;
  }
  
  for (auto it = hapsmap2.begin(); it != hapsmap2.end(); ++it) {
    hapfreq2 = (double)it -> second / (double)19;
    if (hapfreq2 > hapfreqmost2) {
      hapfreqmost2 = hapfreq2;
    }
    hapheter2 = hapheter2 - hapfreq2 * hapfreq2;
  }

  int snp1 = 0;                  // unique snps in population #1
  int snp2 = 0;
  int snpshared = 0;
  int s1 = 0;                    // number of polymorphic sites within population, >= 2 to calculate Tajima's D
  int s2 = 0;
  
  // pi = (2 * j * (n - j)) / (n * (n - 1)), n is number samples (chromosomes), j is count for alternative allele, j = AC & n = AN in vcf file
  double pi1 = 0;
  double theta1 = 0;
  // a1 = 1/1 + 1/2 + 1/3 ... + 1/S, S is number of segregating sites
  double a1 = 0;
  double a2 = 0;
  double b1 = 0;
  double b2 = 0;
  double c1 = 0;
  double c2 = 0;
  double e1 = 0;
  double e2 = 0;
  double tajimasd1 = 0;
  
  // pi = (2 * j * (n - j)) / (n * (n - 1)), n is number samples (chromosomes), j is count for alternative allele, j = AC & n = AN in vcf file
  double pi2 = 0;
  double theta2 = 0;
  // a1 = 1/1 + 1/2 + 1/3 ... + 1/S, S is number of segregating sites
  double a12 = 0;
  double a22 = 0;
  double b12 = 0;
  double b22 = 0;
  double c12 = 0;
  double c22 = 0;
  double e12 = 0;
  double e22 = 0;
  double tajimasd2 = 0;

  double pairb = 0;
  double pairw1 = 0;
  double pairw2 = 0;
  double fst = 0;

  string line1 = "";
  string line2 = "";
  line1 = haplotypes.at(0);
  std::vector<std::vector<int> > snps1(line1.length(), std::vector<int>(16));
  std::vector<std::vector<int> > snps2(line1.length(), std::vector<int>(28));
  
  vector<int> snps1values(line1.length());
  vector<int> snps2values(line1.length());
  
  double xy1 = 0;
  double xsquare1 = 0;
  double ysquare1 = 0;
  double rsquare1 = 0;
  double xy2 = 0;
  double xsquare2 = 0;
  double ysquare2 = 0;
  double rsquare2 = 0;
  int id = 0;
  int id2 = 0;
  double d1 = 0;
  double d2 = 0;
  double d3 = 0;
  double d4 = 0;

  for (unsigned i = 0; i < haplotypes.size(); i++) {
    line1 = haplotypes.at(i);
    for (unsigned k = 0; k < line1.length(); k++) {
      if (i < 8) {
	snps1.at(k).at(2 * i) = line1.at(k) - '0';
	snps1.at(k).at(2 * i + 1) = line1.at(k) - '0';
	snps1values.at(k) += line1.at(k) - '0';
	snps1values.at(k) += line1.at(k) - '0';
      }
      
      if (i >= 8) {
	if (i <= 17) {
	  snps2.at(k).at(i - 8) = line1.at(k) - '0';
	  snps2values.at(k) += line1.at(k) - '0';
	}
	else {
	  snps2.at(k).at(i - 8 + (i - 18)) = line1.at(k) - '0';
	  snps2.at(k).at(i - 8 + (i - 18) + 1) = line1.at(k) - '0';
	  snps2values.at(k) += line1.at(k) - '0';
	  snps2values.at(k) += line1.at(k) - '0';
	}
      }
    }
  }

  for (unsigned i = 0; i <  snps1.size(); i++) {
    if (snps1values.at(i) != 0 && snps1values.at(i) != 16) {
      s1++;
    }

    if (snps2values.at(i) != 0 && snps2values.at(i) != 28) {
      s2++;
    }

    if (snps1values.at(i) == 0 && snps2values.at(i) != 0 && snps2values.at(i) != 28) {
      snp2++;
    }

    if (snps2values.at(i) == 0 && snps1values.at(i) != 0 && snps1values.at(i) != 16) {
      snp1++;
    }

    if (snps1values.at(i) != 0 && snps2values.at(i) != 0) {
      snpshared++;
    }
  }

  if (s1 > 1 && s2 > 1) {
    for (unsigned i = 0; i < haplotypes.size() - 1; i++) {
      line1 = haplotypes.at(i);
      for (unsigned j = i + 1; j < haplotypes.size(); j++) {
	line2 = haplotypes.at(j);
	for (unsigned k = 0; k < line1.length(); k++) {
	  if (line1.at(k) != line2.at(k)) {
	    if( i < 8 && j >= 8) {
	      pairb += 1;
	    }
	    else {
	      if(i < 8 && j < 8) {
		pairw1 += 1;
	      }
	      else if (i >= 8 && j >= 8) {
		pairw2 += 1;
	      }
	    }
	  }
	}
      }
    }
    
    // Tajima's D
    for (int i = 1; i < 8; i++) {
      a1 += (double)1 / (double)i;
      a2 += (double)1 / (double)(i * i);
    }
    
    for (int i = 1; i < 19; i++) {
      a12 += (double)1 / (double)i;
      a22 += (double)1 / (double)(i * i);
    }
    
    b1 = (double)(8 + 1) / (double)(3 * 7); 
    b2 = (double)2 * (double)(8 * 8 + 8 + 3) / ((double)9 * (double)8 * (double)(8 - 1)); 
    c1 = b1 - (double)1 / a1;
    c2 = b2 - (double)(8 + 2) / (a1 * (double)8) + a2 / (a1 * a1); 
    e1 = c1 / a1; 
    e2 = c2 / (a1 * a1 + a2);
    pi1 = pairw1 * (double)2 / (double)(8 * 7);
    theta1 = (double)s1 / a1;
    
    tajimasd1 = (pi1 - theta1) / sqrt(abs(e1 * s1 + e2 * s1 * (s1 - (double)1)));
    
    b12 = (double)(19 + 1) / (double)(3 * 18); 
    b22 = (double)2 * (double)(19 * 19 + 19 + 3) / ((double)9 * (double)19 * (double)(19 - 1)); 
    c12 = b12 - (double)1 / a12;
    c22 = b22 - (double)(19 + 2) / (a12 * (double)19) + a22 / (a12 * a12); 
    e12 = c12 / a12; 
    e22 = c22 / (a12 * a12 + a22);
    pi2 = pairw2 * (double)2 / (double)(19 * 18);
    theta2 = (double)s2 / a12;

    tajimasd2 = (pi2 - theta2) / sqrt(abs(e12 * s2 + e22 * s2 * (s2 - (double)1)));

    pairb = pairb / (double)152;
    pairw1 = pairw1 / (double)28;
    pairw2 = pairw2 / (double)171;
    
    if (pairb > 0) {
      fst = (pairb - ((double)8 * pairw1 / (double)27 + (double)19 * pairw2 / (double)27)) / pairb;
    }
    else {
      fst = 0;
    }

    // Pearson R squared
    for (unsigned i = 0; i < snps1.size() - 1; i++) {
      for (unsigned j = i + 1; j < snps1.size(); j++) {
	xy1 = 0;
	xsquare1 = 0;
	ysquare1 = 0;
	xy2 = 0;
	xsquare2 = 0;
	ysquare2 = 0;

	for (unsigned k = 0; k < snps1.at(0).size() / 2; k++) {
	  xy1 += ((double)snps1.at(i).at(2*k) + (double)snps1.at(i).at(2*k+1) - (double)(snps1values.at(i) / (double)8)) * ((double)snps1.at(j).at(2*k) + (double)snps1.at(j).at(2*k+1) - (double)(snps1values.at(j) /  (double)8));
	  xsquare1 += ((double)snps1.at(i).at(2*k) + (double)snps1.at(i).at(2*k+1) - (double)(snps1values.at(i) /  (double)8)) * ( (double)snps1.at(i).at(2*k) + (double)snps1.at(i).at(2*k+1)- (double)(snps1values.at(i) / (double)8));
	  ysquare1 += ((double)snps1.at(j).at(2*k) + (double)snps1.at(j).at(2*k+1) - (double)(snps1values.at(j) /  (double)8)) * ((double)snps1.at(j).at(2*k) + (double)snps1.at(j).at(2*k+1) - (double)(snps1values.at(j) / (double)8));
	}

	for (unsigned k = 0; k < snps2.at(0).size() / 2; k++) {
	  xy2 += ((double)snps2.at(i).at(2*k) + (double)snps2.at(i).at(2*k+1) - (double)(snps2values.at(i) / (double)14)) * ((double)snps2.at(j).at(2*k) + (double)snps2.at(j).at(2*k+1) - (double)(snps2values.at(j) / (double)14));
	  xsquare2 += ((double)snps2.at(i).at(2*k) + (double)snps2.at(i).at(2*k+1) - (double)(snps2values.at(i) / (double)14)) * ((double)snps2.at(i).at(2*k) + (double)snps2.at(i).at(2*k+1) - (double)(snps2values.at(i) / (double)14));
	  ysquare2 += ((double)snps2.at(j).at(2*k) + (double)snps2.at(j).at(2*k+1) - (double)(snps2values.at(j) / (double)14)) * ((double)snps2.at(j).at(2*k) + (double)snps2.at(j).at(2*k+1) - (double)(snps2values.at(j) / (double)14));
	}
	
	if (xsquare1 * ysquare1 != 0) {
	  rsquare1 += xy1 * xy1 / (xsquare1 * ysquare1);
	  id++;
	}
	
	if (xsquare2 * ysquare2 != 0) {
	  rsquare2 += xy2 *xy2 / (xsquare2 * ysquare2);
	  id2++;
	}
      }
    }

    if (id > 0 ) {
      rsquare1 = rsquare1 / (double)id;
    }
    else {
      rsquare1 = -1000000;
    }
    if (id2 > 0) {
      rsquare2 = rsquare2 / (double)id2;
    }
    else {
      rsquare2 = -1000000;
    }

    stats.at(0) = snp1;
    stats.at(1) = snp2;
    stats.at(2) = snpshared;
    stats.at(3) = pi1;
    stats.at(4) = pi2 ;
    stats.at(5) = theta1;
    stats.at(6) = theta2;
    stats.at(7) = tajimasd1;
    stats.at(8) = tajimasd2;
    stats.at(9) = hapnum1;
    stats.at(10) = hapnum2;
    stats.at(11) = hapfreqmost1;
    stats.at(12) = hapfreqmost2;
    stats.at(13) = hapheter1;
    stats.at(14) = hapheter2;
    stats.at(15) = rsquare1;
    stats.at(16) = rsquare2;
    stats.at(17) = fst;
  }

  else {
    for (unsigned i = 0; i < stats.size(); i++)
      stats.at(i) = -1000000;
  }

  vector<vector<int>>().swap(snps1);
  vector<vector<int>>().swap(snps2);
  vector<int>().swap(snps1values);
  vector<int>().swap(snps2values);
  vector<string>().swap(haplotypes);
  unordered_map<string, int>().swap(hapsmap1);
  unordered_map<string, int>().swap(hapsmap2);
  
  return stats;
}


int main(int argc, char *argv[])
{
  if (argc != 3) {
    cout << "Usage: copy infile outfile" << endl;
    return -1;
  }
  
  ifstream inFile(argv[1]);
  if (!inFile)
    cout << argv[1] << "could not be opened" << endl;

  ofstream outFile(argv[2]);
  if (!outFile)
    cout << argv[2] << "could not be opened" << endl;
  
  string line = "";
  string oneline = "";
  unordered_map<string,int> haps1;
  unordered_map<string,int> haps2;
  vector<string> lines; 
  std::vector<std::vector<double> > sumStats(364, std::vector<double>(36));
  std::vector<double> statsMeans(18);
  std::vector<double> statsVariances(18);
  
  int id = 0;
  int rep = 0;
  
  while (getline(inFile, oneline)) {
    if (oneline.find_first_not_of("01") == std::string::npos && !oneline.empty()) {
      lines.push_back(oneline);
      if (id % 27 < 8) {
	haps1[oneline]++;
      }
      else {
	haps2[oneline]++;
      }

      if (id % 27 == 26) {
	sumStats.at(rep++) = statsCal(lines, haps1, haps2);
	vector<string>().swap(lines);
	unordered_map<string,int>().swap(haps1);
	unordered_map<string,int>().swap(haps2);
      }
      id++;
    }   
  }

  id = rep;
  rep = 0;

  for (unsigned i = 0; i < statsMeans.size(); i++) {
    statsMeans.at(i) = 0;
    statsVariances.at(i) = 0;
  }
  
  for (unsigned i = 0; i < id; i++) {
    if (sumStats.at(i).at(15) != -1000000 && sumStats.at(i).at(16) != -1000000) {
      rep++;
      for (unsigned j = 0; j < sumStats.at(0).size(); j++) {
	statsMeans.at(j) += sumStats.at(i).at(j);
      }
    }
  }

  for (unsigned i = 0; i < statsMeans.size(); i++) {
    statsMeans.at(i) = statsMeans.at(i) / (double)rep;
  }
  
  for (unsigned i = 0; i < id; i++) {
    if (sumStats.at(i).at(15) != -1000000 && sumStats.at(i).at(16) != -1000000) {
      for (unsigned j = 0; j < sumStats.at(0).size(); j++) {
	statsVariances.at(j) += (sumStats.at(i).at(j) - statsMeans.at(j)) * (sumStats.at(i).at(j) - statsMeans.at(j));
      }
    }
  }

  for (unsigned i = 0; i < statsMeans.size(); i++) {
    statsVariances.at(i) = statsVariances.at(i) / (double)rep;
    outFile << statsMeans.at(i) << "\t" << statsVariances.at(i) << "\t";
  }
  
  outFile << rep << endl;
  
  inFile.close();
  outFile.close();
  return 0;
}
