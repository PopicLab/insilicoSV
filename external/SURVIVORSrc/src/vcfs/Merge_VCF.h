/*
 * Merge_VCF.h
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#ifndef MERGE_VCF_H_
#define MERGE_VCF_H_
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include "../structs.h"
std::vector<strvcfentry> parse_vcf(std::string filename);
strcoordinate parse_stop(const char * buffer);
void merge_vcf(std::string filenames, int max_dist, std::string outputfile);
short get_type(std::string type) ;
int overlap(strvcfentry tmp, std::vector<strvcfentry> & final_vcf,int max_dist);
strcoordinate parse_stop(const char * buffer);
std::pair <bool,bool>parse_strands(const char * buffer);
std::vector<std::string> parse_filename(std::string filename);

#endif /* MERGE_VCF_H_ */
