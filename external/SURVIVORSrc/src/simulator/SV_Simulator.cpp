/*
 * SV_Simulator.cpp
 *
 *  Created on: Jan 30, 2016
 *      Author: fsedlaze
 */

#include "SV_Simulator.h"
#include <vector>
#include <stdint.h>


bool is_valid(char base) {
  return (((base == 'A' || base == 'C') || (base == 'R' || base == 'X')) || ((base == 'T' || base == 'G') || (base == 'N' || base == 'M')));
}

void check_genome(std::map<std::string, std::string> &genome, std::string msg) {
  std::cout << msg << " Genome checking:" << std::endl;

  for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
    for (size_t j = 1; j < (*i).second.size() + 1; j++) {
      if (!is_valid((*i).second[j - 1])) {
        std::cout << "err! " << (*i).second[j - 1] << std::endl;
      }
    }
  }
}
int parse_value(char* buffer, size_t buffer_size) { //required for parameter!
  int count = 0;
  for (size_t i = 1; i < buffer_size && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
    if (count == 1) {
      return atoi(&buffer[i]);
    }
    if (buffer[i] == ' ') {
      count++;
    }
  }
  return -1;
}

void simulate(std::map<std::string, std::string> genome, std::vector<struct_var> svs, std::string output_prefix) {

  std::cout << "apply SV" << std::endl;

}

parameter parse_param(std::string filename) {
  parameter tmp;
  size_t buffer_size = 200000;
  char*buffer = new char[buffer_size];
  std::ifstream myfile;

  myfile.open(filename.c_str(), std::ifstream::in);
  if (!myfile.good()) {
    std::cout << "Annotation Parser: could not open file: " << filename.c_str() << std::endl;
    exit(0);
  }

  myfile.getline(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.dup_min = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.dup_max = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.dup_num = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  std::cerr<<"DUP: "<< tmp.dup_min << " " << tmp.dup_max<< " " << tmp.dup_num << std::endl;

  tmp.indel_min = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.indel_max = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.indel_num = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
   std::cerr<<"INDEL: "<< tmp.indel_min << " " << tmp.indel_max<< " " << tmp.indel_num << std::endl;

  tmp.translocations_min = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.translocations_max = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.translocations_num = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);

  tmp.inv_min = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.inv_max = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.inv_num = parse_value(buffer, buffer_size);
  myfile.getline(buffer, buffer_size);
  tmp.inv_del_num = 0;

  if (!myfile.eof()) {
    tmp.inv_del_min = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.inv_del_max = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.inv_del_num = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
  }
  tmp.inv_dup_num = 0;
  if (!myfile.eof()) {
    tmp.inv_dup_min = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.inv_dup_max = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.inv_dup_num = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
  }
  tmp.intrachr_num = 0;
  if (!myfile.eof()) {
    tmp.intrachr_min = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.intrachr_max = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.intrachr_num = parse_value(buffer, buffer_size);
    //std::cout<<"NUM: "<<tmp.intrachr_num<<std::endl;
    myfile.getline(buffer, buffer_size);
  }
  tmp.small_dup_num = 0;
  if (!myfile.eof()) {
    tmp.small_dup_min = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.small_dup_max = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.small_dup_num = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
  }
  tmp.small_indel_num = 0;
  if (!myfile.eof()) {
    tmp.small_indel_min = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.small_indel_max = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.small_indel_num = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
  }
  tmp.small_inv_num = 0;
  if (!myfile.eof()) {
    tmp.small_inv_min = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.small_inv_max = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
    tmp.small_inv_num = parse_value(buffer, buffer_size);
    myfile.getline(buffer, buffer_size);
  }
  tmp.line_num = 0;
  if (!myfile.eof()) {
    tmp.line_num = parse_value(buffer, buffer_size);
    //myfile.getline(buffer, buffer_size);
  }
  std::cerr<<"NUM: "<< " SMALL_INV: " << tmp.small_inv_min << " " << tmp.small_inv_max << " " <<  tmp.small_inv_num << " SMALL_INDEL: " << tmp.small_indel_num << " SMALL_DUP: " << tmp.small_dup_num << " LINE: " << tmp.line_num << std::endl;
  myfile.close();
  return tmp;
}

std::map<std::string, std::string> read_fasta(std::string ref_file, int min_length) {
  std::string buffer;
  std::ifstream myfile;

  myfile.open(ref_file.c_str(), std::ifstream::in);
  if (!myfile.good()) {
    std::cout << "Annotation Parser: could not open file: " << ref_file.c_str() << std::endl;
    exit(0);
  }

  getline(myfile,buffer);
  std::map<std::string, std::string> genome;
  std::string seq;
  std::string name;
  while (!myfile.eof()) {
    if (buffer[0] == '>') {
      if (seq.size() > min_length) {
        genome[name] = seq;
      }
      name.clear();
      seq.clear();

      for (size_t i = 1; i < buffer.size() && buffer[i] != '\n' && buffer[i] != '\0' && buffer[i] != ' '; i++) {
        name += (buffer[i]);
      }
    } else {
      for (size_t i = 0; i < buffer.size() && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
        seq += toupper(buffer[i]);
      }
    }
    getline(myfile,buffer);
  }
  for (size_t i = 0; i < buffer.size() && buffer[i] != '\n' && buffer[i] != '\0'; i++) {
    seq += toupper(buffer[i]);
  }
  myfile.close();
  if (seq.size() > min_length) {
    genome[name] = seq;
  }
  seq.clear();
  name.clear();
  std::cout << "# Chrs passed size threshold:" << genome.size() << std::endl;
  return genome;
}

void sort_svs(std::vector<struct_var> svs) {
  std::map<std::string, std::vector<struct_var> > svs_tmp;
  for (size_t i = 0; i < svs.size(); i++) {
    svs_tmp[svs[i].pos.chr].push_back(svs[i]); //sort by chr:
  }

  for (std::map<std::string, std::vector<struct_var> >::iterator i = svs_tmp.begin(); i != svs_tmp.end(); i++) {
    std::vector<struct_var> tmp;
    if (!(*i).second.empty()) {
      tmp.push_back((*i).second[0]);
      for (size_t j = 0; j < (*i).second.size(); i++) {
        std::cout << (*i).second[j].pos.chr.c_str() << " " << (*i).second[j].pos.start << std::endl;
        size_t t = 0;
        while (tmp[t].pos.start < (*i).second[j].pos.start) {
          t++;
        }
        tmp.insert(tmp.begin() + t, (*i).second[j]);
      }
    }
    for (size_t j = 0; j < tmp.size(); j++) {
      std::cout << tmp[j].pos.chr.c_str() << " " << tmp[j].pos.start << std::endl;
    }
  }
}
float percent_N(std::string seq) {
  double n = 0;
  double size = (double) seq.size();
  for (size_t i = 0; i < seq.size(); i++) {
    if (seq[i] == 'N') {
      n++;
    }
  }
  //	std::cout<<"Percent: "<<n/size << " "<<n <<" " <<size<<std::endl;
  return n / size;
}
position get_pos(std::map<std::string, std::string> genome, int min_pos, int max_pos) {
  //std::cout<<max_pos<<std::endl;
  position tmp;
  std::string seq = "N";
  //std::cout<<percent_N(seq)<<std::endl;
  //std::cout<<"Start:"<<std::endl;
  for (size_t i = 0; i < 100 && percent_N(seq) > 0.05; i++) {
    std::map<std::string, std::string>::iterator chr = genome.begin();
    int id = rand() % genome.size(); //chose rand chr.
    int count = 0;
    while (chr != genome.end() && id != count) { //fast forward
      count++;
      chr++;
    }
    //std::cout<<"Place1"<<std::endl;
    tmp.chr = (*chr).first;
    //std::cout<<"Place1.1"<<std::endl;
    tmp.start = rand() % (((*chr).second.size() - max_pos)); //choose random start pos within chr length
    //std::cout<<"Place1.2"<<std::endl;
    if (max_pos == -1) {
      //insertion, translocation:
      //std::cout<<"Place1.3"<<std::endl;
      tmp.stop = max_pos;
      //std::cout<<"Place1.4"<<std::endl;
    } else {
      //std::cout<<"Place1.5"<<std::endl;
      tmp.stop = tmp.start + min_pos + (rand() % (max_pos - min_pos)); // choose stop location
      //std::cout<<"Place1.6"<<std::endl;
    }
    //std::cout<<"Place2"<<std::endl;
    int num = 0;

    while ((*chr).second.size() < tmp.stop && num < 100) { //choose chr,start, stop such that the mutation fits. Allow max 50 iterations!
      tmp.start = rand() % (((*chr).second.size() - max_pos)); //choose random start pos within chr length
      tmp.stop = tmp.start + min_pos + (rand() % (max_pos - min_pos)); // choose stop location
      num++;
    }
    //std::cout<<"Place3"<<std::endl;
    if (num == 100) {
      std::cerr << "Simulations are hindered by the two small chr size. " << std::endl;
      tmp.stop = -2;
    }
    if (max_pos != -1) {
      seq = (*chr).second.substr(tmp.start, tmp.stop - tmp.start);
    }

    //std::cerr<<i << " " << (*chr).first<<" "<<tmp.start<<" "<<tmp.stop<<std::endl;
  }
  //std::cout<<"end:"<<std::endl;
  return tmp;
}

position get_line_pos(std::vector<position>& line_records, int line_id) {
   if (line_id == -1) {
     // pick a random LINE entry from file
     line_id = rand() % line_records.size();
   }
   return line_records[line_id];
}


bool is_overlapping(position curr, std::vector<struct_var> svs) {

  for (size_t i = 0; i < svs.size(); i++) {
    if (strcmp(svs[i].pos.chr.c_str(), curr.chr.c_str()) == 0) {
      if (svs[i].pos.stop >= curr.start && svs[i].pos.start <= curr.stop) {
        return true;
      }
    }
  }
  return false;
}
position choose_pos(std::map<std::string, std::string> genome, int min, int max, std::vector<struct_var>& svs) {
  position pos = get_pos(genome, min, max);
  int num = 0;
  while (is_overlapping(pos, svs) && num < 100) {
    pos = get_pos(genome, min, max);
    num++;
  }
  if (num == 100) {
    std::cerr << "Terminate program as it could not find a non overlapping region: " << int(svs.size()) << " svs" << std::endl;
    exit(0);
  }
  return pos;
}
std::vector<struct_var> generate_mutations(std::string parameter_file, std::map<std::string, std::string> genome) {
  parameter par = parse_param(parameter_file);
  std::vector<struct_var> svs;
  //duplications
  struct_var mut;
  for (int i = 0; i < par.dup_num; i++) {
    mut.type = 0;
    //get_start location;
    mut.pos = choose_pos(genome, par.dup_min, par.dup_max, svs);
    //get_opposit location
    svs.push_back(mut);
  }
  //indels
  for (int i = 0; i < par.indel_num; i++) {
    //std::cout << "indel" << std::endl;
    if (rand() % 100 <= 50) {
      mut.type = 1; //insertion
    } else {
      mut.type = 4; //deletion
    }
    mut.pos = choose_pos(genome, par.indel_min, par.indel_max, svs);
    mut.target = mut.pos;
    svs.push_back(mut);
  }
  //inv
  for (int i = 0; i < par.inv_num; i++) {
    //	std::cout << "inv" << std::endl;
    mut.type = 2;
    mut.pos = choose_pos(genome, par.inv_min, par.inv_max, svs);
    mut.target = mut.pos;
    mut.target.start = mut.pos.stop;
    mut.target.stop = mut.pos.start;
    svs.push_back(mut);
  }

  //tra inter
  std::cout << par.intrachr_num << std::endl;
  for (int i = 0; i < par.intrachr_num; i++) {
    mut.type = 6;
    mut.pos = choose_pos(genome, par.translocations_min, par.translocations_max, svs);
    //std::cout<<i<<": "<<mut.pos.chr<<" "<<mut.pos.start<<" size: "<<mut.pos.stop-mut.pos.start<<std::endl;
    mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
    while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) != 0) {
      mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
    }

    //I need to be sure about the same lenght of the tra!:
    int size1 = mut.pos.stop - mut.pos.start;
    int size2 = mut.target.stop - mut.target.start;

    mut.pos.stop = mut.pos.start + std::min(size1, size2);
    mut.target.stop = mut.target.start + std::min(size1, size2);
    svs.push_back(mut);
  }
  //tra
  for (int i = 0; i < par.translocations_num; i++) {
    //	std::cout << "tra" << std::endl;
    mut.type = 3;
    mut.pos = choose_pos(genome, par.translocations_min, par.translocations_max, svs);
    //std::cout<<"size: "<<mut.pos.stop-mut.pos.start<<std::endl;
    mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
    while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) == 0) {
      mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
    }

    //I need to be sure about the same lenght of the tra!:
    int size1 = mut.pos.stop - mut.pos.start;
    int size2 = mut.target.stop - mut.target.start;

    mut.pos.stop = mut.pos.start + std::min(size1, size2);
    mut.target.stop = mut.target.start + std::min(size1, size2);
    svs.push_back(mut);
  }
  //complex inv_del
  for (int i = 0; i < par.inv_del_num; i++) {
    //1. sim:
    mut.type = 2;
    mut.pos = choose_pos(genome, par.inv_del_min, par.inv_del_max, svs);
    //2. determin size of del:
    int len = (mut.pos.stop - mut.pos.start) / 10; //dels are ~20% ofthe size!
    mut.pos.start += len;
    mut.pos.stop -= len;
    mut.target = mut.pos;
    mut.target.start = mut.pos.stop;
    mut.target.stop = mut.pos.start;

    svs.push_back(mut);

    struct_var del;
    //the del infront:
    del.type = 4;
    del.pos.chr = mut.pos.chr;
    del.pos.stop = mut.pos.start;
    del.pos.start = del.pos.stop - len;
    del.target = del.pos;
    svs.push_back(del);

    //the del behind:
    del.pos.start = mut.pos.stop;
    del.pos.stop = del.pos.start + len;
    del.target = del.pos;
    svs.push_back(del);
  }
  //inv dup
  for (int i = 0; i < par.inv_dup_num; i++) {
    mut.type = 5;
    //get_start location;
    mut.pos = choose_pos(genome, par.inv_dup_min, par.inv_dup_max, svs);
    //get_opposit location
    svs.push_back(mut);
  }
  //	sort_svs(svs);
  return svs;
}

std::vector<struct_var> generate_mutations_ref(std::string parameter_file, std::map<std::string, std::string> genome) {
  parameter par = parse_param(parameter_file);
  std::vector<struct_var> svs;
  //duplications
  struct_var mut;
  //indels
  for (int i = 0; i < par.indel_num; i++) {
    //std::cout << "indel" << std::endl;
    if (rand() % 100 <= 50) {
      mut.type = 1; //insertion
    } else {
      mut.type = 4; //deletion
    }
    mut.pos = choose_pos(genome, par.indel_min, par.indel_max, svs);
    mut.target = mut.pos;
    svs.push_back(mut);
  }
  //inv
  for (int i = 0; i < par.inv_num; i++) {
    //	std::cout << "inv" << std::endl;
    mut.type = 2;
    mut.pos = choose_pos(genome, par.inv_min, par.inv_max, svs);
    mut.target = mut.pos;
    mut.target.start = mut.pos.stop;
    mut.target.stop = mut.pos.start;
    svs.push_back(mut);
  }
  //tra
  for (int i = 0; i < par.translocations_num; i++) {
    //	std::cout << "tra" << std::endl;
    mut.type = 3;
    mut.pos = choose_pos(genome, par.translocations_min, par.translocations_max, svs);
    //std::cout<<"size: "<<mut.pos.stop-mut.pos.start<<std::endl;
    mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
    while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) == 0) {
      mut.target = choose_pos(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
    }

    //I need to be sure about the same lenght of the tra!:
    int size1 = mut.pos.stop - mut.pos.start;
    int size2 = mut.target.stop - mut.target.start;

    mut.pos.stop = mut.pos.start + std::min(size1, size2);
    mut.target.stop = mut.target.start + std::min(size1, size2);
    svs.push_back(mut);
  }

  return svs;
}

void store_sorted(std::vector<struct_var> &svs, struct_var tmp) {
  std::vector<struct_var>::iterator i = svs.begin();
  while (i != svs.end() && strcmp(tmp.target.chr.c_str(), (*i).target.chr.c_str()) >= 0 && tmp.target.start > (*i).target.start) {
    i++;
  }
  svs.insert(i, tmp);
}
void store_ins(std::vector<insertions> & ins, insertions tmp) {

  std::vector<insertions>::iterator i = ins.begin();
  while (i != ins.end() && tmp.target.start > (*i).target.start) {
    i++;
  }
  ins.insert(i, tmp);
}
char complement(char nuc) {
  switch (nuc) {
    case 'A':
      return 'T';
      break;
    case 'C':
      return 'G';
      break;
    case 'G':
      return 'C';
      break;
    case 'T':
      return 'A';
      break;
    default:
      return nuc;
      break;
  }
  return nuc;
}
void invert(std::string &seq) {
  std::string tmp;
  for (std::string::reverse_iterator i = seq.rbegin(); i != seq.rend(); i++) {
    tmp += complement((*i));
  }
  seq.clear();

  seq = tmp;
}
std::string rand_seq(int length) {
  std::string tmp;
  //tmp.resize((size_t) length);
  for (int i = 0; i < length; i++) {
    switch (rand() % 4) {
      case 0:
        tmp += 'A';
        break;
      case 1:
        tmp += 'C';
        break;
      case 2:
        tmp += 'G';
        break;
      case 3:
        tmp += 'T';
        break;
    }
  }
  return tmp;
}
void apply_mutations(std::map<std::string, std::string> &genome, std::vector<struct_var>& svs) {
  srand(time(NULL));
  std::vector<insertions> ins;
  insertions in;
  std::string seq1;
  std::string seq2;
  int pos;
  std::vector<struct_var> new_svs; //Thanks to the invdup we need this.
  //all mutations that do not change the coordinates are later applied.
  //all others are directly applied (e.g. INV, TRA)
  for (size_t i = 0; i < svs.size(); i++) {
    std::string tmp;
    switch (svs[i].type) {
      case 0:
        //duplication
        svs[i].seq = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
        in.seq = svs[i].seq;
        in.target = svs[i].pos;	//check
        store_ins(ins, in);
        break;
      case 1: { 
        //insertion:
        //svs[i].seq = rand_seq(svs[i].target.stop - svs[i].target.start);
        // simulate dispersed dup as an insertion
        int dup_start_pos = svs[i].pos.start;
        int dup_len = svs[i].target.stop - svs[i].target.start;
        int pos_delta = 3000 + rand() % (40000);
        int origin_start_pos = dup_start_pos - pos_delta;
        if (origin_start_pos >= 0) {
            std::cout << "INV-DDUP " << dup_start_pos << " " << dup_len << " " << pos_delta << " " << origin_start_pos << std::endl;
            svs[i].seq = genome[svs[i].pos.chr].substr(origin_start_pos, dup_len);
        
        } else {
            svs[i].seq = rand_seq(svs[i].target.stop - svs[i].target.start);
        }
        in.seq = svs[i].seq;
        in.target = svs[i].target;
        store_ins(ins, in);
        break; }
      case 2:
        //inversion
        tmp = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
        //std::cout<<"INV: "<<tmp.size()<<std::endl;
        invert(tmp);
        genome[svs[i].pos.chr].erase(svs[i].pos.start, tmp.size());
        genome[svs[i].pos.chr].insert(svs[i].pos.start, tmp);
        break;
      case 3:
        //translocations
        seq1 = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
        seq2 = genome[svs[i].target.chr].substr(svs[i].target.start, (svs[i].target.stop - svs[i].target.start));
        //	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;

        pos = 0;
        for (int j = svs[i].target.start; j < svs[i].target.stop; j++) {
          genome[svs[i].target.chr][j] = seq1[pos];
          pos++;
        }
        pos = 0;
        for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
          genome[svs[i].pos.chr][j] = seq2[pos];
          pos++;
        }
        break;
      case 4: //deletion: //just mark those regions
        //std::cout<<"DEL: "<<svs[i].pos.chr<<" "<<svs[i].pos.start <<" "<<svs[i].pos.stop<<" g: "<< genome[svs[i].pos.chr].size()<<std::endl;
        //	std::cout << "size: " << genome[svs[i].pos.chr].size() << " " << svs[i].pos.start << " " << (svs[i].pos.stop - svs[i].pos.start) << std::endl;
        for (size_t j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
          genome[svs[i].pos.chr][j] = 'X';
        }
        break;

      case 5:
        //invduplication first a duplication and then an inversion
        svs[i].seq = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
        in.seq = svs[i].seq;
        invert(in.seq);
        in.target = svs[i].pos;	//check
        store_ins(ins, in);
        svs[i].type=0;//dup
        new_svs.push_back(svs[i]);
        svs[i].type=2;//inv
        break;
      case 6:
        //inter tra:
        svs[i].seq = genome[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
        // = genome[svs[i].target.chr].substr(svs[i].target.start, (svs[i].target.stop - svs[i].target.start));
        //	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;

        in.seq = svs[i].seq;
        in.target = svs[i].target;
        store_ins(ins, in);
        pos = 0;
        for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
          genome[svs[i].pos.chr][j] = 'X';
          pos++;
        }
        break;
      default:
        break;
    }
  }

  for (std::vector<insertions>::reverse_iterator i = ins.rbegin(); i != ins.rend(); i++) {
    genome[(*i).target.chr].insert((*i).target.start, (*i).seq);
  }
  for(size_t i =0;i<new_svs.size();i++){
    svs.push_back(new_svs[i]);
  }

}

void apply_mutations_ref(std::map<std::string, std::string> &genome, std::vector<struct_var> &svs) {
  std::cout << "apply mut ref!" << std::endl;
  srand(time(NULL));
  std::vector<insertions> ins;
  insertions in;
  std::string seq1;
  std::string seq2;
  int pos;
  //all mutations that do not change the coordinates are later applied.
  //all others are directly applied (e.g. INV, TRA)
  //for (size_t i = 0; i < svs.size(); i++) {
  std::vector<struct_var> tmp_svs;
  for (size_t i = 0; i < svs.size(); i++) {
    store_sorted(tmp_svs, svs[i]);
  }
  svs = tmp_svs;

  //iterator
  for (std::vector<struct_var>::reverse_iterator i = svs.rbegin(); i != svs.rend(); i++) {
    if ((*i).type == 4) {
      //insertions: (simulated deletions) move always further away..???
      (*i).type = 1;

    } else if ((*i).type == 1) {
      (*i).type = 4;
    }
  }
  for (std::vector<struct_var>::iterator i = svs.begin(); i != svs.end(); i++) {
    std::cout << (*i).pos.start << " " << (*i).type << std::endl;
    std::string tmp;
    switch ((*i).type) {
      case 1:
        genome[(*i).pos.chr].erase((*i).pos.start, (*i).pos.stop - (*i).pos.start);
        break;
      case 2:
        //inversion
        tmp = genome[(*i).pos.chr].substr((*i).pos.start, ((*i).pos.stop - (*i).pos.start));
        //std::cout<<"INV: "<<tmp.size()<<std::endl;
        invert(tmp);
        genome[(*i).pos.chr].erase((*i).pos.start, tmp.size());
        genome[(*i).pos.chr].insert((*i).pos.start, tmp);
        break;
      case 3:
        //translocations
        seq1 = genome[(*i).pos.chr].substr((*i).pos.start, ((*i).pos.stop - (*i).pos.start));
        seq2 = genome[(*i).target.chr].substr((*i).target.start, ((*i).target.stop - (*i).target.start));
        //	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;
        pos = 0;
        for (int j = (*i).target.start; j < (*i).target.stop; j++) {
          genome[(*i).target.chr][j] = seq1[pos];
          pos++;
        }
        pos = 0;
        for (int j = (*i).pos.start; j < (*i).pos.stop; j++) {
          genome[(*i).pos.chr][j] = seq2[pos];
          pos++;
        }
        break;
      case 4:
        //deletions: (simulated insertions)
        (*i).seq = rand_seq((*i).target.stop - (*i).target.start);
        in.seq = (*i).seq;
        in.target = (*i).target;
        genome[in.target.chr].insert(in.target.start, in.seq);
        break;

      default:
        break;
    }
  }
}
void write_fasta(std::string output_prefix, std::map<std::string, std::string> genome) {
  std::string out = output_prefix;
  out += ".fasta";
  FILE *file2;
  file2 = fopen(out.c_str(), "w");
  if (file2 == NULL) {
    std::cout << "Error in printing: The file or path that you set " << output_prefix.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
    exit(0);
  }

  for (std::map<std::string, std::string>::iterator i = genome.begin(); i != genome.end(); i++) {
    fprintf(file2, "%c", '>');
    fprintf(file2, "%s", (*i).first.c_str());
    fprintf(file2, "%c", '\n');
    int len = 0;
    for (size_t j = 1; j < (*i).second.size() + 1; j++) {
      if (!is_valid((*i).second[j - 1])) {
        std::cout << "err! " << (*i).second[j - 1] << std::endl;
      }
      if ((*i).second[j - 1] != 'X') {
        fprintf(file2, "%c", (*i).second[j - 1]);
        len++;
      }
      if (len % 100 == 0) {
        fprintf(file2, "%c", '\n');
      }

    }
    if (len % 100 != 0) {
      fprintf(file2, "%c", '\n');
    }
  }
  //std::cout << std::endl;
  fclose(file2);
}
void write_sv(std::string output_prefix, std::vector<struct_var> svs) {
  std::string out = output_prefix;
  out += ".bed";
  FILE *file2;
  file2 = fopen(out.c_str(), "w");
  if (file2 == NULL) {
    std::cout << "Error in printing: The file or path that you set " << out.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
    exit(0);
  }

  out = output_prefix;
  out += ".insertions.fa";
  FILE *file;
  file = fopen(out.c_str(), "w");
  if (file == NULL) {
    std::cout << "Error in printing: The file or path that you set " << out.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
    exit(0);
  }

  for (size_t i = 0; i < svs.size(); i++) {
    if (svs[i].type == 1) { //write inserted sequeces to fasta file!
      fprintf(file, "%c", '>');
      fprintf(file, "%s", svs[i].pos.chr.c_str());
      fprintf(file, "%c", '_');
      fprintf(file, "%i", svs[i].pos.start);
      fprintf(file, "%c", '\n');
      fprintf(file, "%s", svs[i].seq.c_str());
      fprintf(file, "%c", '\n');
    }
    //write pseudo bed:
    if (svs[i].type == 3 || svs[i].type == 6) {
      fprintf(file2, "%s", svs[i].pos.chr.c_str());
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%i", svs[i].pos.start);
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%s", svs[i].target.chr.c_str());
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%i", svs[i].target.start);
      fprintf(file2, "%c", '\t');
      if (svs[i].type == 3) {
        fprintf(file2, "%s", "TRA\n");
      } else {
        fprintf(file2, "%s", "INTRATRA\n");
      }
      fprintf(file2, "%s", svs[i].pos.chr.c_str());
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%i", svs[i].pos.stop);
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%s", svs[i].target.chr.c_str());
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%i", svs[i].target.stop);
      fprintf(file2, "%c", '\t');
    } else {
      fprintf(file2, "%s", svs[i].pos.chr.c_str());
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%i", svs[i].pos.start);
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%s", svs[i].pos.chr.c_str());
      fprintf(file2, "%c", '\t');
      fprintf(file2, "%i", svs[i].pos.stop);
      fprintf(file2, "%c", '\t');
    }
    switch (svs[i].type) {
      case 0:
        fprintf(file2, "%s", "DUP");
        break;
      case 1:
        fprintf(file2, "%s", "INS");
        break;
      case 2:
        fprintf(file2, "%s", "INV");
        break;
      case 3:
        fprintf(file2, "%s", "TRA");
        break;
      case 4:
        fprintf(file2, "%s", "DEL");
        break;
      case 5:
        fprintf(file2, "%s", "INVDUP");
        break;
      case 6:
        fprintf(file2, "%s", "INTRATRA");
        break;
      default:
        break;
    }
    fprintf(file2, "%c", '\n');
  }
  fclose(file);
  fclose(file2);
}

void simulate_SV(std::string ref_file, std::string parameter_file, bool coordinates, std::string output_prefix) {

  //read in list of SVs over vcf?
  //apply vcf to genome?
  srand(time(NULL));
  parameter par = parse_param(parameter_file);
  int min_chr_len = std::max(std::max(par.dup_max, par.indel_max), std::max(par.inv_max, par.translocations_max));
  std::map<std::string, std::string> genome = read_fasta(ref_file, min_chr_len);
  check_genome(genome, "First:");
  std::cout << "generate SV" << std::endl;
  std::vector<struct_var> svs;
  if (coordinates) {
    //simulate reads
    svs = generate_mutations(parameter_file, genome);
    check_genome(genome, "Sec:");
    apply_mutations(genome, svs);	//problem: We need two different coordinates. Simulate once for one and then for the other???
  } else {
    svs = generate_mutations_ref(parameter_file, genome);
    check_genome(genome, "Sec:");
    apply_mutations_ref(genome, svs);	//problem: We need two different coordinates. Simulate once for one and then for the other???
  }
  check_genome(genome, "Last:");
  std::cout << "write genome" << std::endl;
  write_fasta(output_prefix, genome);
  std::cout << "write SV" << std::endl;
  write_sv(output_prefix, svs);
}

void generate_parameter_file(std::string parameter_file) {
  FILE *file2;
  file2 = fopen(parameter_file.c_str(), "w");
  if (file2 == NULL) {
    std::cerr << "Error in printing: The file or path that you set " << parameter_file.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
    exit(0);
  }
  fprintf(file2, "%s", "PARAMETER FILE: DO JUST MODIFY THE VALUES AND KEEP THE SPACES!\n");
  fprintf(file2, "%s", "DUPLICATION_minimum_length: 100\n");
  fprintf(file2, "%s", "DUPLICATION_maximum_length: 10000\n");
  fprintf(file2, "%s", "DUPLICATION_number: 3\n");

  fprintf(file2, "%s", "INDEL_minimum_length: 20\n");
  fprintf(file2, "%s", "INDEL_maximum_length: 500\n");
  fprintf(file2, "%s", "INDEL_number: 1\n");

  fprintf(file2, "%s", "TRANSLOCATION_minimum_length: 1000\n");
  fprintf(file2, "%s", "TRANSLOCATION_maximum_length: 3000\n");
  fprintf(file2, "%s", "TRANSLOCATION_number: 2\n");

  fprintf(file2, "%s", "INVERSION_minimum_length: 600\n");
  fprintf(file2, "%s", "INVERSION_maximum_length: 800\n");
  fprintf(file2, "%s", "INVERSION_number: 4\n");

  fprintf(file2, "%s", "INV_del_minimum_length: 600\n");
  fprintf(file2, "%s", "INV_del_maximum_length: 800\n");
  fprintf(file2, "%s", "INV_del_number: 2\n");

  fprintf(file2, "%s", "INV_dup_minimum_length: 600\n");
  fprintf(file2, "%s", "INV_dup_maximum_length: 800\n");
  fprintf(file2, "%s", "INV_dup_number: 2\n");

  fprintf(file2, "%s", "INTRA_TRANS_minimum_length: 600\n");
  fprintf(file2, "%s", "INTRA_TRANS_maximum_length: 800\n");
  fprintf(file2, "%s", "INTRA_TRANS_number: 2\n");
  fclose(file2);
}

/*******************************************************************/

void store_ins(std::vector<insertions_diploid> & ins, insertions_diploid tmp) { //change input type

  std::vector<insertions_diploid>::iterator i = ins.begin();
  while (i != ins.end() && tmp.target.start > (*i).target.start) {
    i++;
  }
  ins.insert(i, tmp);
}

bool is_overlapping_diploid(position curr, std::vector<struct_var_diploid> svs) { //change input type

  for (size_t i = 0; i < svs.size(); i++) {
    if (strcmp(svs[i].pos.chr.c_str(), curr.chr.c_str()) == 0) {
      if (svs[i].pos.stop >= curr.start && svs[i].pos.start <= curr.stop) {
        return true;
      }
    }
  }
  return false;
}
position choose_pos_diploid(std::map<std::string, std::string> genome, int min, int max, std::vector<struct_var_diploid>& svs) { //change input type
  position pos = get_pos(genome, min, max);
  int num = 0;
  while (is_overlapping_diploid(pos, svs) && num < 100) {
    pos = get_pos(genome, min, max);
    num++;
  }
  if (num == 100) {
    std::cerr << "Terminate program as it could not find a non overlapping region " << svs.size() << " svs" << std::endl;
    exit(0);
  }
  return pos;
}


position choose_line_pos_diploid(std::vector<position>& line_records, std::vector<struct_var_diploid>& svs, bool all_line) {
  int line_id = -1;
  if (all_line) {
      std::cout << " size " << svs.size() << " line " << line_id << "\n"; 
      line_id = svs.size();
  }
  position pos = get_line_pos(line_records, line_id);
  int num = 0;
  while (is_overlapping_diploid(pos, svs) && num < 100) {
    pos = get_line_pos(line_records, line_id);
    num++;
  }
  if (num == 100) {
    std::cerr << "Terminate program as it could not find a non overlapping region for a LINE event "  << svs.size() << " svs" << std::endl;
    exit(0);
  }
  return pos;
}



void apply_SNP(std::map<std::string, std::string>& genomeA,std::map<std::string, std::string> &genomeB, int bases_per_snp, FILE *file2A, FILE *file2B, FILE *file2AB) { //NEW

  //iterate thru chromosomes and apply SNPs
  FILE *tmpFile;
  for (std::map<std::string, std::string>::iterator i = genomeA.begin(); i != genomeA.end(); i++) { 
    std::string &strA = (*i).second;
    std::string chrname = (*i).first;
    std::string &strB = genomeB[chrname];
    for (size_t j = 1; j < strA.size() + 1; j++) {
      if ((!is_valid(strA[j - 1])) || (!is_valid(strB[j - 1])) || (strA[j - 1] != strB[j - 1])) {
        std::cout << "err! " << strA[j - 1] << strB[j - 1] << std::endl;
      }
      if (rand() % bases_per_snp == 0) { //introduce mutation here
        char snp;
        while ((snp = rand_seq(1)[0]) == strA[j - 1]){
          continue;
        }
        int choice = rand() % 3;
        if (choice == 0) { //modify A only
          strA[j - 1] = snp;
          tmpFile = file2A;
        } else if (choice == 1) { //modify B only
          strB[j - 1] = snp;
          tmpFile = file2B;
        } else { //homozygous choice == 2
          strA[j - 1] = snp;
          strB[j - 1] = snp;
          tmpFile = file2AB;
        }
        //std::cout << "snp " << strA[j - 1] << strB[j - 1] << std::endl;

        fprintf(tmpFile, "%s", chrname.c_str());
        fprintf(tmpFile, "%c", '\t');
        fprintf(tmpFile, "%zu", j);
        fprintf(tmpFile, "%c", '\t');
        fprintf(tmpFile, "%s", chrname.c_str());
        fprintf(tmpFile, "%c", '\t');
        fprintf(tmpFile, "%zu", (j+1) );
        fprintf(tmpFile, "%c", '\t');
        fprintf(tmpFile, "%s", "SNP\n");
      }    
    }
  }
}

int determine_haplotype(int* num_a, int* num_b, int limit) { //NEW
  int hap;
  if (*num_a == limit) { //if only B needs more events
    hap = 1;
    *num_b = *num_b + 1;
  } else if (*num_b == limit) { //if only A needs more events
    hap = 0;
    *num_a = *num_a + 1;
  } else { //event has opportunity to be homozygous
    if (rand() % 3 == 0) { //both with 1/3 chance
      hap = 2;
      *num_a = *num_a + 1;
      *num_b = *num_b + 1;
    } else if (rand() % 2 == 0) { //B only with 1/3 chance
      hap = 1;
      *num_b = *num_b + 1;
    } else { //A only with 1/3 chance
      hap = 0;
      *num_a = *num_a + 1;
    }
  }
  //std::cout << *num_a << *num_b << std::endl; 
  return hap;
}

std::vector<struct_var_diploid> generate_mutations_diploid(std::string parameter_file, std::map<std::string, std::string> genome, std::vector<position> &line_records) { //EDITED
  parameter par = parse_param(parameter_file);
  std::vector<struct_var_diploid> svs;

  //variables for throughout function
  struct_var_diploid mut;
  int num_a = 0;
  int num_b = 0;
  std::cerr << "generate_mutations_diploid function" << std::endl;
  //duplications
  while (num_a < par.dup_num or num_b < par.dup_num) {
    mut.type = 0;
    //get_start location;
    mut.pos = choose_pos_diploid(genome, par.dup_min, par.dup_max, svs);
    //std::cout << "chose position duplication" << std::endl;
    //get_opposit location
    //determine haplotype
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.dup_num);
    //std::cout << "got haplotype duplication" << std::endl;
    svs.push_back(mut);
    //std::cerr << num_a << " " << num_b << std::endl; 
  }
  std::cerr << "DUP " << num_a << num_b << std::endl; 
  num_a = 0;
  num_b = 0;
  while (num_a < par.small_dup_num or num_b < par.small_dup_num) {
    mut.type = 0;
    //get_start location;
    mut.pos = choose_pos_diploid(genome, par.small_dup_min, par.small_dup_max, svs);
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.small_dup_num);
    svs.push_back(mut);
  }
  std::cerr << "SMALL DUP " << num_a << num_b << std::endl; 
  num_a = 0;
  num_b = 0;
  //indels
  while (num_a < par.indel_num or num_b < par.indel_num) {
    //std::cout << "indel" << std::endl;
    if (rand() % 100 <= 20) {
      mut.type = 1; //insertion
    } else {
      mut.type = 4; //deletion
    }
    mut.pos = choose_pos_diploid(genome, par.indel_min, par.indel_max, svs);
    mut.target = mut.pos;
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.indel_num);
    svs.push_back(mut);
  }
  std::cerr << "INDEL " << num_a << num_b << std::endl; 

  num_a = 0;
  num_b = 0;
  while (num_a < par.small_indel_num or num_b < par.small_indel_num) {
    if (rand() % 100 <= 20) { //!!! update for dels
      mut.type = 1; //insertion
    } else {
      mut.type = 4; //deletion
    }
    mut.pos = choose_pos_diploid(genome, par.small_indel_min, par.small_indel_max, svs);
    mut.target = mut.pos;
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.small_indel_num);
    svs.push_back(mut);
  }
  std::cerr << "SMALL INDEL " << num_a << num_b << std::endl; 

  num_a = 0;
  num_b = 0;
  //inv
  while (num_a < par.inv_num or num_b < par.inv_num) {
    //std::cout << "inv" << std::endl;
    mut.type = 2;
    mut.pos = choose_pos_diploid(genome, par.inv_min, par.inv_max, svs);
    mut.target = mut.pos;
    mut.target.start = mut.pos.stop;
    mut.target.stop = mut.pos.start;
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.inv_num);
    svs.push_back(mut);
    //std::cerr << num_a << " " << num_b << std::endl;
  }
  std::cerr << "INV " << num_a << num_b << std::endl;

  num_a = 0;
  num_b = 0;
  while (num_a < par.small_inv_num or num_b < par.small_inv_num) {
    mut.type = 2;
    mut.pos = choose_pos_diploid(genome, par.small_inv_min, par.small_inv_max, svs);
    mut.target = mut.pos;
    mut.target.start = mut.pos.stop;
    mut.target.stop = mut.pos.start;
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.small_inv_num);
    svs.push_back(mut);  
  }
  std::cerr << "SMALL INV " << num_a << num_b << std::endl;

  std::cerr << " Simulating LINEs " << par.line_num << std::endl;
  // produce delections at specific positions in the genome 
  num_a = 0;
  num_b = 0;
  int n_total = 0;
  while ((num_a < par.line_num or num_b < par.line_num) and (n_total < line_records.size())) {
    mut.type = 4; //deletion 
    mut.pos = choose_line_pos_diploid(line_records, svs, false); //true);
    mut.target = mut.pos;
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.line_num);
    svs.push_back(mut);
    n_total += 1;
  }
  std::cerr << " Simulated LINEs " << num_a << num_b << std::endl;

  num_a = 0;
  num_b = 0;
  //tra inter
  //std::cout << par.intrachr_num << std::endl;
  while (num_a < par.intrachr_num or num_b < par.intrachr_num) {
    mut.type = 6;
    mut.pos = choose_pos_diploid(genome, par.translocations_min, par.translocations_max, svs);
    //std::cout<<num_a<<": "<<mut.pos.chr<<" "<<mut.pos.start<<" size: "<<mut.pos.stop-mut.pos.start<<std::endl;
    mut.target = choose_pos_diploid(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
    while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) != 0) {
      mut.target = choose_pos_diploid(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
    }

    //I need to be sure about the same lenght of the tra!:
    int size1 = mut.pos.stop - mut.pos.start;
    int size2 = mut.target.stop - mut.target.start;

    mut.pos.stop = mut.pos.start + std::min(size1, size2);
    mut.target.stop = mut.target.start + std::min(size1, size2);

    mut.haplotype=determine_haplotype(&num_a,&num_b,par.intrachr_num);
    svs.push_back(mut);
  }

  num_a = 0;
  num_b = 0;
  //tra
  while (num_a < par.translocations_num or num_b < par.translocations_num) {
    //std::cout << "tra" << std::endl;
    mut.type = 3;
    mut.pos = choose_pos_diploid(genome, par.translocations_min, par.translocations_max, svs);
    //std::cout<<"size: "<<mut.pos.stop-mut.pos.start<<std::endl;
    mut.target = choose_pos_diploid(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs); //TRA has to be of the same size!
    while (strcmp(mut.target.chr.c_str(), mut.pos.chr.c_str()) == 0) {
      mut.target = choose_pos_diploid(genome, mut.pos.stop - mut.pos.start, (mut.pos.stop - mut.pos.start) + 1, svs);
    }

    //I need to be sure about the same lenght of the tra!:
    int size1 = mut.pos.stop - mut.pos.start;
    int size2 = mut.target.stop - mut.target.start;

    mut.pos.stop = mut.pos.start + std::min(size1, size2);
    mut.target.stop = mut.target.start + std::min(size1, size2);

    mut.haplotype=determine_haplotype(&num_a,&num_b,par.translocations_num);
    svs.push_back(mut);
  }
  std::cerr << "TRA " << num_a << num_b << std::endl;

  num_a = 0;
  num_b = 0;
  //complex inv_del
  while (num_a < par.inv_del_num or num_b < par.inv_del_num) {
    //1. sim:
    mut.type = 2;
    mut.pos = choose_pos_diploid(genome, par.inv_del_min, par.inv_del_max, svs);
    //2. determin size of del:
    int len = (mut.pos.stop - mut.pos.start) / 10; //dels are ~20% ofthe size!
    mut.pos.start += len;
    mut.pos.stop -= len;
    mut.target = mut.pos;
    mut.target.start = mut.pos.stop;
    mut.target.stop = mut.pos.start;

    mut.haplotype=determine_haplotype(&num_a,&num_b,par.inv_del_num);
    svs.push_back(mut);

    struct_var_diploid del;
    del.haplotype = mut.haplotype;
    //the del infront:
    del.type = 4;
    del.pos.chr = mut.pos.chr;
    del.pos.stop = mut.pos.start;
    del.pos.start = del.pos.stop - len;
    del.target = del.pos;
    svs.push_back(del);

    //the del behind:
    del.pos.start = mut.pos.stop;
    del.pos.stop = del.pos.start + len;
    del.target = del.pos;

    svs.push_back(del);
  }
  std::cerr << "INVDEL " << num_a << num_b << std::endl;

  num_a = 0;
  num_b = 0;
  //inv dup
  for (int i = 0; i < par.inv_dup_num; i++) {
    mut.type = 5;
    //get_start location;
    mut.pos = choose_pos_diploid(genome, par.inv_dup_min, par.inv_dup_max, svs);
    //get_opposit location
    mut.haplotype=determine_haplotype(&num_a,&num_b,par.inv_dup_num);
    svs.push_back(mut);
  }
  std::cerr << "INVDUP " << num_a << num_b << std::endl;
  //	sort_svs(svs);
  return svs;
}

void write_fasta_diploid(std::string output_prefix, std::map<std::string, std::string> genomeA,std::map<std::string, std::string> genomeB) { //EDITED
  std::string outA = output_prefix;
  std::string outB = output_prefix;
  outA += "A.fasta";
  outB += "B.fasta";
  FILE *file2A;
  FILE *file2B;
  file2A = fopen(outA.c_str(), "w");
  file2B = fopen(outB.c_str(), "w");
  if (file2A == NULL || file2B == NULL) {
    std::cout << "Error in printing: The file or path that you set " << output_prefix.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
    exit(0);
  }

  for (std::map<std::string, std::string>::iterator i = genomeA.begin(); i != genomeA.end(); i++) { 
    fprintf(file2A, "%c", '>');
    fprintf(file2A, "%s", (*i).first.c_str());
    fprintf(file2A, "%c", '\n');

    fprintf(file2B, "%c", '>');
    fprintf(file2B, "%s", (*i).first.c_str());
    fprintf(file2B, "%c", '\n');

    std::string strA = (*i).second;
    std::string strB = genomeB[(*i).first];

    int len = 0;
    for (size_t j = 1; j < strA.size() + 1; j++) {
      if (!is_valid(strA[j - 1])) {
        std::cout << "err! " << strA[j - 1] << std::endl;
      }
      if (strA[j - 1] != 'X') {
        fprintf(file2A, "%c", strA[j - 1]);
        len++;
      }
      if (len % 100 == 0) {
        fprintf(file2A, "%c", '\n');
      }

    }
    if (len % 100 != 0) {
      fprintf(file2A, "%c", '\n');
    }

    len = 0;
    for (size_t j = 1; j < strB.size() + 1; j++) {
      if (!is_valid(strB[j - 1])) {
        std::cout << "err! " << strB[j - 1] << std::endl;
      }
      if (strB[j - 1] != 'X') {
        fprintf(file2B, "%c", strB[j - 1]);
        len++;
      }
      if (len % 100 == 0) {
        fprintf(file2B, "%c", '\n');
      }

    }
    if (len % 100 != 0) {
      fprintf(file2B, "%c", '\n');
    }
  }
  //std::cout << std::endl;
  fclose(file2A);
  fclose(file2B);
}



void apply_mutations_diploid(std::map<std::string, std::string>& genomeA, std::map<std::string, std::string>& genomeB, std::vector<struct_var_diploid>& svs) { //EDITED
  srand(time(NULL));
  std::vector<insertions_diploid> ins;
  insertions_diploid in;
  std::string seq1;
  std::string seq2;
  int pos;
  std::vector<struct_var_diploid> new_svs; //Thanks to the invdup we need this.
  //all mutations that do not change the coordinates are later applied.
  //all others are directly applied (e.g. INV, TRA)
  for (size_t i = 0; i < svs.size(); i++) {
    std::string tmp;
    switch (svs[i].type) {
      case 0:
        //duplication
        //separate homozygous for each genome b/c of SNPs
        if (svs[i].haplotype == 0 or svs[i].haplotype == 2) { //implement on A
          svs[i].seq = genomeA[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          in.haplotype = 0;
          in.seq = svs[i].seq;
          in.target = svs[i].pos;	//check
          store_ins(ins, in);
        }
        if (svs[i].haplotype == 1 or svs[i].haplotype == 2) { //implement on B
          svs[i].seq = genomeB[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          in.haplotype = 1;
          in.seq = svs[i].seq;
          in.target = svs[i].pos;	//check
          store_ins(ins, in);
        }

        break;
      case 1: {
        /*int dup_start_pos = svs[i].pos.start;
        int dup_len = svs[i].target.stop - svs[i].target.start;
        int pos_delta = 3000 + rand() % (40000);
        int origin_start_pos = dup_start_pos - pos_delta;
        if (origin_start_pos >= 0) {
            std::cout << "INV-DDUP " << dup_start_pos << " " << dup_len << " " << pos_delta << " " << origin_start_pos << std::endl;
            if (svs[i].haplotype == 0 or svs[i].haplotype == 2) {
                tmp = genomeA[svs[i].pos.chr].substr(origin_start_pos, dup_len);
                invert(tmp);
                svs[i].seq = tmp;
            } 
            if (svs[i].haplotype == 1 or svs[i].haplotype == 2) {
                svs[i].seq = genomeB[svs[i].pos.chr].substr(origin_start_pos, dup_len);
            }
        } else {
            std::cout << "INS " << std::endl;
            svs[i].seq = rand_seq(svs[i].target.stop - svs[i].target.start);
        }
        in.haplotype = svs[i].haplotype;
        in.seq = svs[i].seq;
        in.target = svs[i].target;
        store_ins(ins, in);
        break; */ 
        //insertion:
        svs[i].seq = rand_seq(svs[i].target.stop - svs[i].target.start);
        in.haplotype = svs[i].haplotype;
        in.seq = svs[i].seq;
        in.target = svs[i].target;
        store_ins(ins, in);
        break;
        }
      case 2:
        //inversion
        //std::cout<<"INV: "<<tmp.size()<<std::endl;
        //separate homozygous for each genome b/c of SNPs
        if (svs[i].haplotype == 0 or svs[i].haplotype == 2) { //implement on A
          tmp = genomeA[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          invert(tmp);
          genomeA[svs[i].pos.chr].erase(svs[i].pos.start, tmp.size());
          genomeA[svs[i].pos.chr].insert(svs[i].pos.start, tmp);
        } 
        if (svs[i].haplotype == 1 or svs[i].haplotype == 2) { //implement on B
          tmp = genomeB[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          invert(tmp);
          genomeB[svs[i].pos.chr].erase(svs[i].pos.start, tmp.size());
          genomeB[svs[i].pos.chr].insert(svs[i].pos.start, tmp);
        }
        break;
      case 3:
        //translocations
        //	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;

        //separate homozygous for each genome b/c of SNPs
        if (svs[i].haplotype == 0 or svs[i].haplotype == 2) { //implement on A
          seq1 = genomeA[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          seq2 = genomeA[svs[i].target.chr].substr(svs[i].target.start, (svs[i].target.stop - svs[i].target.start));

          pos = 0;
          for (int j = svs[i].target.start; j < svs[i].target.stop; j++) {
            genomeA[svs[i].target.chr][j] = seq1[pos];
            pos++;
          }
          pos = 0;
          for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
            genomeA[svs[i].pos.chr][j] = seq2[pos];
            pos++;
          }
        } 
        if (svs[i].haplotype == 1 or svs[i].haplotype == 2) { //implement on B
          seq1 = genomeA[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          seq2 = genomeA[svs[i].target.chr].substr(svs[i].target.start, (svs[i].target.stop - svs[i].target.start));

          pos = 0;
          for (int j = svs[i].target.start; j < svs[i].target.stop; j++) {
            genomeB[svs[i].target.chr][j] = seq1[pos];
            pos++;
          }
          pos = 0;
          for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
            genomeB[svs[i].pos.chr][j] = seq2[pos];
            pos++;
          }
        }
        break;
      case 4: //deletion: //just mark those regions
        //std::cout<<"DEL: "<<svs[i].pos.chr<<" "<<svs[i].pos.start <<" "<<svs[i].pos.stop<<" g: "<< genome[svs[i].pos.chr].size()<<std::endl;
        //	std::cout << "size: " << genome[svs[i].pos.chr].size() << " " << svs[i].pos.start << " " << (svs[i].pos.stop - svs[i].pos.start) << std::endl;
        if (svs[i].haplotype == 0 or svs[i].haplotype == 2) { //implement on A
          for (size_t j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
            genomeA[svs[i].pos.chr][j] = 'X';
          }
        } 
        if (svs[i].haplotype == 1 or svs[i].haplotype == 2) { //implement on B
          for (size_t j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
            genomeB[svs[i].pos.chr][j] = 'X';
          }
        }
        break;

      case 5:
        //invduplication first a duplication and then an inversion
        //separate homozygous for each genome b/c of SNPs
        if (svs[i].haplotype == 0 or svs[i].haplotype == 2) { //implement on A
          svs[i].seq = genomeA[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          in.haplotype = 0;
          in.seq = svs[i].seq;
          invert(in.seq);
          in.target = svs[i].pos;	//check
          store_ins(ins, in);
        }
        if (svs[i].haplotype == 1 or svs[i].haplotype == 2) { //implement on B
          svs[i].seq = genomeB[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          in.haplotype = 1;
          in.seq = svs[i].seq;
          invert(in.seq);
          in.target = svs[i].pos;	//check
          store_ins(ins, in);
        }
        //haplotype for svs[i] is set appropriately
        svs[i].type=0;//dup
        new_svs.push_back(svs[i]);
        svs[i].type=2;//inv
        break;
      case 6:
        //inter tra:
        //	std::cout<<"TRA: "<<seq1.size()<<" "<<seq2.size()<<std::endl;

        //separate homozygous insertions for each genome b/c of SNPs
        if (svs[i].haplotype == 0 or svs[i].haplotype == 2) { //implement on A
          svs[i].seq = genomeA[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          in.haplotype = 0;
          in.seq = svs[i].seq;
          in.target = svs[i].target;
          store_ins(ins, in);

          pos = 0;
          for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
            genomeA[svs[i].pos.chr][j] = 'X';
            pos++;
          }    
        } 
        if (svs[i].haplotype == 1 or svs[i].haplotype == 2) { //implement on B
          svs[i].seq = genomeB[svs[i].pos.chr].substr(svs[i].pos.start, (svs[i].pos.stop - svs[i].pos.start));
          in.haplotype = 1;
          in.seq = svs[i].seq;
          in.target = svs[i].target;
          store_ins(ins, in);

          pos = 0;
          for (int j = svs[i].pos.start; j < svs[i].pos.stop; j++) {
            genomeB[svs[i].pos.chr][j] = 'X';
            pos++;
          }
        }

        break;
      default:
        break;
    }
  }

  for (std::vector<insertions_diploid>::reverse_iterator i = ins.rbegin(); i != ins.rend(); i++) {
    if ((*i).haplotype == 0 or (*i).haplotype == 2) { //implement on A
      genomeA[(*i).target.chr].insert((*i).target.start, (*i).seq);
    } 
    if ((*i).haplotype == 1 or (*i).haplotype == 2) { //implement on B
      genomeB[(*i).target.chr].insert((*i).target.start, (*i).seq);
    }
  }
  for(size_t i =0;i<new_svs.size();i++){
    svs.push_back(new_svs[i]); //have correct haplotypes already
  }

}

void write_sv_diploid(std::string output_prefix, std::vector<struct_var_diploid> svs, FILE *file2A, FILE *file2B, FILE *file2AB) { //EDITED

  std::string outA = output_prefix;
  std::string outB = output_prefix;
  std::string outAB = output_prefix;
  outA += ".hetA.insertions.fa";
  outB += ".hetB.insertions.fa";
  outAB += ".homAB.insertions.fa";
  FILE *fileA;
  FILE *fileB;
  FILE *fileAB;
  fileA = fopen(outA.c_str(), "w");
  fileB = fopen(outB.c_str(), "w");
  fileAB = fopen(outAB.c_str(), "w");
  if (fileA == NULL || fileB == NULL || fileAB == NULL) {
    std::cout << "Error in printing: The file or path that you set " << outA.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
    exit(0);
  }

  FILE *tmpFile;
  for (size_t i = 0; i < svs.size(); i++) {
    if (svs[i].haplotype == 0) { //appropriate insertion file
      tmpFile = fileA;
    } else if (svs[i].haplotype == 1) {
      tmpFile = fileB;
    } else {
      tmpFile = fileAB;
    }
    if (svs[i].type == 1) { //write inserted sequeces to fasta file!  
      fprintf(tmpFile, "%c", '>');
      fprintf(tmpFile, "%s", svs[i].pos.chr.c_str());
      fprintf(tmpFile, "%c", '_');
      fprintf(tmpFile, "%i", svs[i].pos.start);
      fprintf(tmpFile, "%c", '\n');
      fprintf(tmpFile, "%s", svs[i].seq.c_str());
      fprintf(tmpFile, "%c", '\n');

    }
    if (svs[i].haplotype == 0) { //appropriate bed file
      tmpFile = file2A;
    } else if (svs[i].haplotype == 1) {
      tmpFile = file2B;
    } else {
      tmpFile = file2AB;
    }
    //write pseudo bed:
    if (svs[i].type == 3 || svs[i].type == 6) {
      fprintf(tmpFile, "%s", svs[i].pos.chr.c_str());
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%i", svs[i].pos.start);
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%s", svs[i].target.chr.c_str());
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%i", svs[i].target.start);
      fprintf(tmpFile, "%c", '\t');
      if (svs[i].type == 3) {
        fprintf(tmpFile, "%s", "TRA\n");
      } else {
        fprintf(tmpFile, "%s", "INTRATRA\n");
      }
      fprintf(tmpFile, "%s", svs[i].pos.chr.c_str());
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%i", svs[i].pos.stop);
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%s", svs[i].target.chr.c_str());
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%i", svs[i].target.stop);
      fprintf(tmpFile, "%c", '\t');
    } else {
      fprintf(tmpFile, "%s", svs[i].pos.chr.c_str());
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%i", svs[i].pos.start);
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%s", svs[i].pos.chr.c_str());
      fprintf(tmpFile, "%c", '\t');
      fprintf(tmpFile, "%i", svs[i].pos.stop);
      fprintf(tmpFile, "%c", '\t');
    }
    switch (svs[i].type) {
      case 0:
        fprintf(tmpFile, "%s", "DUP");
        break;
      case 1:
        fprintf(tmpFile, "%s", "INS");
        break;
      case 2:
        fprintf(tmpFile, "%s", "INV");
        break;
      case 3:
        fprintf(tmpFile, "%s", "TRA");
        break;
      case 4:
        fprintf(tmpFile, "%s", "DEL");
        break;
      case 5:
        fprintf(tmpFile, "%s", "INVDUP");
        break;
      case 6:
        fprintf(tmpFile, "%s", "INTRATRA");
        break;
      default:
        break;
    }
    fprintf(tmpFile, "%c", '\n');    
  }
  fclose(fileA);
  fclose(file2A);
  fclose(fileB);
  fclose(file2B);
  fclose(fileAB);
  fclose(file2AB);
}

std::vector<position> load_line_bed_file() {
  std::ifstream in_file("/athena/ihlab/scratch/vpopic/SVNet/data/simulations/metadata/LINE.hg38.bed");
  std::string line;
  int n_lines = 0;
  std::vector<position> records;
  while (std::getline(in_file, line)) {
    std::istringstream line_parser(line);
    position rec;
    if (!(line_parser >> rec.chr >> rec.start >> rec.stop)) { 
      std::cout << "Couldn't parse BED file line " << n_lines << std::endl;
    } 
    n_lines++;
    if (rec.stop - rec.start >= 5000) {
      //std::cout << "Parsed " << chr_name << " " << start_pos << " " << end_pos << " " << end_pos - start_pos << std::endl;  
      records.push_back(rec);
      std::cout << "Number of LINE records: " << records.size() << "\n";
    }
  }
  return records;
}


void simulate_SV_diploid(std::string ref_file, std::string parameter_file, bool coordinates, std::string output_prefix, int bases_per_snp) { //EDITED
  std::vector<position> line_records = load_line_bed_file();


  //read in list of SVs over vcf?
  //apply vcf to genome?
  srand(time(NULL));
  parameter par = parse_param(parameter_file);
  int min_chr_len = std::max(std::max(par.dup_max, par.indel_max), std::max(par.inv_max, par.translocations_max));

  //Create two copies of the reference, A and B
  std::map<std::string, std::string> genomeA = read_fasta(ref_file, min_chr_len);
  std::map<std::string, std::string> genomeB = read_fasta(ref_file, min_chr_len);

  //create bed files
  std::string outA = output_prefix;
  std::string outB = output_prefix;
  std::string outAB = output_prefix;
  outA += ".hetA.bed";
  outB += ".hetB.bed";
  outAB += ".homAB.bed";
  FILE *file2A;
  FILE *file2B;
  FILE *file2AB;
  file2A = fopen(outA.c_str(), "w");
  file2B = fopen(outB.c_str(), "w");
  file2AB = fopen(outAB.c_str(), "w");
  if (file2A == NULL || file2B == NULL || file2AB == NULL) {
    std::cout << "Error in printing: The file or path that you set " << outA.c_str() << " is not valid. It can be that there is no disc space available." << std::endl;
    exit(0);
  }

  std::cerr << "SNPS" << std::endl;
  apply_SNP(genomeA,genomeB,bases_per_snp,file2A,file2B,file2AB);
  
  std::cerr << "CHECK" << std::endl;
  check_genome(genomeA, "First:");
  check_genome(genomeB, "First:");

  std::cerr << "generate SV" << std::endl;
  std::vector<struct_var_diploid> svs; //struct_var_diploid has flags for whether it's in copy A and in copy B
  if (coordinates || !coordinates) {
    //simulate reads
    svs = generate_mutations_diploid(parameter_file, genomeA, line_records); //which genome doesn't matter, just chooses coordinates 
    check_genome(genomeA, "Sec:");
    check_genome(genomeB, "Sec:");
    apply_mutations_diploid(genomeA, genomeB, svs);

  } //else { //ignore this case
  //svs = generate_mutations_ref(parameter_file, genome);
  //check_genome(genome, "Sec:");
  //apply_mutations_ref(genome, svs);	//problem: We need two different coordinates. Simulate once for one and then for the other???
  //}

  check_genome(genomeA, "Last, A:");
  check_genome(genomeB, "Last, B:");
  std::cout << "write genome" << std::endl;
  write_fasta_diploid(output_prefix, genomeA,genomeB);
  std::cout << "write SV" << std::endl;
  write_sv_diploid(output_prefix, svs,file2A,file2B,file2AB);
}
