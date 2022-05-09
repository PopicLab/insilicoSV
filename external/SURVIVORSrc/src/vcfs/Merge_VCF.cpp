/*
 * Merge_VCF.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: fsedlaze
 */

#include "Merge_VCF.h"

//read in all the vcf filenames:
std::vector<std::string> parse_filename(std::string filename) {
	std::vector<std::string> names;
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "File Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		names.push_back(std::string(buffer));
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();

	return names;
}

strcoordinate parse_stop(const char * buffer) {
	size_t i = 0;
	bool chr_flag = false;
	strcoordinate pos;
	pos.chr = "";
	pos.pos = -1;
	while (buffer[i] != '\t' && (buffer[i] != '\n' && buffer[i] != '\0')) {
		if (strncmp(&buffer[i], ";END=", 5) == 0) {
			pos.pos = atoi(&buffer[i + 5]);
		}
		if ((strncmp(&buffer[i], "END=", 4) == 0 && i == 0)) {
			pos.pos = atoi(&buffer[i + 4]);
		}
		if (strncmp(&buffer[i], "CHR2=", 5) == 0) {
			i = i + 5;
			chr_flag = true;
		}
		if (buffer[i] == ';') {
			chr_flag = false;
		}
		if (chr_flag) {
			pos.chr += buffer[i];
		}

		i++;
	}
	//std::cout<<"end: "<<pos.chr<<" "<<pos.pos<<std::endl;
	return pos;
}
std::pair<bool, bool> parse_strands(const char * buffer) {
	std::pair<bool, bool> strands;
	size_t i = 0;
	while (buffer[i] != '\t' && (buffer[i] != '\n' && buffer[i] != '\0')) {
		if (strncmp(&buffer[i], "3to5", 4) == 0) {
			strands.first = false;
			strands.second = true;
			return strands;
		}
		if (strncmp(&buffer[i], "3to3", 4) == 0) {
			strands.first = false;
			strands.second = false;
			return strands;
		}
		if (strncmp(&buffer[i], "5to3", 4) == 0) {
			strands.first = true;
			strands.second = false;
			return strands;
		}
		if (strncmp(&buffer[i], "5to5", 4) == 0) {
			strands.first = true;
			strands.second = true;
			return strands;
		}
		i++;
	}
	return strands;
}
std::vector<int> parse_callers(char* buffer) {
	size_t i = 0;
	std::vector<int> entries;
	//std::cout<<buffer[i]<<std::endl;
	entries.push_back(atoi(&buffer[i]));
	while (buffer[i] != ';' && buffer[i] != '\0') {
		if (buffer[i] == ',') {
			entries.push_back(atoi(&buffer[i + 1]));
		}
		i++;
	}
	//std::cout<<entries.size()<<std::endl;
	return entries;
}

short get_type(std::string type) {
	if (strncmp(type.c_str(), "DEL", 3) == 0) {
		return 0;
	} else if (strncmp(type.c_str(), "DUP", 3) == 0) {
		return 1;
	} else if (strncmp(type.c_str(), "INV", 3) == 0) {
		return 2;
	} else if (strncmp(type.c_str(), "TRA", 3) == 0) {
		return 3;
	} else if (strncmp(type.c_str(), "INS", 3) == 0) {
		return 4;
	} else if (strncmp(type.c_str(), "BND", 3) == 0) { //can be inv/inter/tra
		return 5;
	} else {
		std::cerr << "Unknown type!" << type << std::endl;

	}
	return -1;
}

strcoordinate parse_pos(char * buffer) {
	std::string tmp = std::string(buffer);
	size_t found = tmp.find(':');
	strcoordinate pos;
	pos.chr = tmp.substr(0, found);
	found++;
	pos.pos = atoi(&tmp[found]);
	//std::cout << pos.chr << "| " << pos.pos << "|" << std::endl;
	return pos;
}

std::pair<int, int> parse_manta(char * buffer) {
	std::pair<int, int> res;
	res.first = 0;
	res.second = 0;
	size_t i = 0;
	while (buffer[i] != '\t') {
		i++;
	}
	//std::cout<<buffer<<std::endl;
	// 0/1:PASS:170:220,0,561:5,3:9,4
	int count = 0;
	while (buffer[i] != '\n' && buffer[i] != '\0') {
		if (count > 3) {
			if (buffer[i - 1] == ':') {
				res.first += atoi(&buffer[i]);
				//	std::cout<<"first: "<<atoi(&buffer[i])<<std::endl;
			}
			if (buffer[i - 1] == ',') {
				res.second += atoi(&buffer[i]);
				//	std::cout<<"second: "<<atoi(&buffer[i])<<std::endl;
			}
		}
		if (buffer[i] == ':') {
			count++;
		}
		i++;
	}
	return res;
}

std::pair<int, int> parse_delly(char * buffer) {

	// GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-9.02876,0,-32.6842:90:PASS:0:22219093:0:-1:14:3:0:0

	std::pair<int, int> res;
	res.first = 0;
	res.second = 0;
	size_t i = 0;
	while (buffer[i] != '\t') {
		i++;
	}
	//std::cout<<buffer<<std::endl;
	int count = 0;
	while (buffer[i] != '\n' && buffer[i] != '\0') {
		if ((count == 8 || count == 9) && buffer[i - 1] == ':') {
			res.first += atoi(&buffer[i]);
			//std::cout<<"first: "<< atoi(&buffer[i])<<std::endl;
		}
		if ((count == 10 || count == 11) && buffer[i - 1] == ':') {
			res.second += atoi(&buffer[i]);
			//	std::cout<<"second: "<< atoi(&buffer[i])<<std::endl;
		}
		if (buffer[i] == ':') {
			count++;
		}
		i++;
	}
	return res;
}

int parse_support(char * buffer) {
	size_t i = 0;
	int support = 0;
	while (buffer[i] != '\t' && buffer[i] != '\0') {

		if (strncmp(&buffer[i], "VT_AC=", 6) == 0) {
			support = atoi(&buffer[i + 6]);
		}
		if ((strncmp(&buffer[i], ";SU=", 4) == 0 || strncmp(&buffer[i], ";RE=", 4) == 0) || (strncmp(&buffer[i], ";PE=", 4) == 0 || strncmp(&buffer[i], ";SR=", 4) == 0)) { // SU: Lumpy, RE: Sniffles
			support += atoi(&buffer[i + 4]);
		}
//TOOD extned for the tags that other caller use!
		i++;
	}
	return support;
}
std::pair<bool, bool> parse_strands_lumpy(char * buffer) {
	std::pair<bool, bool> strand;
	size_t i = 0;
	bool is_first = true;
	while (buffer[i] != '\t') {
		if (buffer[i] == '[') {
			if (is_first) {
				strand.first = false;
				is_first = false;
			} else {
				strand.second = false;
			}
		} else if (buffer[i] == ']') {
			if (is_first) {
				strand.first = true;
				is_first = false;
			} else {
				strand.second = true;
			}
		}
		i++;
	}
	return strand;
}
std::string get_most_effect(std::string alt, int ref) {
	size_t i = 0;
	std::string most_alt = "";

	std::string tmp = "";
	while (i < alt.size()) {
		if (alt[i] == ',') {
			int size = most_alt.size();
			int curr = tmp.size();
			if (abs(ref - curr) > abs(ref - size)) {
				most_alt = tmp;
			}
			tmp.clear();
		} else {
			tmp += alt[i];
		}
		i++;
	}
	return most_alt;
}

//for each file parse the entries
std::vector<strvcfentry> parse_vcf(std::string filename) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(filename.c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "VCF Parser: could not open file: " << filename.c_str() << std::endl;
		exit(0);
	}

	std::vector<strvcfentry> calls;
	myfile.getline(buffer, buffer_size);

	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			strvcfentry tmp;
			tmp.sup_lumpy = 0;
			tmp.stop.pos = -1;
			tmp.type = -1;
			bool set_strand = false;
			std::string ref;
			std::string alt;
			tmp.genotype = "./.";
			tmp.strands.first = true;
			tmp.strands.second = true;
			tmp.num_reads.first = 0;
			tmp.num_reads.second = 0;
			//std::cout<<buffer<<std::endl;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {

				if (count == 0 && buffer[i] != '\t') {
					tmp.start.chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {
					tmp.start.pos = atoi(&buffer[i]);
					//std::cout<<tmp.start.pos<<std::endl;
				}
				if (count == 3 && buffer[i] != '\t') {
					ref += buffer[i];
				}
				if (count == 4 && buffer[i] != '\t') {
					alt += buffer[i];
				}
				if (count == 4 && buffer[i - 1] == '\t') {
					tmp.strands = parse_strands_lumpy(&buffer[i]);
				}
				if (tmp.stop.pos == -1 && (count == 7 && buffer[i - 1] == '\t')) {
					tmp.stop = parse_stop(&buffer[i]);
					//std::cout<<tmp.stop.chr<<std::endl;
				}
				if (count == 7 && strncmp(&buffer[i], "SVTYPE=", 7) == 0) {
					tmp.type = get_type(std::string(&buffer[i + 7]));
				}
				if (count == 7 && strncmp(&buffer[i], ";SU=", 4) == 0) { //for lumpy!
					tmp.num_reads.second = atoi(&buffer[i + 4]);
				}
				if (count == 7 && strncmp(&buffer[i], ";CT=", 4) == 0) {
					//parse strand delly:
					set_strand = true;
					tmp.strands.first = (bool) (buffer[i + 4] != '5');
					tmp.strands.second = (bool) (buffer[i + 7] != '5');
				}
				if (count == 7 && strncmp(&buffer[i], ";STRANDS=", 9) == 0) {
					set_strand = true;
					tmp.strands.first = (bool) (buffer[i + 9] == '+');
					tmp.strands.second = (bool) (buffer[i + 10] == '+');
				}

				if (count == 9 && buffer[i - 1] == '\t') { //parsing genotype;
					size_t j = i;
					tmp.genotype = "";
					while (buffer[j] != '\0' && buffer[j] != ':') {
						tmp.genotype += buffer[j];
						j++;
					}
					//	std::cout<<"GO: "<<tmp.genotype<<std::endl;
				}
				if (count == 8 && strncmp(&buffer[i], "PR:SR", 5) == 0) {
					//manta
					tmp.num_reads = parse_manta(&buffer[i]);
					//std::cout<<"HIT MANTA "<<tmp.num_reads.first<<" "<<tmp.num_reads.second<<std::endl;
				}
				if (count == 8 && strncmp(&buffer[i], "DR:DV:RR:RV", 11) == 0) {
					//delly
					tmp.num_reads = parse_delly(&buffer[i]);
				}

				if (count == 4 && buffer[i - 1] == '<') {
					tmp.type = get_type(std::string(&buffer[i]));
				}
				if (tmp.stop.pos == -1 && (count == 4 && (buffer[i - 1] == '[' || buffer[i - 1] == ']'))) {
					tmp.stop = parse_pos(&buffer[i]);
				}

				if (count == 9 && buffer[i - 1] == '\t') {
					tmp.calls[filename] = std::string(&buffer[i]);
					break;
				}

				if (count < 9) {
					tmp.header += buffer[i];
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			if (!set_strand) {
				if (tmp.type == 0 || tmp.type==4) {
					tmp.strands.first = true;
					tmp.strands.second = false;
				} else if (tmp.type == 1) {
					tmp.strands.first = false;
					tmp.strands.second = true;
				}else { //should not happen??
					tmp.strands.first = true;
					tmp.strands.second = true;
				}
			}
			if (tmp.stop.pos == -1) {

				std::size_t found = alt.find(",");
				if (found != std::string::npos) {
					alt = get_most_effect(alt, (int) ref.size());
				}
				tmp.stop.chr = tmp.start.chr;
				int len = (int) ref.size() - (int) alt.size();
				tmp.stop.pos = tmp.start.pos + abs(len);
				if (len > 0) {
					tmp.type = 0;
				} else if (len < 0) {
					tmp.type = 1;
				}
			}
			if (tmp.stop.chr.empty()) {
				tmp.stop.chr = tmp.start.chr;
			}
			calls.push_back(tmp);
			tmp.calls.clear();
		} else {

		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
//std::cout << calls.size() << std::endl;
	return calls;
}

int overlap(strvcfentry tmp, std::vector<strvcfentry> & final_vcf, int max_dist) {
	for (size_t i = 0; i < final_vcf.size(); i++) {
		//check type:
		if (final_vcf[i].type == tmp.type) {
			//check chrs:
			if (strcmp(final_vcf[i].stop.chr.c_str(), tmp.stop.chr.c_str()) == 0 && strcmp(final_vcf[i].start.chr.c_str(), tmp.start.chr.c_str()) == 0) {
				//check coordinates:
				if (abs(final_vcf[i].stop.pos - tmp.stop.pos) < max_dist && abs(final_vcf[i].start.pos - tmp.start.pos) < max_dist) {
					return i;
				}
			}
		}
	}
	return -1;
}

//detect overlap and merge:
void merge_entries(std::string filename, int max_dist, std::vector<strvcfentry> & final_vcf) {
//get new entries
	std::vector<strvcfentry> new_entries = parse_vcf(filename);
//merge entires:
	for (size_t i = 0; i < new_entries.size(); i++) {
		int id = overlap(new_entries[i], final_vcf, max_dist);
		if (id > -1) {
			//std::cout<<"add " <<new_entries[i].calls[filename]<<std::endl;
			final_vcf[id].calls[filename] = new_entries[i].calls[filename]; //add call to entries;
		} else {
			//std::cout<<"push " <<new_entries[i].calls[filename]<<std::endl;
			final_vcf.push_back(new_entries[i]);
		}
	}
}

std::string get_header(std::vector<std::string> names) {
	size_t buffer_size = 2000000;
	char*buffer = new char[buffer_size];
	std::ifstream myfile;

	myfile.open(names[0].c_str(), std::ifstream::in);
	if (!myfile.good()) {
		std::cout << "Annotation Parser: could not open file: " << names[0].c_str() << std::endl;
		exit(0);
	}
	std::string header;
	myfile.getline(buffer, buffer_size);
	while (!myfile.eof()) {
		if (buffer[0] == '#' && buffer[1] == '#') {
			header += std::string(buffer);
			header += '\n';
		} else if (buffer[0] == '#') {
			int count = 0;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count < 9) {
					header += buffer[i];
				}
				if (count == 9) {
					break;
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
			break;
		}
		myfile.getline(buffer, buffer_size);
	}

	for (size_t i = 0; i < names.size(); i++) {
		header += '\t';
		header += names[i];
	}
	header += '\n';
	myfile.close();
	return header;
}

void print_merged_vcf(std::string outputfile, std::string header, std::vector<strvcfentry> &final_vcf, std::vector<std::string> names) {
	FILE *file;
	file = fopen(outputfile.c_str(), "w");

	fprintf(file, "%s", header.c_str());
	header.clear();

	for (size_t i = 0; i < final_vcf.size(); i++) {
		fprintf(file, "%s", final_vcf[i].header.c_str());
		for (size_t t = 0; t < names.size(); t++) {
			fprintf(file, "%c", '\t');
			if (final_vcf[i].calls.find(names[t]) != final_vcf[i].calls.end()) { //found an entry
				fprintf(file, "%s", final_vcf[i].calls[names[t]].c_str());
			} else {
				fprintf(file, "%s", "./.:0,0.0,0.0:0:NotDetected:0:0:0:0:0");
			}
		}
		fprintf(file, "%c", '\n');
	}
}

//main:
void merge_vcf(std::string filenames, int max_dist, std::string outputfile) {

	std::vector<std::string> names = parse_filename(filenames);
	std::cout << "found in file: " << names.size() << std::endl;
	std::vector<strvcfentry> final_vcf;

	for (size_t i = 0; i < names.size(); i++) {
		merge_entries(names[i], max_dist, final_vcf);
		std::cout << "merged: " << final_vcf.size() << std::endl;
	}
	std::cout << "get header:" << std::endl;
	std::string header = get_header(names);
	std::cout << "print:" << std::endl;
	print_merged_vcf(outputfile, header, final_vcf, names);
}
