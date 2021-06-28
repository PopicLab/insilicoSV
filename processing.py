import yaml
import sys
import tracemalloc

class Error(Exception):
    pass

class Formater():
    def __init__(self):
        self.bedpe_counter = 1

    def fasta_to_dict(self, filename):
        '''
        filename: str, location of fasta file
        -> returns a dictionary mapping ID to genome sequence
        '''

        gen_dict = dict()

        with open(filename, "r") as fin:
            lines = fin.read().strip().splitlines()
            assert(len(lines) % 2 == 0)

            for x in range(len(lines)//2):
                gen_dict[lines[2*x]] = lines[2*x + 1]    # Key = ID, value = genome sequence
        
        return gen_dict

    def dict_to_fasta(self, filename, gen_dict):
        '''
        filename: str, location of fasta file to write genome to
        gen_dict: (ID, genome)
        -> None
        '''
        pass

    def yaml_to_var_list(self, filename):
        par_list = yaml.full_load(open(filename))
        variant_list = []
        for variant in par_list:
            var_ranges = []
            if isinstance(variant["min_length"], list):
                var_ranges = [(variant["min_length"][x], variant["max_length"][x]) for x in range(len(variant["min_length"]))]
            else:
                var_ranges = [(variant["min_length"], variant["max_length"])]
            variant_list.append([variant["type"], variant["number"], var_ranges])
        return variant_list
    
    def export_to_bedpe(self, svs, ishomogeneous, id, bedfile, reset_file = True):
        if reset_file:
            FastaFile.reset_file(bedfile)
        
        with open(bedfile, "a") as fout:
            for x in range(len(svs)):
                sv = svs[x]
                place = sv.start
                for y in range(len(sv.lengths)):
                    label = sv.piece_des[y]
                    row = [str(id),
                            str(place),
                            str(place + 1),
                            str(id),
                            str(place + sv.lengths[y]),
                            str(place + sv.lengths[y] + 1),
                            label,
                            str(sv.lengths[y]),
                            str(int(ishomogeneous[x])) + "/1",
                            sv.type.name,
                            str(self.bedpe_counter)]

                    fout.write("\t".join(row) + "\n")
                    place += sv.lengths[y]
                self.bedpe_counter += 1



    def export_to_bed12(self, svs, id, fasta_ref, fasta_out, bedfile, reset_file = True):
        if reset_file:
            with open(bedfile, "w") as f_reset:
                f_reset.truncate()

        # bed file: chr_id, start, end, name, score, strand, thickStart, thickEnd, itemRGB, blockCount, blockSizes, blockStarts
        with open(bedfile,"a") as fout:
            for sv in svs:
                block_sizes = ""
                for length in sv.lengths:
                    block_sizes += str(length) + ","
                '''
                row1_list = ["Chr_{}_{}".format(str(fasta_ref), str(id)),
                            str(int(sv.start)),
                            str(sv.end),
                            sv.type.name,
                            '0',
                            '+',
                            '0',
                            '1',
                            '(255,0,0)',
                            str(len(sv.block_starts)),
                            block_sizes,
                            sv.block_starts + ","]'''
                row1_list = ["Chr_{}_{}".format(str(fasta_ref), str(id)),
                            str(int(sv.start)),
                            str(sv.end),
                            sv.type.name,
                            str(len(sv.lengths)),
                            block_sizes,
                            sv.block_starts + ","]
                
                target_block_sizes = ""
                for target_length in sv.target_lengths:
                    target_block_sizes += str(target_length) + ","
                '''
                row2_list = ["Chr_{}_{}".format(str(fasta_out), str(id)),
                            str(int(sv.target_start)),
                            str(sv.target_end),
                            sv.type.name,
                            '0',
                            '+',
                            '0',
                            '1',
                            '(255,0,0)',
                            str(len(sv.target_block_starts)),
                            target_block_sizes,
                            sv.target_block_starts + ","]'''
                
                row2_list = ["Chr_{}_{}".format(str(fasta_out), str(id)),
                            str(int(sv.target_start)),
                            str(sv.target_end),
                            sv.type.name,
                            str(len(sv.target_lengths)),
                            target_block_sizes,
                            sv.target_block_starts + ","]
                
                #print("row1: ", row1_list)
                #print("row2: ", row2_list)
                
                #fout.write('\t'.join(row1_list) + "\n")
                fout.write("\t".join(row2_list) + "\n")



class FastaFile():
    def __init__(self, filein):
        self.len_dict = dict()     # maps chromosome ID to length of chromosome
        self.order_ids = []
        count = 0
        with open(filein,"r") as fin:
            store = True
            store_var = ""
            char_count = 0
            while True:
                count += 1
                if count % 10**6 == 0:
                    print("{} Count achieved".format(count))
                c = fin.read(1)
                #print("c: ", c)
                if not c and len(store_var) == 0:    # we have reached the end of the file
                    break
                
                elif (c == ">" or not c) and len(store_var) > 0:    # when we have reached the end of a chromosome and entered a new ID
                    store = True
                    self.len_dict[store_var.strip()[1:]] = char_count
                    self.order_ids.append(store_var.strip()[1:])
                    char_count = 0
                    store_var = c
                
                elif (c == "\n" or (not c)) and store:   # when we have reached end of line with ID - we assume that the ID will be separated from the chromosome with a newline
                    store = False
                
                elif store:
                    store_var += c
                
                elif not store and c != "\n" and c != " ":
                    char_count += 1

        self.pointer = -1
        try:
            self.f = open(filein, "r")
            self.file_in = filein
            self.f.readline()

            self.fin_export1 = open(filein, "r")   # needed to read original genome when exporting
            self.fin_export2 = open(filein, "r")   # needed to read original genome when exporting
        except FileNotFoundError:
            print("Error: File \"{}\" Not Found".format(filein))
        self.fout_export = False
        #print("File handler taking up {} bytes of storage".format(sys.getsizeof(self.f)))

    def __str__(self):
        return "Len_dict: {}".format(self.len_dict)

    def fetch(self,start,end):
        # including start position, EXCLUDING end position

        if start < self.pointer:
            raise LookupError("Requested Index no longer can be accessed")

        piece = ""
        while self.pointer < end - 1:
            c = self.f.read(1)
            if c == "\n" or c == " ":
                continue
            self.pointer += 1
            if c == ">":
                raise IndexError("Index {} out of range".format(self.pointer))
            if self.pointer >= start and self.pointer < end:
                piece += c
        
        #print("Start: {}, End: {}, Piece: {}".format(start, end, piece.upper()))
        return piece.upper()
    
    def next(self):
        self.pointer = -1
        while True:
            c = self.f.read(1)
            if c == ">":
                break
            elif not c:
                return False
        self.f.readline()                # to move past the ID line
        return True

    def export_piece(self, variants, fasta_out, index, verbose = False):
        '''
        only appends to file fasta_out the edited chromosomes given in variants

        index: 0 or 1, indicates which file handlers to use
        variants: dictionary (key = ID, value = list of list [start position, end position, edited piece]})
        reset_file: boolean to indicate whether to reset the fasta_out file at the beginning
        '''

        assert (index == 0 or index == 1)
        if index == 0:
            fin_export = self.fin_export1
        elif index == 1:
            fin_export = self.fin_export2

        self.fout_export = open(fasta_out, "a")

        for id in variants:
            if id not in self.len_dict:
                raise KeyError("ID {} not found in inputted fasta file".format(id))
        
        for id in variants:
            if verbose:
                print("New ID: ", id)
            fin_export.readline()
            self.fout_export.write(">" + str(id) + "\n")
            id_variants = variants[id]

            pos = -1
            id_variants.sort()
            c = None
            for variant in id_variants:
                if verbose:
                    print("Current Variant: ", variant)
                
                while pos < variant[1] - 1 or (variant[0] == variant[1] and pos < variant[1]):
                    c = fin_export.read(1).upper()

                    if c == ">" or not c:
                        raise IndexError("Position {} out of range for file {}".format(pos, fasta_out))
                    
                    if c.isalpha():
                        pos += 1
                        if pos == variant[0]:

                            if variant[0] == variant[1]:   # handles special case for insertion
                                self.fout_export.write(c)
                                print("Writing {} at position {}".format(c, pos))

                            self.fout_export.write(variant[2])
                            if verbose:
                                print("Writing2 {} at position {}".format(variant[2], pos))
                        elif pos < variant[0]:
                            self.fout_export.write(c)
                            if verbose:
                                print("Writing3 {} at position {}".format(c, pos))
                
            if verbose:
                print("Finishing the rest of lines")
                print("============================")
            skip = True
            while c != ">" and c:    # you use <= to catch the newline after each chromosome
                if not skip and c.isalpha():
                    self.fout_export.write(c)
                    if verbose:
                        print("Writing {} at position {}".format(c,pos))
                skip = False

                c = fin_export.read(1).upper()
                if c.isalpha():
                    pos += 1
            assert(pos == self.len_dict[id] - 1)
            self.fout_export.write("\n")
        
        self.fout_export.close()

    @staticmethod
    def reset_file(filename):
        print("Overwritting File {}...".format(filename))
        with open(filename, "w") as f_reset:
            f_reset.truncate()

    def close(self):
        self.f.close()
        self.fin_export.close()


if __name__ == "__main__":
    
    tracemalloc.start()

    x = 10
    y = 15
    filein = "reference/test.fna"
    fileout1 = "reference/test1_out.fna"
    fileout2 = "reference/test2_out.fna"
    fasta = FastaFile(filein)
    print(fasta)
    print(fasta.fetch(x,y))

    print(fasta.fetch(x + 10,y + 10))
    #fasta.next()
    #print(fasta.fetch(80,85))
    variants = {"Chromosome19": [[2,5,"TTTTTTTT"],[7,10, "AAAAAAAAAAAAAAA"], [15, 20, "TC"]]}
    variants2 = {"Chromosome19": [[2,5,"UUUUUU"]], "Chromosome21": [[4,225,"AAA"]]}
    FastaFile.reset_file(fileout1)
    FastaFile.reset_file(fileout2)
    fasta.export_piece(variants, fileout1, 0, verbose = False)
    fasta.export_piece(variants2, fileout2, 1, verbose = False)
    #fasta.close()

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()

    '''tracemalloc.start()
    with open("test.txt", "r") as fin:
        a = fin.readline()
        fin.readline()

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()'''


