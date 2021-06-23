import yaml
import sys
import tracemalloc

class Formater():
    def __init__(self):
        pass

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

    def yaml_to_dict(self, filename):
        par_list = yaml.full_load(open(filename))
        variant_list = []
        for variant in par_list:
            var_ranges = []
            if isinstance(variant["min_length"], list):
                var_ranges = [(variant["min_length"][x], variant["max_length"][x]) for x in range(len(variant["min_length"]))]
            else:
                var_ranges = [(variant["min_length"], variant["max_length"])]
            variant_list.append([variant["type"], var_ranges])
        return variant_list
        #return filename

class FastaFile():
    def __init__(self, filein):
        self.len_dict = dict()     # maps chromosome ID to length of chromosome
        with open(filein,"r") as fin:
            store = True
            store_var = ""
            char_count = 0
            while True:
                c = fin.read(1)
                if not c and len(store_var) == 0:    # we have reached the end of the file
                    break
                elif c == "\n" or (not c):           # second condition accounts for the last chromosome
                    if store == False:
                        self.len_dict[store_var] = char_count
                        char_count = 0
                        store_var = ""

                    store = not store
                elif store:
                    store_var += c
                elif not store:
                    char_count += 1

        self.pointer = -1
        try:
            self.f = open(filein, "r")
            self.file_in = filein
            self.f.readline()

            self.fin_export = open(filein, "r")   # needed to read original genome when exporting
        except FileNotFoundError:
            print("Error: File \"{}\" Not Found".format(filein))
        self.fout_export = False
        #print("File handler taking up {} bytes of storage".format(sys.getsizeof(self.f)))


    def fetch(self,start,end):
        if start < self.pointer:
            raise LookupError("Requested Index no longer can be accessed")

        piece = ""
        while self.pointer < end:
            c = self.f.read(1)
            self.pointer += 1
            if c == "\n":
                raise IndexError("Ending index out of range")
            if self.pointer >= start and self.pointer < end:
                piece += c
        return piece
    
    def next(self):
        self.pointer = -1
        while True:
            c = self.f.read(1)
            if c == "\n":
                break
            elif not c:
                return False
        self.f.readline()                # to move past the ID line
        return True

    def export_piece(self, variants, fasta_out, bed_out = None, reset_file = True):
        '''
        variants: dictionary (key = ID, value = list of list [start position, end position, edited piece]})
        reset_file: boolean to indicate whether to reset the fasta_out file at the beginning
        '''

        if reset_file:
            self.reset_file(fasta_out)

        self.fout_export = open(fasta_out, "a")

        for id in variants:
            if id not in self.len_dict:
                raise KeyError("ID {} not found in inputted fasta file".format(id))
        
        for id in variants:
            #print("New ID: ", id)
            self.fout_export.write(str(id) + "\n")
            self.fin_export.readline()
            id_variants = variants[id]

            pos = -1
            id_variants.sort()
            for variant in id_variants:
                #print("Current Variant: ", variant)

                while pos < variant[1]-1:
                    c = self.fin_export.read(1)
                    pos += 1
                    #print(pos)
                    if pos == variant[0]:
                        self.fout_export.write(variant[2])
                        #print("Writing: ", variant[2])
                    elif pos < variant[0]:
                        self.fout_export.write(c)
                        #print("Writing: ", c)
                        if c == "\n" or not c:
                            raise IndexError("Position {} out of range".format(pos))

            while pos <= self.len_dict[id] and c != "\n" and c:    # you use <= to catch the newline after each chromosome
                #print(pos)
                c = self.fin_export.read(1)
                self.fout_export.write(c)
                #print("Writing: ", c)
                pos += 1
                
            #print("Len_dict: ", self.len_dict)
            #print(pos)
            #print(fin.readline())
            #print("C left off: ", c)
        
        self.fout_export.close()
        

    def reset_file(self, filename):
        print("Resetting File...")
        with open(filename, "w") as f_reset:
            f_reset.truncate()

    def close(self):
        self.f.close()
        self.fin_export.close()


if __name__ == "__main__":
    
    '''tracemalloc.start()

    x = 3
    y = 5
    fasta = FastaFile("test.txt")
    variants = {"1": [[2,5,"TTTTTTTT"],[7,10, "AAAAAAAAAAAAAAA"]]}
    variants2 = {"2": [[4,10,"AAA"]]}
    fasta.export_piece(variants, "test_out.txt")
    fasta.export_piece(variants2, "test_out.txt")
    #fasta.close()
    #fasta.reset_file("test_out.txt")

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()'''

    '''tracemalloc.start()
    with open("test.txt", "r") as fin:
        a = fin.readline()
        fin.readline()

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()'''


