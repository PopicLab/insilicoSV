# Additional Utilities

### Additional utils
One of the features of `utils.py` is a function `get_original_base_pos(query_vcf, query_loc, query_chrom)` which takes in a
query position and vcf describing the SVs populating a given mutated genome (e.g., a vcf of the form that is output by insilicoSV
after simulation) and calculates the location of the corresponding base in the original reference. This function can be called
directly from the command line with a call of the form:
