import os
import pandas as pd
import gzip
import zipfile
import subprocess
from dataclasses import dataclass

@dataclass
class VCFReader:
    """class for reading vcf file"""
    vcf_path: str
    genotypes: list = None
    nrows: int = None
     
    def get_extension(self):
        """Determines the file extension and returns it"""
        _, extension = os.path.splitext(self.vcf_path)
        return extension.lower()
     
    def read_vcf(self, chunksize: int = 10000):
        """Reads the VCF file and returns a DataFrame of variants"""
        extension = self.get_extension()
        self.n_rows = self.calculate_nrows()
        self.n_chunks = self.n_rows // chunksize
        
        if extension == ".gz" or extension == ".gzip":
            with gzip.open(self.vcf_path, "rt") as ifile:
                vcf_names = self._get_vcf_header(ifile)
                compression = 'gzip'
            ifile.close()
        
        elif extension == ".zip":
            with zipfile.ZipFile(self.vcf_path, 'r') as zfile:
                with zfile.open(zfile.namelist()[0]) as ifile:
                    vcf_names = self._get_vcf_header(ifile)
                    compression = 'zip'
                ifile.close()
            zfile.close()
            
        elif extension == ".vcf":
            with open(self.vcf_path, 'r') as ifile:
                vcf_names = self._get_vcf_header(ifile)
                compression = None
            ifile.close()

        self.genotypes = vcf_names[9:]
        
        vcf_df = pd.read_csv(self.vcf_path, 
                             names=vcf_names,
                             compression=compression,
                             comment='#', 
                             delim_whitespace=True, 
                             chunksize=chunksize)
        return vcf_df
    
    def _get_vcf_header(self, file):
        """Extracts the VCF header (column names) from the file."""
        for line in file:
            if line.startswith("#CHROM"):
                vcf_names = line.strip('#\n').split('\t')
                return vcf_names
    
    def calculate_nrows(self):
        """Returns the number of rows in the VCF file"""
        result = subprocess.run(['wc', '-l', self.vcf_path], stdout=subprocess.PIPE, text=True)
        n_rows = int(result.stdout.split()[0])
        return n_rows
    
