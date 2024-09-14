import os
import re
import pandas as pd
from typing import List, Dict
from dataclasses import dataclass

from utils.run_cmd import run_cmd

@dataclass
class VCFComparator:
    vcf_file_path: str

    def __post_init__(self):
        self.data_folder = os.path.dirname(self.vcf_file_path)
        self.output_file = os.path.join(self.data_folder, "pairwise_comparisons.csv")

    def run_gtcheck(self):
        cmd = f"bcftools gtcheck -E 0 {self.vcf_file_path}"
        result = run_cmd(cmd, shell=True, capture_output=True, text=True)
        return result.stdout
    
    def parse_line(self, line: str) -> Dict[str, any]:
        """Parse a single line of results."""
        parts = re.split(r'\s+', line.strip())
        return {
            "Sample1": parts[1],
            "Sample2": parts[2],
            "Discordance": float(parts[3]),
            "Average_neg_log_P_HWE": float(parts[4]),
            "Number_of_sites_compared": int(parts[5]),
            "Number_of_matching_genotypes": int(parts[6])
        }

    def parse_results(self, results):
        lines = results.strip().split('\n')
        parsed_results = [
            self.parse_line(line) for line in lines 
            if line.startswith("DCv2") and len(re.split(r'\s+', line.strip())) == 7
            ]
        return parsed_results

    def save_results_to_csv(self, parsed_df):
        parsed_df.to_csv(self.output_file, index=False)

        print(f"Результат сравнения образцов из .vcf записан в {self.output_file}")

    def compare(self, save_result: bool = False):
        result = self.run_gtcheck()
        parsed_results = self.parse_results(result)
        parsed_df = pd.DataFrame(parsed_results)
        if save_result:
            self.save_results_to_csv(parsed_df)
        return parsed_df