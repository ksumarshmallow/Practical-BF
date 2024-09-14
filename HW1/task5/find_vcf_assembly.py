import subprocess
import re
from dataclasses import dataclass

@dataclass
class VCFAssemblyFinder:
    """Класс для поиска сборки в .vcf.gz файле."""
    vcf_file_path: str

    def _get_vcf_header(self):
        """Выполняет команду для извлечения заголовка VCF файла."""
        cmd = ["bcftools", "view", "-h", self.vcf_file_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout

    def find_assembly(self):
        """Ищет информацию о сборке в заголовке VCF файла."""
        header = self._get_vcf_header()
        match = re.search(r'assembly=([^>]+)', header)
        
        if match:
            assembly_info = match.group(1)
            print(f"Найдена сборка (по тегу assembly): {assembly_info}")
            return assembly_info
        
        print("Сборка (по тегу assembly) не найдена")
        return None
    
    def find_reference(self):
        """Ищет информацию о ссылке на сборку в заголовке VCF файла."""
        header = self._get_vcf_header()
        match = re.search(r'^##reference=([^>]+)', header, re.MULTILINE)
        
        if match:
            reference_info = match.group(1)
            print(f"Найдено reference (по тегу ##reference): {reference_info}")
            return reference_info
        
        print("Ссылка на сборку (по тегу ##reference) не найдена")
        return None

    def find_assembly_info(self, return_info=False):
        """Ищет информацию о сборке, сначала пытается найти сборку по тегу assembly, затем по тегу ##reference."""
        result = self.find_assembly()
        
        if result is None:
            result = self.find_reference()
        
        if return_info:
            return result