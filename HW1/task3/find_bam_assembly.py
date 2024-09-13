import re
from dataclasses import dataclass

from utils.run_cmd import run_cmd

@dataclass
class BamAssemblyFinder:
    """Класс для поиска сборки .bam файла"""
    bam_file_path: str

    def __post_init__(self):
        self.header = self._get_bam_header()
    
    def _get_bam_header(self):
        """Выполняет команду для извлечения заголовка BAM файла"""
        cmd = ["samtools", "view", "-H", self.bam_file_path]
        result = run_cmd(cmd).stdout
        return result

    def find_reference_genome(self):
        """Ищет URI референсного генома в заголовке BAM файла"""
        sq_lines = [line for line in self.header.splitlines() if line.startswith("@SQ")]
        
        reference_genome_path = None
        for line in sq_lines:
            match = re.search(r'UR:(\S+)', line)
            if match:
                reference_genome_path = match.group(1)
                break

        if reference_genome_path:
            print(f"Найден путь к референсному файлу: {reference_genome_path}")
            return reference_genome_path
        else:
            print(f"Путь к референсному файлу не найден")
            return None
    
    def find_program_info(self):
        """Ищет информацию о программе выравнивания (@PG) в заголовке BAM файла."""
        pg_lines = [line for line in self.header.splitlines() if line.startswith("@PG")]
        
        if pg_lines:
            print("Найдена информация о программе выравнивания (@PG):")
            for line in pg_lines:
                print(line)
            
            return pg_lines
        else:
            print("Информации о программе запуска (@PG) нет")
            return None
    
    def find_assembly(self, return_info=False):
        """Ищет сборку, сначала пытается найти референсный геном, затем программу выравнивания (@PG)"""
        result = self.find_reference_genome()
        
        if not result:
            result = self.find_program_info()
        
        if return_info:
            return result
