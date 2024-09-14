import re
from dataclasses import dataclass
from utils.run_cmd import run_cmd

@dataclass
class BamAssemblyFinder:
    """Класс для поиска сборки в .bam файле."""
    bam_file_path: str

    def __post_init__(self):
        self.header = self._get_bam_header()
    
    def _get_bam_header(self):
        """Выполняет команду для извлечения заголовка BAM файла."""
        cmd = ["samtools", "view", "-H", self.bam_file_path]
        result = run_cmd(cmd).stdout
        return result

    def find_assembly(self):
        """Ищет информацию о сборке в заголовке BAM файла."""
        sq_lines = [line for line in self.header.splitlines() if line.startswith("@SQ")]
        
        for line in sq_lines:
            match = re.search(r'AS:(\S+)', line)
            if match:
                reference_genome_path = match.group(1)
                print(f"Найдена сборка (по тегу @SQ в поле AS): {reference_genome_path}")
                return reference_genome_path
        
        print("Сборка (AS в теге @SQ) не найдена")
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
            print("Информации о программе выравнивания (@PG) нет")
            return None

    def find_assembly_info(self, return_info=False):
        """Ищет информацию о сборке, сначала пытается найти сборку в @SQ, затем информацию о программе выравнивания (@PG)."""
        result = self.find_assembly()
        
        if result is None:
            result = self.find_program_info()
        
        if return_info:
            return result