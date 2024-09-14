import os
import mygene
from subprocess import PIPE
from dataclasses import dataclass

from utils.run_cmd import run_cmd

@dataclass
class GeneCoverageAnalyzer:
    """Класс для поиска доли гена в .bam файле с покрытием больше порогового"""
    bam_file_path: str
    gene: str   # can be NCBI/RefSeq id or gene name, or gene symbol
    coverage_threshold: float = 10
    genome_assembly: str = "hg19"  # or hg38
    species: str = "human"

    def __post_init__(self):
        self.data_folder = os.path.dirname(self.bam_file_path)
        self.gene_coord_file_path = os.path.join(self.data_folder, f"{self.gene}_{self.species}.bed")

        stem = os.path.basename(self.bam_file_path).split('.')[0]
        self.filtered_bed_file_path = os.path.join(self.data_folder, f"{stem}_x{self.coverage_threshold}.bed")
        self.intersected_file_path = os.path.join(self.data_folder, f"{stem}_{self.gene}_{self.species}.bed")
    
    def compute(self, return_info=False):
        """Главный метод для поиска доли гена в .bam файле с покрытием больше порогового"""

        # Получаем координаты гена и сохраняем их в .bed
        if not os.path.isfile(self.gene_coord_file_path):
            gene_pos = self.fetch_gene_info()
            self.save_to_bed(gene_pos=gene_pos)
        else:
            print(f"Файл с координатами гена {self.gene} уже существует по пути {self.gene_coord_file_path}")

        # Получаем .bam-файл с покрытием больше порогового
        if not os.path.isfile(self.filtered_bed_file_path):
            self.filter_bam_by_coverage()
        else:
            print(f"Файл с отфильтрованными регионами по покрытию уже существует по пути {self.filtered_bed_file_path}")

        # Запускаем команду для пересечения .bed-файла с .bam-файлом
        if not os.path.isfile(self.intersected_file_path):
            self.intersect()
        else:
            print(f"Файл пересечения уже существует по пути {self.intersected_file_path}")

        # Ищем долю гена с покрытием больше порогового
        intersected_length = self.get_bed_length(self.intersected_file_path)
        gene_length = self.get_bed_length(self.gene_coord_file_path)
        coverage_fraction = intersected_length / gene_length if gene_length != 0 else 0
        
        print(f"Длина гена {self.gene}: {gene_length} п.о.")
        print(f"Длина пересечения гена {self.gene} и .bam-файла: {intersected_length} п.о.")
        print(f"Доля гена {self.gene} с покрытием больше {self.coverage_threshold}x: {100 * coverage_fraction:.2f}%")
        
        if return_info:
            return {"intersected_length": intersected_length,
                    "gene_length": gene_length,
                    "coverage_fraction": coverage_fraction}

    def intersect(self):
        cmd = ["bedtools", "intersect", "-a", self.filtered_bed_file_path, "-b", self.gene_coord_file_path]
        result = run_cmd(cmd)

        with open(self.intersected_file_path, 'w') as intersect_file:
            intersect_file.write(result.stdout)

        print(f"Пересечение {self.filtered_bed_file_path} и {self.gene_coord_file_path} сохранено по пути {self.intersected_file_path}")

    def filter_bam_by_coverage(self):
        """Функция для фильтрации .bam файла в соответствии с покрытием больше порогового"""
        cmd_genomecov = ["bedtools", "genomecov", "-ibam", self.bam_file_path, "-bg"]
        genomecov_result = run_cmd(cmd_genomecov).stdout

        cmd_awk = f"awk '$4 >= {self.coverage_threshold}'"
        awk_result = run_cmd(cmd_awk, input_data=genomecov_result, shell=True)
        
        with open(self.filtered_bed_file_path, 'w') as bed_file:
            bed_file.write(awk_result.stdout)

        print(f"Регионы с покрытием больше x{self.coverage_threshold} сохранены по пути {self.filtered_bed_file_path}")
    
    def fetch_gene_info(self):
        """Функция для получения информации о гене (координат старта и конца) в указанной сборке"""
         # Создаем объект MyGeneInfo
        mg = mygene.MyGeneInfo()

        # Поиск гена
        gene_info = mg.query(self.gene, species=self.species, fields=f"genomic_pos_{self.genome_assembly}")

        if not gene_info['hits']:
            raise ValueError(f"Ген '{self.gene}' не найден")

        # Получаем координаты гена для указанной сборки
        gene_pos = gene_info['hits'][0].get(f'genomic_pos_{self.genome_assembly}', None)

        if not gene_pos:
            raise ValueError(f"Координаты для '{self.gene}' не найдены")

        return gene_pos
 
    def get_bed_length(self, bed_file_path):
        cmd = f"awk '{{sum += $3 - $2}} END {{print sum}}' {bed_file_path}"
        result = run_cmd(cmd, shell=True).stdout

        return int(result.strip())

    def save_to_bed(self, gene_pos):
        """
        Функция для сохранения информации о гене в формате .bed
        
        :param gene_pos: информация о гене, полученная из MyGene
        """
        # Создаем строку в формате .bed
        bed_lines = f'{gene_pos["chr"]}\t{gene_pos["start"]}\t{gene_pos["end"]}\n'
        
        # Записываем в файл
        with open(self.gene_coord_file_path, 'w') as f:
            f.write(bed_lines)
        
        print(f"Координаты гена {self.gene} сохранены в {self.gene_coord_file_path}")

