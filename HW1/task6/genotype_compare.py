import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from dataclasses import dataclass

sns.set_style("darkgrid", {"axes.facecolor": ".9"})
warnings.simplefilter('ignore')

from utils.run_cmd import run_cmd

@dataclass
class GenotypeComparator:
    vcf_file_path: str

    def __post_init__(self):
        self.data_folder = os.path.dirname(self.vcf_file_path)
        self.output_file = os.path.join(self.data_folder, "sample_genotypes.csv")
        self.gt_columns = ["POS", "SAMPLE", "REF", "ALT", "GT", "QUAL"]

    def run_gt_calc(self):
        """
        Запускает команду bcftools для извлечения генотипов из VCF файла и сохраняет их в CSV файл.
        """
        cmd = f"bcftools query -f '[%POS %SAMPLE %REF %ALT %GT %QUAL\n]' {self.vcf_file_path} > {self.output_file}"
        run_cmd(cmd, shell=True, capture_output=True, text=True)
        print(f"Генотипы из {self.vcf_file_path} записаны в {self.output_file}")

    def get_genotypes(self):
        """
        Загружает данные о генотипах из CSV файла или извлекает их из VCF файла, если CSV файл не существует.
        
        Returns:
            pd.DataFrame: DataFrame с данными о генотипах.
        """
        if not os.path.isfile(self.output_file):
            self.run_gt_calc()

        self.gt_data = pd.read_csv(self.output_file, sep=" ", header=None, names=self.gt_columns)
        return self.gt_data
    
    def preprocessing(self):
        """
        Подготавливает данные для анализа, создавая внутренние индексы для позиций и генотипов.
        """
        if not hasattr(self, 'gt_data'):
            self.get_genotypes()
        
        # 1. Для каждой позиции делаем внутренние индексы
        self.positions = self.gt_data.POS.unique()
        self.positions_idx = np.arange(len(self.positions))

        # словарь геномный индекс : внутренний индекс
        self.pos2idx = dict(zip(self.positions, self.positions_idx))

        # 2. Делаем label encoding для генотипов
        # Сначала заменим 1|0 на 0|1
        self.gt_data.GT = self.gt_data.GT.replace("1|0", "0|1")
        self.genotypes = self.gt_data.GT.unique()
        self.genotypes_idx = np.arange(len(self.genotypes)) 

        # словарь генотип : label
        self.gt2idx = dict(zip(self.genotypes, self.genotypes_idx))
    
    def process_sample(self, sample):
        """
        Обрабатывает данные для одного образца и создает матрицу вероятностей.
        
        Args:
            sample (str): Идентификатор образца.
        
        Returns:
            np.ndarray: Матрица вероятностей размером (n_positions, n_genotypes).
        """
        # Фильтруем данные для текущего образца
        data_sample = self.gt_data[self.gt_data.SAMPLE == sample]
        
        # Преобразуем позиции и генотипы в индексы
        data_sample["POS_IDX"] = data_sample["POS"].map(self.pos2idx)
        data_sample["GT_IDX"] = data_sample["GT"].map(self.gt2idx)
        
        # Подсчитываем количество наблюдений
        tm = data_sample[["POS_IDX", "GT_IDX"]].value_counts()
        pos_idx_vector = tm.index.get_level_values(0).values
        gt_idx_vector = tm.index.get_level_values(1).values
        counts_vector = tm.values

        # Создаем матрицу частот, [n_positions, n_genotypes]
        sample_matrix = np.zeros((len(self.positions_idx), len(self.genotypes_idx)))
        sample_matrix[pos_idx_vector, gt_idx_vector] = counts_vector
        
        # Рассчитываем вероятности
        row_sums = sample_matrix.sum(axis=1, keepdims=True)                             # [n_positions, ]
        probability_matrix = np.divide(sample_matrix, row_sums, where=row_sums != 0)    # [n_positions, n_genotypes]
        
        return probability_matrix
    
    def calculate_genotype_frequency(self):
        """
        Рассчитывает частоты генотипов для всех образцов и сохраняет их в словаре.
        
        Returns:
            dict: Словарь, где ключи — имена образцов, а значения — матрицы вероятностей.
        """
        self.preprocessing()
        self.probability_matrices = {}      # словарь {sample : probability_matrix}
        for sample in tqdm(self.gt_data.SAMPLE.unique(), desc="Обработка образцов..."):
            self.probability_matrices[sample] = self.process_sample(sample)
        
        return self.probability_matrices
    
    def compare(self):
        """
        Сравнивает частоты генотипов между всеми парами образцов.
        
        Returns:
            dict: Словарь, где ключи — пары имен образцов, а значения — матрицы сравнения.
        """
        if not hasattr(self, 'probability_matrices'):
            self.calculate_genotype_frequency()
        
        pairs = list(combinations(self.gt_data.SAMPLE.unique(), 2))
        self.paired_matrices = {f'{pair[0]}_{pair[1]}': None for pair in pairs}
        
        for sample1, sample2 in tqdm(pairs, desc="Сравнение образцов..."):
            matrix1 = self.probability_matrices[sample1]
            matrix2 = self.probability_matrices[sample2]
            itog_matrix = matrix1.T @ matrix2
            itog_matrix = pd.DataFrame(itog_matrix, index=self.genotypes, columns=self.genotypes)
            self.paired_matrices[f'{sample1}_{sample2}'] = itog_matrix
        
        return self.paired_matrices

    def plot_sample_heatmap(self, matrix, ax, title):
        """
        Строит тепловую карту для матрицы вероятностей.
        
        Args:
            matrix (np.ndarray): Матрица вероятностей для построения.
            ax (matplotlib.axes.Axes): Объект Axes для построения графика.
            title (str): Заголовок тепловой карты.
        """
        sns.heatmap(matrix, annot=True, fmt=".2%", cmap="YlGnBu", ax=ax, cbar=False, linewidths=0.5)
        ax.set_title(title)
        ax.set_xlabel("Генотип")
        ax.set_ylabel("Генотип")

    def plot_heatmaps(self):
        """
        Строит тепловые карты для всех пар сравнений образцов.
        """
        if not hasattr(self, 'paired_matrices'):
            self.compare()
        
        n = len(self.paired_matrices)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), sharex=True, sharey=True)
        
        if n == 1:
            axes = [axes]
        
        for (pair, matrix), ax in zip(self.paired_matrices.items(), axes):
            matrix_normalized = matrix / matrix.sum().sum()
            self.plot_sample_heatmap(matrix_normalized, ax, pair)
        
        plt.tight_layout()
        plt.show()