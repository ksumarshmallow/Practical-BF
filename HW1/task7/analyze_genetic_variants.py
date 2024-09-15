import os
import pandas as pd

from tqdm import tqdm
from itertools import combinations
import random
from dataclasses import dataclass

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
import seaborn as sns
import warnings
from utils.run_cmd import run_cmd

random.seed(17)
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
warnings.simplefilter('ignore')

@dataclass
class GeneticVariantAnalyzer:
    population_data_path: str
    frequency_data_path: str

    def __post_init__(self):
        self.data_folder = os.path.dirname(self.population_data_path)
        self.population_data = pd.read_csv(self.population_data_path, sep='\t')
        self.populations = self.population_data['pop'].unique()
        self.af_files = {}
        self.slopes = {}
        self.intercepts = {}
        self.r_squares = {}
        self._log_population_info()

    def _log_population_info(self):
        print("Количество образцов в анализируемых популяциях:")
        print(self.population_data['pop'].value_counts())
    
    def calculate_af(self, save_path: str = None) -> pd.DataFrame:
        """Главный метод для подсчета генетических частот для каждого варианта среди всех популяций"""
        col_names = ["CHROM:POS", "REF", "ALT", "AF", "Population"]
        df_frequencies = pd.DataFrame(columns=col_names)
        for pop_ in tqdm(self.populations, desc='Подсчет генетических частот...', colour='GREEN'):
            af_file = self._get_or_calculate_af_file(pop_)
            df_frequencies = self._update_frequencies(df_frequencies, af_file, pop_)

        df_frequencies[["CHROM", "POS"]] = df_frequencies["CHROM:POS"].str.split(':', expand=True).values
        df_frequencies = df_frequencies[["CHROM", "POS"] + col_names[1:]]
        df_frequencies = self._infer_types(df_frequencies)
        if save_path:
            df_frequencies.to_csv(save_path, index=False)
        return df_frequencies

    def _infer_types(self, df: pd.DataFrame) -> pd.DataFrame:
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='ignore', downcast='float')
        return df

    def _get_or_calculate_af_file(self, pop_):
        af_file = f"{self.data_folder}/AF_{pop_}.txt"
        if not (os.path.exists(af_file) and os.path.getsize(af_file) > 0):
            self._calculate_af_for_population(pop_)
        return af_file

    def _calculate_af_for_population(self, pop_):
        """Метод для подсчета Allele frequency для заданной популяции"""
        # 1. Подготавливаем файлы с id образцов для каждой популяции 
        samples_file = self._create_samples_file(pop_)
        
        # 2. Создаем VCF-файлы для каждой популяции (bcftools view) и 
        # вычисляем частоты аллелей для каждой популяции (bcftools +fill-tags)
        frequency_data_pop_path = f"{self.data_folder}/data_{pop_}_AF.vcf.gz"
        bcftools_command = f"bcftools view -S {samples_file} {self.frequency_data_path} | bcftools +fill-tags -- -t AF | bcftools view -Oz -o {frequency_data_pop_path}"
        run_cmd(bcftools_command, shell=True, check=True)

        # 3. Извлекаем частоты аллелей
        af_file = f"{self.data_folder}/AF_{pop_}.txt"
        bcftools_af_command = f"bcftools query -f '%CHROM:%POS %REF %ALT %AF\n' {frequency_data_pop_path} > {af_file}"
        run_cmd(bcftools_af_command, shell=True, check=True)
        self.af_files[pop_] = af_file

    def _create_samples_file(self, pop_):
        """Создает файлы с sample_id для популяции 'pop_'"""
        samples_file = f"{self.data_folder}/samples_{pop_}.txt"
        awk_command = f"awk '$2 == \"{pop_}\"' {self.population_data_path} | cut -f1 > {samples_file}"
        run_cmd(awk_command, shell=True)
        return samples_file
    
    def _update_frequencies(self, df_frequencies, af_file, pop_):
        with open(af_file, 'r') as f:
            data = [line.split() for line in f.readlines()]
        df = pd.DataFrame(data, columns=df_frequencies.columns[:-1])
        df['Population'] = pop_
        return pd.concat([df_frequencies, df], ignore_index=True)
    
    def _get_common_positions(self, dfs: dict, populations: list) -> list:
        """Находит общие позиции для всех популяций."""
        positions = set(dfs[populations[0]].POS.values)
        for pop_ in populations[1:]:
            positions &= set(dfs[pop_].POS.values)
        return sorted(positions)

    def _split_and_filter_dfs(self, df_frequencies, populations):
        """Разбивает единую таблицу с AF (со всеми популяциями) на несколько таблиц (одна таблица на одну популяцию)"""
        dfs = {pop: df_frequencies[df_frequencies['Population'] == pop] for pop in populations}
        positions = self._get_common_positions(dfs, populations)
        return {pop: df[df.POS.isin(positions)].drop_duplicates('POS').sort_values("POS") for pop, df in dfs.items()}

    def plot_scatter(self, df_frequencies: pd.DataFrame, populations: list = None):
        """Главный метод для построения графика рассеяния"""
        df_frequencies, populations = self._check_data(df_frequencies, populations)
        dfs = self._split_and_filter_dfs(df_frequencies, populations)
        self._plot_scatter_for_populations(dfs, populations)

    def _plot_scatter_for_populations(self, dfs, populations):
        pairs = list(combinations(populations, 2))
        fig, axs = plt.subplots(1, len(pairs), figsize=(7 * len(pairs), 8))

        if len(pairs) == 1:
            axs = [axs]

        # Предполагаем, что вы можете собрать все позиции из всех пар графиков
        positions = dfs[populations[0]]['POS'].values

        # Определение нормализации по глобальному минимуму и максимуму позиций
        norm = plt.Normalize(vmin=positions.min(), vmax=positions.max())
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])  # Для пустого массива, чтобы colorbar работал

        for i, (pop1, pop2) in enumerate(pairs):
            self._plot_single_scatter(axs[i], dfs[pop1], dfs[pop2], pop1, pop2, sm)

        # Один colorbar для всех графиков
        cbar = fig.colorbar(sm, ax=axs, orientation='vertical', pad=0.02)
        cbar.ax.tick_params(labelsize=14)  # Установка размера шрифта на colorbar
        cbar.ax.yaxis.set_ticks([])  # Убрать цифры на colorbar

        plt.figlegend(["Генетический вариант"], loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=20, frameon=False)
        plt.subplots_adjust(wspace=.2)
        plt.show()

    def _plot_single_scatter(self, ax, df1, df2, pop1, pop2, sm):
        """Строит график рассеяния для одной пары популяций"""
        af_pop1 = df1['AF'].values.astype(float)
        af_pop2 = df2['AF'].values.astype(float)
        pos = df1['POS'].values.astype(float)  # Предположим, что в df1 и df2 есть столбец POS

        scatter = ax.scatter(af_pop1, af_pop2, c=pos, cmap='viridis', edgecolor='gray', s=50, alpha=0.5, norm=sm.norm)
        ax.set_xlabel(f"AF {pop1}", fontsize=18)
        ax.set_ylabel(f"AF {pop2}", fontsize=18)

        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        return scatter

    def _check_data(self, df_frequencies, populations):
        if df_frequencies is None:
            df_frequencies = self.calculate_af()
        if populations is None:
            populations = self.populations
        return df_frequencies, populations