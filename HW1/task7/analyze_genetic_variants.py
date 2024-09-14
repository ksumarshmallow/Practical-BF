import os
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.stats import pearsonr

from tqdm import tqdm
from itertools import combinations
import random
from dataclasses import dataclass

import matplotlib.pyplot as plt
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
        # Это будет временная папка, так как мы из .zip-архива туда все вытащили
        # Так что результаты bcftools будут удалены после всего
        self.data_folder = os.path.dirname(self.population_data_path)
        self.population_data = pd.read_csv(self.population_data_path, sep='\t')
        self.populations = self.population_data['pop'].unique()

        print("Количество образцов в анализируемых популяциях:")
        print(self.population_data['pop'].value_counts())

        # словарь для хранения файлов с Allele Frequencies (AF)
        self.af_files = {}

    def calculate_af(self, save_path: str = None) -> pd.DataFrame:
        col_names = ['AF_indicator', 'CHR', 'Allele_Frequency', 
                     'COUNT_ALL', 'COUNT_REF', 'COUNT_ALT', 
                     'POS_REF', 'REF', 'ALT', 'POS_ALT']
        df_frequencies = pd.DataFrame(columns=col_names + ['Population'])

        with tqdm(total=len(self.populations), desc='Подсчитываем AF для каждой популяции...', unit='pop', leave=True) as pbar:
            for pop_ in self.populations:
                af_file = f"{self.data_folder}/af_{pop_}.txt"
                if os.path.exists(af_file) and os.path.getsize(af_file) > 0:
                    print(f"Частоты генетических вариантов (AF) подсчитаны для популяции {pop_} и записаны по пути {af_file}")
                else:
                    self.calculate_af_pop(pop_)

                with open(af_file, 'r') as f:
                    lines = f.readlines()
                data = [line.split() for line in lines if line.startswith('AF')]
                df = pd.DataFrame(data, columns=col_names)
                df['Population'] = pop_
                df_frequencies = pd.concat((df_frequencies, df), ignore_index=True)
                                
                pbar.update(1)
        
        df_frequencies = self.infer_types(df_frequencies)

        if save_path:
            df_frequencies.to_csv(save_path, index=False)

        return df_frequencies

    def plot_scatter(self, df_frequencies: pd.DataFrame = None, populations: list = None, show_coeffs:bool = False):        
        if df_frequencies is None:
            print("'df_frequencies' не определен. Считаем частоты генетических вариантов.")
            df_frequencies = self.calculate_af()
        
        if populations is None:
            populations = self.populations
        
        dfs = self._split_dataframe_by_population(df_frequencies, populations)
        positions = self._get_common_positions(dfs, populations)
        dfs = self._filter_and_deduplicate(dfs, positions)
        self._plot_scatterplots(dfs, populations)

    def infer_types(self, df: pd.DataFrame) -> pd.DataFrame:
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='ignore', downcast='float')
        return df

    def calculate_af_pop(self, pop_: str):
        samples_file = f"{self.data_folder}/samples_{pop_}.txt"
        af_file = f"{self.data_folder}/af_{pop_}.txt"
        
        # Создание файла с образцами для текущей популяции
        awk_command = f"awk '$2 == \"{pop_}\"' {self.population_data_path} | cut -f1 > {samples_file}"
        run_cmd(awk_command, shell=True)
            
        # Расчет частот альтернативных аллелей для текущей популяции
        bcftools_command = f"bcftools view -S {samples_file} {self.frequency_data_path} | bcftools stats - | grep '^AF' > {af_file}"
        run_cmd(bcftools_command, shell=True)

        # Пополняем словарь файлов с AF
        self.af_files[pop_] = af_file

        print(f'Частоты генетических вариантов (AF) подсчитаны для популяции {pop_} и записаны по пути {af_file}')

    def _split_dataframe_by_population(self, df: pd.DataFrame, populations: list) -> dict:
        """Разделяет DataFrame по популяциям."""
        return {pop: df[df['Population'] == pop] for pop in populations}

    def _get_common_positions(self, dfs: dict, populations: list) -> list:
        """Находит общие позиции для всех популяций."""
        positions = set(dfs[populations[0]].POS_REF.values)
        for pop_ in populations[1:]:
            positions &= set(dfs[pop_].POS_REF.values)
        return sorted(positions)

    def _filter_and_deduplicate(self, dfs: dict, positions: list) -> dict:
        """Фильтрует DataFrame, оставляя только общие позиции и удаляет дубликаты"""
        for pop_ in dfs:
            dfs[pop_] = dfs[pop_][dfs[pop_].POS_REF.isin(positions)].sort_values('POS_REF').drop_duplicates('POS_REF')
        return dfs

    def _plot_scatterplots(self, dfs: dict, populations: list):
        """Строит scatterplots для попарных комбинаций популяций."""        
        pairs = list(combinations(populations, 2))
        num_plots = pairs.__len__()
        fig, axs = plt.subplots(1, num_plots, figsize=(7 * num_plots, 8))
        
        if num_plots == 1:
            axs = [axs]  # для случая, когда есть только 2 популяции

        kwargs = {"color": "skyblue", "edgecolors": 'gray', "s": 50, "alpha": 0.7}
        
        self.correlations = {}
        self.slopes = {}
        self.intercepts = {}

        for i, (pop1, pop2) in enumerate(pairs):
            af_pop1 = dfs[pop1]['Allele_Frequency'].values.astype(float).reshape(-1, 1)
            af_pop2 = dfs[pop2]['Allele_Frequency'].values.astype(float)

            # Рассчитываем корреляцию Пирсона
            corr, _ = pearsonr(af_pop1.flatten(), af_pop2)
            self.correlations[f'{pop1}_{pop2}'] = corr.item()

             # Линейная регрессия
            predictions, r_squared, slope, intercept = self.fit_linear_regression(af_pop1, af_pop2)
            self.slopes[f'{pop1}_{pop2}'] = slope.item()
            self.intercepts[f'{pop1}_{pop2}'] = intercept.item()

            # Строим scatterplot
            sns.scatterplot(x=af_pop1.flatten(), y=af_pop2, ax=axs[i], **kwargs)

            # Добавляем линию регрессии
            axs[i].plot(af_pop1, predictions, color="salmon", lw=1.5)

            # Добавляем легенду
            self.add_legend(axs[i], r_squared)

            axs[i].set_xlabel(f"AF {pop1}", fontsize=14)
            axs[i].set_ylabel(f"AF {pop2}", fontsize=14)

        plt.tight_layout()
        plt.show()
    
    def fit_linear_regression(self, af_pop1, af_pop2):
        """Функция для вычисления линейной регрессии и метрик"""
        model = LinearRegression()
        model.fit(af_pop1, af_pop2)
        predictions = model.predict(af_pop1)
        r_squared = r2_score(af_pop2, predictions)
        slope = model.coef_[0]
        intercept = model.intercept_
        
        return predictions, r_squared, slope, intercept
    
    def add_legend(self, ax, r_squared):
        """Функция для добавления текста в легенду"""
        legend_text = [f"Генетический вариант", f"Линейная регрессия ($R^2={r_squared:.2f}$)"]
        ax.legend(
            legend_text, 
            loc='upper center', 
            bbox_to_anchor=(0.5, 1.1), 
            fontsize='large', 
            frameon=False, 
            title=None, 
            title_fontsize='large',
            ncol=1
        )