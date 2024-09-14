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
        df_frequencies = pd.DataFrame(columns=[
            'AF_indicator', 'CHR', 'Allele_Frequency', 'COUNT_ALL', 
            'COUNT_REF', 'COUNT_ALT', 'POS_REF', 'REF', 'ALT', 'POS_ALT', 'Population'
        ])
        for pop_ in tqdm(self.populations, desc='Подсчет AF...', unit='pop'):
            af_file = self._get_or_calculate_af_file(pop_)
            df_frequencies = self._update_frequencies(df_frequencies, af_file, pop_)

        df_frequencies = self._infer_types(df_frequencies)
        if save_path:
            df_frequencies.to_csv(save_path, index=False)
        return df_frequencies
    
    def _infer_types(self, df: pd.DataFrame) -> pd.DataFrame:
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='ignore', downcast='float')
        return df

    def _get_or_calculate_af_file(self, pop_):
        af_file = f"{self.data_folder}/af_{pop_}.txt"
        if not (os.path.exists(af_file) and os.path.getsize(af_file) > 0):
            self._calculate_af_for_population(pop_)
        return af_file

    def _calculate_af_for_population(self, pop_):
        samples_file = self._create_samples_file(pop_)
        af_file = f"{self.data_folder}/af_{pop_}.txt"
        bcftools_command = f"bcftools view -S {samples_file} {self.frequency_data_path} | bcftools stats - | grep '^AF' > {af_file}"
        run_cmd(bcftools_command, shell=True)
        self.af_files[pop_] = af_file

    def _create_samples_file(self, pop_):
        samples_file = f"{self.data_folder}/samples_{pop_}.txt"
        awk_command = f"awk '$2 == \"{pop_}\"' {self.population_data_path} | cut -f1 > {samples_file}"
        run_cmd(awk_command, shell=True)
        return samples_file

    def _update_frequencies(self, df_frequencies, af_file, pop_):
        with open(af_file, 'r') as f:
            data = [line.split() for line in f.readlines() if line.startswith('AF')]
        df = pd.DataFrame(data, columns=df_frequencies.columns[:-1])
        df['Population'] = pop_
        return pd.concat([df_frequencies, df], ignore_index=True)

    def _get_common_positions(self, dfs: dict, populations: list) -> list:
        """Находит общие позиции для всех популяций."""
        positions = set(dfs[populations[0]].POS_REF.values)
        for pop_ in populations[1:]:
            positions &= set(dfs[pop_].POS_REF.values)
        return sorted(positions)

    def get_correlation(self, df_frequencies: pd.DataFrame, populations: list = None):
        df_frequencies, populations = self._check_data(df_frequencies, populations)
        dfs = self._split_and_filter_dfs(df_frequencies, populations)
        return self._calculate_correlation(dfs, populations)

    def _split_and_filter_dfs(self, df_frequencies, populations):
        dfs = {pop: df_frequencies[df_frequencies['Population'] == pop] for pop in populations}
        positions = self._get_common_positions(dfs, populations)
        return {pop: df[df.POS_REF.isin(positions)].drop_duplicates('POS_REF').sort_values("POS_REF") for pop, df in dfs.items()}

    def _calculate_correlation(self, df_frequencies: pd.DataFrame, populations: list = None):
        df_frequencies, populations = self._check_data(df_frequencies, populations)
        dfs = self._split_and_filter_dfs(df_frequencies, populations)
        
        self.correlations = {}
        for pop1, pop2 in combinations(populations, 2):
            corr, _ = pearsonr(dfs[pop1]['Allele_Frequency'], dfs[pop2]['Allele_Frequency'])
            self.correlations[f'{pop1}_{pop2}'] = corr.item()
        return self.correlations

    def plot_scatter(self, df_frequencies: pd.DataFrame, populations: list = None):
        df_frequencies, populations = self._check_data(df_frequencies, populations)
        dfs = self._split_and_filter_dfs(df_frequencies, populations)
        self._plot_scatter_for_populations(dfs, populations)

    def _plot_scatter_for_populations(self, dfs, populations):
        pairs = list(combinations(populations, 2))
        fig, axs = plt.subplots(1, len(pairs), figsize=(7 * len(pairs), 8))

        if len(pairs) == 1:
            axs = [axs]

        for i, (pop1, pop2) in enumerate(pairs):
            self._plot_single_scatter(axs[i], dfs[pop1], dfs[pop2], pop1, pop2)

        plt.figlegend(["Генетический вариант", "Линейная регрессия"], loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=20, frameon=False)
        plt.subplots_adjust(wspace=.2)
        plt.show()

    def _plot_single_scatter(self, ax, df1, df2, pop1, pop2):
        af_pop1 = df1['Allele_Frequency'].values.astype(float).reshape(-1, 1)
        af_pop2 = df2['Allele_Frequency'].values.astype(float)

        predictions, r_squared, slope, intercept = self.fit_linear_regression(af_pop1, af_pop2)

        self.slopes[f'{pop1}_{pop2}'] = slope.item()
        self.intercepts[f'{pop1}_{pop2}'] = intercept.item()
        self.r_squares[f'{pop1}_{pop2}'] = r_squared

        sns.scatterplot(x=af_pop1.flatten(), y=af_pop2, ax=ax, color="skyblue", edgecolor='gray', s=50, alpha=0.9)
        ax.plot(af_pop1, predictions, color="salmon", lw=1.5)
        ax.text(0.05, 0.95, f"$R^2={r_squared:.2f}$", transform=ax.transAxes, fontsize=18, verticalalignment='top')
        ax.set_xlabel(f"AF {pop1}", fontsize=18)
        ax.set_ylabel(f"AF {pop2}", fontsize=18)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.set_xlim(0.01, 1.1)
        ax.set_ylim(0.01, 1.1)

    def fit_linear_regression(self, af_pop1, af_pop2):
        model = LinearRegression()
        model.fit(af_pop1, af_pop2)
        predictions = model.predict(af_pop1)
        r_squared = r2_score(af_pop2, predictions)

        return predictions, r_squared, model.coef_[0], model.intercept_

    def _check_data(self, df_frequencies, populations):
        if df_frequencies is None:
            df_frequencies = self.calculate_af()
        if populations is None:
            populations = self.populations
        return df_frequencies, populations
