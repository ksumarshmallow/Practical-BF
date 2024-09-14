import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import combinations
from sklearn.linear_model import LinearRegression
from dataclasses import dataclass

@dataclass
class OutliersSearcher:
    df_frequencies: pd.DataFrame
    
    def __post_init__(self):
        self.populations = self.df_frequencies['Population'].unique()
    
    def find_outliers_iqr(self, data):
        """Находит выбросы с помощью метода межквартильного размаха."""
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        return np.where((data < lower_bound) | (data > upper_bound))

    def find_outliers_zscore(self, data, threshold=3):
        """Находит выбросы на основе z-оценок."""
        mean = np.mean(data)
        std = np.std(data)
        z_scores = (data - mean) / std
        return np.where(np.abs(z_scores) > threshold)

    def find_outliers_residuals(self, af_pop1, af_pop2, method, threshold):
        """Находит выбросы на основе остатков линейной регрессии."""
        model = LinearRegression()
        model.fit(af_pop1.reshape(-1, 1), af_pop2)
        predictions = model.predict(af_pop1.reshape(-1, 1))
        
        # Вычисляем остатки (residuals)
        residuals = af_pop2 - predictions
        
        # Находим выбросы 
        if method == 'iqr':
            outliers_idx = self.find_outliers_iqr(residuals)
        elif method == 'zscore':
            outliers_idx = self.find_outliers_zscore(residuals, threshold=threshold)
        
        return outliers_idx, residuals

    def search(self, method='iqr', threshold=3):
        dfs = self._split_and_filter_dfs()

        outliers_pos_dict = {}
        for pop1, pop2 in tqdm(combinations(self.populations, 2), desc='Search outliers...', colour='GREEN'):
            af_pop1 = dfs[pop1]['Allele_Frequency'].values.astype(float).reshape(-1, 1)
            af_pop2 = dfs[pop2]['Allele_Frequency'].values.astype(float)
            outliers_idx, _ = self.find_outliers_residuals(af_pop1, af_pop2, method, threshold)
            
            outliers_pos = pd.concat((dfs[pop1].iloc[outliers_idx[0]], dfs[pop2].iloc[outliers_idx[0]]))
            outliers_pos_dict[f'{pop1}_{pop2}'] = outliers_pos
        return outliers_pos_dict

    def _split_and_filter_dfs(self):
        dfs = {pop: self.df_frequencies[self.df_frequencies['Population'] == pop] for pop in self.populations}
        positions = self._get_common_positions(dfs)
        return {pop: df[df.POS_REF.isin(positions)].drop_duplicates('POS_REF').sort_values("POS_REF") for pop, df in dfs.items()}
    
    def _get_common_positions(self, dfs: dict) -> list:
        """Находит общие позиции для всех популяций."""
        positions = set(dfs[self.populations[0]].POS_REF.values)
        for pop_ in self.populations[1:]:
            positions &= set(dfs[pop_].POS_REF.values)
        return sorted(positions)
