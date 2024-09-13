import os
import subprocess
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass

from utils.run_cmd import run_cmd

@dataclass
class KlinefelterSearcher:
    """Class for searching for Klinefelter Syndrom from .bam files"""
    data_folder: str

    def __post_init__(self):
        self.chrom = ("X", "Y")
        self.chrom_norm = [f'N_{chr}' for chr in self.chrom] + ['N_total']

    def get_chromosome_reads(self, bam_file_path: str):
        cmd = ["samtools", "idxstats", bam_file_path]
        result = run_cmd(cmd).stdout.strip().split("\n")

        sample_info = {}.fromkeys(self.chrom_norm, 0)

        C_total, L_total = 0, 0

        for line in result:
            # 1. Извлечение количества прочтений (C), пришедшихся на хромосому, и ее длины (L)
            cols = line.split("\t")
            chr_name, L, C = cols[0], int(cols[1]), int(cols[2])

            if L == 0:
                continue

            # 2. Нормализация числа прочтений на длину хромосомы
            N = C / L

            if chr_name in self.chrom:
                sample_info[f'N_{chr_name}'] = N

            C_total += C
            L_total += L

        sample_info['N_total'] = C_total / L_total
        return sample_info
    
    def get_normalized_reads(self, save_data=False):
        data_karyo = pd.DataFrame(columns = self.chrom_norm)
        
        for bam_file_path in tqdm(os.listdir(self.data_folder), desc='Process .bam files', leave=True, colour='green'):
            if not bam_file_path.endswith(".bam"):
                continue
            sample_name = bam_file_path.split(".")[0]
            bam_file_path = os.path.join(self.data_folder, bam_file_path)
            sample_info = self.get_chromosome_reads(bam_file_path)
            sample_df = pd.DataFrame(sample_info, index=[sample_name])
            data_karyo = pd.concat([data_karyo, sample_df])

        # Нормализуем на (нормализованное) общее число прочтений
        data_karyo['N_X_norm'] = data_karyo['N_X'] / data_karyo['N_total']
        data_karyo['N_Y_norm'] = data_karyo['N_Y'] / data_karyo['N_total']
        
        # Вычисляем Ratio = N_X / N_Y
        data_karyo['Ratio'] = data_karyo['N_X'] / data_karyo['N_Y']

        if save_data:
            save_data_path = os.path.join(self.data_folder, 'karyotypes.csv')
            data_karyo.to_csv(save_data_path, sep=',', index=False)
        return data_karyo

    def plot_karyotypes(self, data_karyo, save_fig=False):
        # сортируем по значениям Ratio (возрастающий порядок)
        data_karyo = data_karyo.sort_values(by='Ratio', ascending=False)

        fig, ax = plt.subplots(figsize=(10, 6))

        # Рисуем N_X_norm и N_Y_norm
        bar1 = sns.barplot(x=data_karyo.index, y=data_karyo['N_X_norm'], 
                           color='salmon', edgecolor='grey',
                           label='X chr', ax=ax)

        bar2 = sns.barplot(x=data_karyo.index, y=data_karyo['N_Y_norm'],
                           estimator=sum, ci=None, 
                           color='skyblue', edgecolor='grey',
                           label = "Y chr", ax=ax)

        # Отображаем Ratio
        for i in range(len(data_karyo.index)):
            plt.text(i, 
                    max(data_karyo['N_X_norm'][i], data_karyo['N_Y_norm'][i]) + 0.05, 
                    round(data_karyo['Ratio'][i], 2),
                    ha = 'center',
                    color='k', 
            )

        # Не забываем подписывать оси
        plt.xlabel("Samples", fontsize=12)
        plt.ylabel(r"$\frac{\text{Read count}}{\text{length}}$", fontsize=16)

        plt.xticks(rotation=45)
        plt.ylim(top=2.5)

        # Добавляем Ratio в легенду
        handles, _ = ax.get_legend_handles_labels()
        ax.legend(handles + [plt.Line2D([0], [0], color='black', lw=1)],
                [r'$N_X$', r'$N_Y$', 'Ratio'], 
                loc='upper right', facecolor='lightgrey', 
                fontsize='large', edgecolor='black')

        plt.show()

        if save_fig:
            save_data_path = os.path.join(self.data_folder, 'karyotypes.png')
            fig.savefig(save_data_path)