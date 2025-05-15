import requests
from bs4 import BeautifulSoup
from Bio.Data import CodonTable
import matplotlib.pyplot as plt
import numpy as np
def fetch_codon_table(species_id: int):
    url = f"https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={species_id}"

    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, 'html.parser')
        pre_tag = soup.find('pre')
        if not pre_tag:
            print(f"Ошибка: Тег <pre> не найден для species ID {species_id}")
            return None

        arr = pre_tag.text.strip().split()

        # Создаем словарь с частотами кодонов
        codon_freq = {}
        for i in range(0, len(arr), 3):
            if i + 1 >= len(arr):
                break
            codon = arr[i]
            freq_str = arr[i + 1]
            try:
                freq = float(freq_str.split('(')[0])
                codon_freq[codon] = freq
            except (ValueError, IndexError):
                continue

        # Получаем таблицу соответствия кодонов и аминокислот
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        codon_to_aa = {}

        # Обрабатываем обычные кодоны (DNA -> RNA)
        for dna_codon, aa in standard_table.forward_table.items():
            rna_codon = dna_codon.replace('T', 'U')
            codon_to_aa[rna_codon] = aa

        # Обрабатываем стоп-кодоны
        for dna_codon in standard_table.stop_codons:
            rna_codon = dna_codon.replace('T', 'U')
            codon_to_aa[rna_codon] = '*'

        # Группируем по аминокислотам
        aa_dict = {}
        for codon, freq in codon_freq.items():
            aa = codon_to_aa.get(codon)

            if not aa or aa == '*':  # Пропускаем неизвестные и стоп-кодоны
                continue

            if aa not in aa_dict:
                aa_dict[aa] = []
            aa_dict[aa].append((codon, freq))

        # Сортируем кодоны внутри аминокислот
        for aa in aa_dict:
            aa_dict[aa].sort(key=lambda x: x[0])

        return aa_dict

    except requests.exceptions.RequestException as e:
        print(f"Ошибка при получении URL: {e}")
        return None
    except Exception as e:
        print(f"Произошла непредвиденная ошибка: {e}")
        return None

def optimize_codon(amino_acid_sequence, codon_table):
    optimized_dna = []
    for aa in amino_acid_sequence:
        codons = codon_table.get(aa, [])
        if not codons:
            raise ValueError(f"Аминокислота {aa} не найдена в таблице кодонов.")
        # Выбор кодона с максимальной частотой
        best_codon = max(codons, key=lambda x: x[1])[0]
        optimized_dna.append(best_codon)
    return ''.join(optimized_dna)











