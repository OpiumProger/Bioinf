from sklearn.base import BaseEstimator
import numpy as np
from Bio.Seq import Seq
import optuna
from ParsOrganise import fetch_codon_table
import subprocess



class GCOptimizer(BaseEstimator):
    def __init__(self, codon_table, species_id=None, target_gc=0.55, n_iter=150, pop_size=45,
                 end_zone=70, min_end_gc=40, max_end_gc=70):
        """
        species_id: ID организма из Kazusa (например, 83333 для E.coli)
        Если None - тогда используется таблци с равномерным распределением кодонов
        """
        self.codon_table = codon_table
        self.species_id = species_id
        self.target_gc = target_gc
        self.n_iter = n_iter
        self.min_end_gc = min_end_gc
        self.max_end_gc = max_end_gc
        self.pop_size = pop_size
        self.end_zone = end_zone
        self.docker_image = "medicinalgenomics/viennarna"
        # Загружаем таблицу частот кодонов
        self.codon_freq = fetch_codon_table(species_id) if species_id else None

    def _evaluate(self, dna, aa_seq):
        """Оценка качества последовательности с учетом GC-состава и CAI"""
        # Проверка соответствия белку
        try:
            assert str(Seq(dna).translate()) == aa_seq
        except:
            return -1

        # Расчет GC-состава
        gc = sum(1 for base in dna if base in ('G', 'C')) / len(dna)
        gc_score = 1 - abs(gc - self.target_gc)

        # Контроль GC на концах
        end5_gc = (dna[:self.end_zone].count('G') + dna[:self.end_zone].count('C')) / self.end_zone
        end3_gc = (dna[-self.end_zone:].count('G') + dna[-self.end_zone:].count('C')) / self.end_zone

        # Штрафы за выход за границы
        end_penalty = (
                              max(0, self.min_end_gc / 100 - end5_gc) +
                              max(0, end5_gc - self.max_end_gc / 100) +
                              max(0, self.min_end_gc / 100 - end3_gc) +
                              max(0, end3_gc - self.max_end_gc / 100)
                      ) * 2

        # Расчет CAI (если имеем данные о частотах)
        cai_score = self._calculate_cai(dna) if self.codon_freq else 1.0





        # Комбинированная оценка
        return 0.7 * gc_score + 0.3 * cai_score - end_penalty


    def _calculate_cai(self, dna):
        """Вычисляет индекс адаптации кодонов (CAI)"""
        cai = 1.0
        count = 0

        for i in range(0, len(dna), 3):
            codon = dna[i:i + 3]
            aa = str(Seq(codon).translate())

            if aa in self.codon_freq:
                # Находим максимальную частоту для этой аминокислоты
                max_freq = max(freq for _, freq in self.codon_freq[aa])

                # Находим частоту текущего кодона
                for c, freq in self.codon_freq[aa]:
                    if c == codon:
                        cai *= (freq / max_freq)
                        count += 1
                        break

        return cai ** (1 / count) if count > 0 else 1.0

    def _mutate(self, dna, aa_seq):
        """Мутация одного кодона с учетом частот (если доступны)"""
        pos = np.random.randint(0, len(aa_seq))
        aa = aa_seq[pos]

        if self.codon_freq and aa in self.codon_freq:
            # Выбираем кодон по частотам использования
            codons, freqs = zip(*self.codon_freq[aa])
            probas = np.array(freqs) / sum(freqs)
            new_codon = codons[np.random.choice(len(codons), p=probas)]
        else:
            # Резервный вариант - равномерный выбор
            codons = self.codon_table[aa]
            new_codon = codons[np.random.randint(0, len(codons))][0]

        return dna[:pos * 3] + new_codon + dna[pos * 3 + 3:]

    def fit(self, initial_dna, aa_seq):
        """Генетический алгоритм с оптимизацией параметров через Optuna"""

        # Оптимизация параметров min_end_gc и max_end_gc
        def objective(trial):
            self.min_end_gc = trial.suggest_float("min_end_gc", 38, 50)
            self.max_end_gc = trial.suggest_float("max_end_gc", 55, 70)

            population = [initial_dna]
            for _ in range(self.pop_size - 1):
                population.append(self._mutate(initial_dna, aa_seq))

            for _ in range(self.n_iter):
                scores = [self._evaluate(dna, aa_seq) for dna in population]
                best = population[np.argmax(scores)]
                new_pop = [best]
                for _ in range(self.pop_size - 1):
                    parent = population[np.random.choice(
                        range(self.pop_size),
                        p=np.array(scores) / sum(scores)
                    )]
                    new_pop.append(self._mutate(parent, aa_seq))
                population = new_pop

            best_dna = max(population, key=lambda x: self._evaluate(x, aa_seq))
            return self._evaluate(best_dna, aa_seq)

        study = optuna.create_study(direction="maximize")
        study.optimize(objective, n_trials=30)

        # Устанавливаем лучшие параметры
        self.min_end_gc = study.best_params["min_end_gc"]
        self.max_end_gc = study.best_params["max_end_gc"]
        print(f"Оптимальные параметры: min_end_gc={self.min_end_gc}, max_end_gc={self.max_end_gc}")

        # Финальный запуск с лучшими параметрами
        population = [initial_dna]
        for _ in range(self.pop_size - 1):
            population.append(self._mutate(initial_dna, aa_seq))

        for _ in range(self.n_iter):
            scores = [self._evaluate(dna, aa_seq) for dna in population]
            best = population[np.argmax(scores)]
            new_pop = [best]
            for _ in range(self.pop_size - 1):
                parent = population[np.random.choice(
                    range(self.pop_size),
                    p=np.array(scores) / sum(scores)
                )]
                new_pop.append(self._mutate(parent, aa_seq))
            population = new_pop

        self.best_dna_ = max(population, key=lambda x: self._evaluate(x, aa_seq))

        # Вывод результатов
        gc_content = (self.best_dna_.count('G') + self.best_dna_.count('C')) / len(self.best_dna_)
        print(f"GC после оптимизации: {gc_content * 100:.1f}%")
        if self.codon_freq:
            print(f"CAI: {self._calculate_cai(self.best_dna_):.3f}")
        print(f"Лучшая последовательность: {self.best_dna_}")

        return self