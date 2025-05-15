import matplotlib.pyplot as plt


def plot_gc_profile(sequence1, sequence2, window=30):
    """Визуализирует GC-состав двух последовательностей на одном изображении"""

    # Рассчитываем GC-состав для обеих последовательностей
    def calculate_gc(sequence, window_size):
        return [(sequence[i:i + window_size].count('G') + sequence[i:i + window_size].count('C')) / window_size * 100
                for i in range(len(sequence) - window_size)]

    gc_values1 = calculate_gc(sequence1, window)
    gc_values2 = calculate_gc(sequence2, window)

    # Создаем фигуру с двумя подграфиками
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # График первой последовательности
    ax1.plot(gc_values1, color='blue')
    ax1.axhline(y=40, color='r', linestyle='--')
    ax1.axhline(y=60, color='r', linestyle='--')
    ax1.set_ylabel('GC %')
    ax1.set_xlabel('Position')
    ax1.set_title('GC Content - Seq 1')

    # График второй последовательности
    ax2.plot(gc_values2, color='green')
    ax2.axhline(y=40, color='r', linestyle='--')
    ax2.axhline(y=60, color='r', linestyle='--')
    ax2.set_ylabel('GC %')
    ax2.set_xlabel('Position')
    ax2.set_title('GC Content - Seq 2')

    # Оптимизируем расположение и показываем графики
    plt.tight_layout()
    plt.show()