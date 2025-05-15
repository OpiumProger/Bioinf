from Bio import Entrez, SeqIO, Seq
from Seq import get_codons_for_amino_acid, allSeq
from ParsOrganise import fetch_codon_table, optimize_codon
from GC_show import plot_gc_profile
from ML import GCOptimizer


Id = 0
CodonList = []
Codons = {}
Entrez.email = input("Введите название почты аккаунта на NCBI: ") # Прописывается почта аккаунта на NCBI, с которого будут производиться запросы
DataBase = input("Введите название бд: ")
SearchID = input("Введите название белка: ")   #"AAA40590.1" пример
try:
    with Entrez.efetch(
            db=DataBase,
            id=SearchID,
            rettype="fasta",
            retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "fasta")
        ID = record.id
        Name = record.name
        Description = record.description
        Seq = record.seq
        amino_acid_seq = Seq
        with open("output.fasta", "w") as f:
            f.write(f">{record.description}\n{record.seq}")
        print(f"Описание: {record.description}")
        print(f"Последовательность: {Seq}")

        for aminos in Seq:
            amino = aminos
            codons = get_codons_for_amino_acid(amino)
            Id += 1
            Codons[Id] = codons
            CodonList.append(codons)

except Exception as e:
    print(f"Ошибка: {e}")


species_id = input("Введите id организма: ")          #83333 пример для E.coli
fetch_codon_table(species_id)
codon_table = fetch_codon_table(species_id)
opt_RNA = optimize_codon(amino_acid_seq, codon_table)

# Прописываем нашу модель и параметры обучения
optimizer = GCOptimizer(
    codon_table=codon_table,
    target_gc=0.51,  # Целевой GC состав
    n_iter=120,     # Число итераций
    pop_size=50    # Размер принимаемых вариантов
)

# Запуск алгоритма оптимизации
optimizer.fit(initial_dna=opt_RNA, aa_seq=amino_acid_seq)

#Ниже прописан пример, который я использовал для сравнения со своим алгоритмом на основе онлайнк калукулятора как опт. послед. с CAI = 0.82 +-+

check_Seq = ''' ATGGCACCTTGGATGCACCTGCTGACTGTTCTGGCGCTGCTGGCCCTGTGGGGTCCGAACTCTGTGCAGGCGTACAGCT
 CTCAGCACCTGTGCGGTAGCAACCTGGTTGAGGCTCTGTACATGACGTGTGGCCGCTCCGGTTTTTACCGTCCACACGA
 CCGTCGTGAACTGGAAGACCTGCAGGTTGAACAGGCCGAACTGGGTCTGGAAGCTGGTGGTCTGCAGCCGAGCGCGCTG
 GAAATGATTCTGCAGAAGCGTGGCATCGTCGACCAGTGTTGCAACAACATCTGCACCTTCAACCAGCTGCAGAACTACT
 GCAACGTTCCA'''

# Результат
optimized_dna = optimizer.best_dna_

plot_gc_profile(check_Seq, optimized_dna) # Изображение GC состава