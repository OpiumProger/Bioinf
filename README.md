# Bioinf
My progects

The GCOptimizer class implements an intelligent genetic sequence optimization system that:
   • Selects optimal codons for a given amino acid sequence
   • Controls the overall GC composition according to the target value
   • Regulates GC composition in the terminal regions of the gene
   • Takes into account codon preferences for specific organisms (if species_id is specified)

   
 # Key Features
  1. Adaptive Genetic Algorithm
    •	Uses a population-based approach with controlled mutations
    •	Automatically optimizes terminal GC zones using Optuna hyperparameter tuning
    •	Combines multiple optimization criteria into a single fitness function
  2. Multi-Objective Sequence Evaluation
    The _evaluate method assesses:
    •	Translation accuracy (strict match to target amino acid sequence)
    •	GC content (proximity to target_gc)
    •	Codon Adaptation Index (CAI) (optimization for host organism preferences)
    •	Terminal GC zones (special requirements for 5' and 3' ends)

  3. Biologically Informed Mutations
    The _mutate method implements:
    •	Smart codon selection based on organism-specific usage frequencies
    •	Fallback to uniform distribution if no frequency data is available
    •	Point mutations that preserve the amino acid sequence

# Algorithm Workflow
  1. Initialization
    •	Loads codon tables and (optionally) codon usage frequencies
    •	Configures optimization parameters (target_gc, n_iter, pop_size, etc.)
  2. Hyperparameter Optimization (Optuna)
    •	Automatically searches for the best GC boundaries in terminal regions
    •	Runs 30 trials to find the global optimum for min_end_gc and max_end_gc
  3. Core Optimization Process
    •	Generates an initial population of DNA variants
    •	Iteratively improves sequences through:
    o	Fitness evaluation
    o	Elite selection
    o	Controlled mutations
    •	Evolves over 150 generations
  4. Results & Output
    •	Returns the optimized DNA sequence
    •	Provides detailed statistics on:
    o	Final GC content
    o	Codon Adaptation Index (CAI)
    o	Terminal GC compliance
   
   ⭐ Flexibility – Works with or without organism-specific preferences
   ⭐ Balanced Optimization – Considers both global GC and terminal stability
   ⭐ Automation – Minimizes manual tuning, maximizes adaptability
   ⭐ Reproducibility – Deterministic results for identical inputs









