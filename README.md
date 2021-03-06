
# Ant Colony Optimization for Motif Searching
## Algorithm
*Coming soon*

## Usage
### Step 1: Compile
```
g++ main.cpp -o main 
```
### Step 2: Execute
```
./main [ARGS] [OPTIONS]

Arguments:
    Key         Value                                           Required/Default
------------------------------------------------------------------------------------
    -i          [input_file]                                    required
    -o          [output_file]                                   standard output
    -w          [motif_length]                                  motif length
    -rho        [pheromone evaporation coefficient]             0.5
    -alpha      [exponential parameter of trail pheromone]      1.0
    -beta       [exponential parameter of trail desirability]   1.0
    -min-pher   [minimum pheromone]                             0.1
    -max-pher   [maximum pheromone]                             1.0
    -ant        [number of ants]                                20
    -gen        [number of generations]                         100
    -offset     [threshold of distribution to pick offsets]     0.05

Options:
    -dna: Using DNA Alphabet
    -recalc-background: Background probabilities of characters will be estimated based on the frequencies on the database. Or else, uniform. 			
```