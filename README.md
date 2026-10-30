# GenLabs.DnaUtils

[![CI](https://github.com/valburkin/GenLabs.DnaUtils/actions/workflows/ci.yml/badge.svg)](https://github.com/valburkin/GenLabs.DnaUtils/actions/workflows/ci.yml)

Genomic sequence analysis for .NET. GenLabs.DnaUtils implements the core algorithms for finding the origin of replication and regulatory motifs in DNA, as presented in Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*. The tests check each algorithm against the book's published sample datasets.

| Area | What it does |
|---|---|
| Sequences | An immutable `Sequence` of bases, with complement, reverse complement, mutation and slicing |
| k-mers | All k-mers with their positions; approximate occurrences with mismatches, on one strand or both |
| Clumps | k-mers that occur at least *t* times within a window of *L* bases |
| GC skew | Cumulative G − C skew and its minimum, which points to the origin of replication |
| Pattern matching | Exact and approximate pattern locations, Hamming distance, d-neighborhoods |
| Motifs | Profile matrices with pseudocounts, score, consensus, entropy, greedy motif search |

## Quick start

Requires .NET 10.

```csharp
using GenLabs.DnaUtils;

Sequence genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGCCGAAC";

genome.GetMinSkewPositions();      // [11, 24]
genome.ReverseComplement();        // GTTCGGCCCGCGTTTACGCCCCTCGTTCTAGCACTCGTGTTGGCCTCTCGGCAGTCTTTA

// The most frequent 4-mers with up to 1 mismatch, counting both strands.
var kMers = new Sequence("ACGTTGCATGTCGCATGATGCATGAGAGCT")
    .FindAllKMers(k: 4, mismatches: 1, treatReverseComplementsAsEqual: true);
var max = kMers.Values.Max(p => p.Count);
kMers.Where(k => k.Value.Count == max).Select(k => k.Key);   // ATGT, ACAT

// Greedy motif search across the upstream regions of co-regulated genes.
var motifs = new MatrixAnalyzer(
[
    "GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG",
]).SearchMotifs(k: 3);                                       // TTC, ATC, TTC, ATC, TTC
```

## Sample: the origin of replication in *E. coli*

[`samples/OriFinder`](samples/OriFinder/Program.cs) runs the library on a whole bacterial genome from [GenLabs.Genome](https://github.com/valburkin/GenLabs.Genome). It finds where the GC skew is lowest and checks how often a known DnaA box (the 9-mer that the replication-initiator protein DnaA binds) occurs nearby:

```console
$ dotnet run -c Release --project samples/OriFinder -- escherichia_coli.txt TTATCCACA
Genome: 4,639,675 bases
Minimum skew at: 3923620, 3923621, 3923622, 3923623
In the 500 bases from 3923620, 34 9-mers tie as most frequent (4 occurrences, up to 1 mismatch, either strand).
TTATCCACA occurs 4 times (reverse complement TGTGGATAA): rank 1 by frequency.
Approximate matches on the forward strand: 3923929, 3924003
Approximate matches on the reverse strand: 3923823, 3923964
```

The skew minimum matches the book's result. The experimentally confirmed DnaA box ties for the most frequent 9-mer near it. It doesn't stand out on frequency alone, which is why the book follows skew and frequency with other evidence.

## Build and test

```bash
dotnet build GenLabs.DnaUtils.slnx
dotnet test GenLabs.DnaUtils.slnx
```

Analyzer and code-style warnings fail the build.

## License

[MIT](LICENSE)
