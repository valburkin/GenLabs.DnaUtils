namespace GenLabs.DnaUtils.Tests;

// Where a test cites a sample dataset, the expected values are the published answers from
// Compeau & Pevzner, "Bioinformatics Algorithms: An Active Learning Approach".
public sealed class SequenceTests
{
    [Fact]
    public void Parses_and_prints_bases_case_insensitively()
    {
        var sequence = new Sequence("acgT");

        Assert.Equal("ACGT", sequence.ToString());
        Assert.Equal([NBase.A, NBase.C, NBase.G, NBase.T], sequence);
    }

    [Fact]
    public void Rejects_characters_that_are_not_bases()
    {
        Assert.Throws<ArgumentException>(() => new Sequence("ACGN"));
    }

    [Fact]
    public void Equal_sequences_are_equal_and_hash_alike()
    {
        Sequence left = "GATTACA";
        Sequence right = "GATTACA";

        Assert.Equal(left, right);
        Assert.True(left == right);
        Assert.True(left.Equals((object)right));
        Assert.Equal(left.GetHashCode(), right.GetHashCode());
        Assert.NotEqual(left, (Sequence)"GATTACT");
    }

    [Fact]
    public void Copies_the_list_it_is_built_from()
    {
        var bases = new[] { NBase.A, NBase.C };
        var sequence = new Sequence(bases);

        bases[0] = NBase.T;

        Assert.Equal("AC", sequence.ToString());
    }

    [Fact]
    public void Reverse_complement_matches_the_sample_dataset()
    {
        Assert.Equal("ACCGGGTTTT", new Sequence("AAAACCCGGT").ReverseComplement().ToString());
    }

    [Fact]
    public void Complement_and_reverse_compose_to_the_reverse_complement()
    {
        Sequence sequence = "GATTACA";

        Assert.Equal("CTAATGT", sequence.Complement().ToString());
        Assert.Equal("ACATTAG", sequence.Reverse().ToString());
        Assert.Equal(sequence.ReverseComplement(), sequence.Complement().Reverse());
    }

    [Fact]
    public void Prepending_and_appending_return_new_sequences()
    {
        Sequence sequence = "CG";

        Assert.Equal("ACGT", (NBase.A + sequence + NBase.T).ToString());
        Assert.Equal("CG", sequence.ToString());
    }

    [Fact]
    public void Mutate_and_part_leave_the_original_unchanged()
    {
        Sequence sequence = "AAAA";

        Assert.Equal("AAGA", sequence.Mutate(2, NBase.G).ToString());
        Assert.Equal("AA", sequence.Part(1, 2).ToString());
        Assert.Equal("AAAA", sequence.ToString());
    }

    [Fact]
    public void Counts_each_base()
    {
        var counts = new Sequence("AACGTTT").Counts;

        Assert.Equal(2, counts[NBase.A]);
        Assert.Equal(1, counts[NBase.C]);
        Assert.Equal(1, counts[NBase.G]);
        Assert.Equal(3, counts[NBase.T]);
        Assert.Equal(7, counts.Total);
        Assert.Equal([(NBase.T, 3)], counts.Max);
    }

    [Fact]
    public void Finds_overlapping_k_mer_occurrences()
    {
        var kMers = new Sequence("GCGCG").FindAllKMers(3);

        Assert.Equal([0, 2], kMers["GCG"]);
        Assert.Equal([1], kMers["CGC"]);
        Assert.Equal(["GCG", "CGC"], kMers.Keys.Select(k => k.ToString()));
    }

    [Fact]
    public void Clump_finding_matches_the_sample_dataset()
    {
        Sequence genome = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA";

        var clumps = genome.FindClumpingKMers(windowSize: 50, k: 5, threshold: 4);

        Assert.Equal(["CGACA", "GAAGA"], clumps.Select(c => c.ToString()).Order());
    }

    [Fact]
    public void A_clump_must_fit_entirely_inside_the_window()
    {
        // "AAA" at 0 and 4: the two occurrences span 7 bases.
        Sequence genome = "AAACAAA";

        Assert.Empty(genome.FindClumpingKMers(windowSize: 6, k: 3, threshold: 2));
        Assert.Equal(["AAA"], genome.FindClumpingKMers(windowSize: 7, k: 3, threshold: 2).Select(c => c.ToString()));
    }

    [Fact]
    public void Minimum_skew_matches_the_sample_dataset()
    {
        Sequence genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGCCGAAC";

        Assert.Equal([11, 24], genome.GetMinSkewPositions());
    }

    [Fact]
    public void Skew_starts_at_zero_and_tracks_g_minus_c()
    {
        Assert.Equal(
            [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2],
            new Sequence("CATGGGCATCGGCCATACGCC").GetSkew());
    }

    [Fact]
    public void Hamming_distance_matches_the_sample_dataset()
    {
        Assert.Equal(3, new Sequence("GGGCCGTTGGT").HammingDistance("GGACCGTTGAC"));
    }

    [Fact]
    public void Hamming_distance_requires_equal_lengths()
    {
        Assert.Throws<ArgumentException>(() => new Sequence("ACG").HammingDistance("AC"));
    }

    [Fact]
    public void Finds_exact_and_approximate_pattern_locations()
    {
        Sequence text = "ACGTACTA";

        Assert.Equal([0], text.FindPatternLocations("ACG"));
        Assert.Equal([0, 4], text.FindPatternLocations("ACG", mismatches: 1));
    }

    [Fact]
    public void Neighbors_match_the_sample_dataset()
    {
        var neighbors = new Sequence("ACG").Neighbors(1).Select(n => n.ToString()).Order();

        Assert.Equal(["AAG", "ACA", "ACC", "ACG", "ACT", "AGG", "ATG", "CCG", "GCG", "TCG"], neighbors);
    }

    [Fact]
    public void Neighbors_within_the_full_length_are_every_sequence()
    {
        Assert.Equal(16, new Sequence("AC").Neighbors(2).Distinct().Count());
    }

    [Fact]
    public void Frequent_words_with_mismatches_match_the_sample_dataset()
    {
        Sequence text = "ACGTTGCATGTCGCATGATGCATGAGAGCT";

        Assert.Equal(["ATGC", "ATGT", "GATG"], MostFrequent(text.FindAllKMers(4, mismatches: 1, treatReverseComplementsAsEqual: false)));
    }

    [Fact]
    public void Frequent_words_with_mismatches_and_reverse_complements_match_the_sample_dataset()
    {
        Sequence text = "ACGTTGCATGTCGCATGATGCATGAGAGCT";

        Assert.Equal(["ACAT", "ATGT"], MostFrequent(text.FindAllKMers(4, mismatches: 1, treatReverseComplementsAsEqual: true)));
    }

    [Fact]
    public void Approximate_k_mer_positions_are_ascending()
    {
        var kMers = new Sequence("ACGTTGCATGTCGCATGATGCATGAGAGCT").FindAllKMers(4, mismatches: 1, treatReverseComplementsAsEqual: true);

        Assert.All(kMers.Values, positions => Assert.Equal(positions.Order(), positions));
    }

    [Fact]
    public void A_position_counts_once_per_strand_the_k_mer_occurs_on()
    {
        // "ACGT" is its own reverse complement, so it occurs at 0 on both strands.
        var kMers = new Sequence("ACGT").FindAllKMers(4, mismatches: 0, treatReverseComplementsAsEqual: true);

        Assert.Equal([0, 0], kMers["ACGT"]);
    }

    [Fact]
    public void Enumerates_every_sequence_of_a_length_in_lexicographic_order()
    {
        var all = Sequence.AllOfLength(2).Select(s => s.ToString()).ToArray();

        Assert.Equal(16, all.Length);
        Assert.Equal(all.Order(StringComparer.Ordinal), all);
        Assert.Equal("AA", all[0]);
        Assert.Equal("TT", all[^1]);
    }

    [Fact]
    public void Random_sequences_are_reproducible_with_a_seed()
    {
        Assert.Equal(NBaseHelper.Random(20, new Random(42)), NBaseHelper.Random(20, new Random(42)));
    }

    private static IEnumerable<string> MostFrequent(Dictionary<Sequence, List<int>> kMers)
    {
        var max = kMers.Values.Max(p => p.Count);
        return kMers.Where(k => k.Value.Count == max).Select(k => k.Key.ToString()).Order();
    }
}
