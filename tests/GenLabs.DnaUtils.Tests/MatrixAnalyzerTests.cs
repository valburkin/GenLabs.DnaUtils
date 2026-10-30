namespace GenLabs.DnaUtils.Tests;

public sealed class MatrixAnalyzerTests
{
    [Fact]
    public void Greedy_motif_search_with_pseudocounts_matches_the_sample_dataset()
    {
        Sequence[] dna =
        [
            "GGCGTTCAGGCA",
            "AAGAATCAGTCA",
            "CAAGGAGTTCGC",
            "CACGTCAATCAC",
            "CAATAATATTCG",
        ];

        var motifs = new MatrixAnalyzer(dna).SearchMotifs(3);

        Assert.Equal(["TTC", "ATC", "TTC", "ATC", "TTC"], motifs.Select(m => m.ToString()));
    }

    [Fact]
    public void Score_consensus_and_entropy_describe_the_columns()
    {
        // Column 1 is all A; column 2 is half A, half C.
        var analyzer = new MatrixAnalyzer([(Sequence)"AA", "AC", "AA", "AC"]);

        Assert.Equal(2, analyzer.Score);
        Assert.Equal([NBase.A, NBase.A], analyzer.Consensus);
        Assert.Equal(1.0, analyzer.Entropy, precision: 10);
    }

    [Fact]
    public void A_fully_conserved_matrix_scores_zero_with_zero_entropy()
    {
        var analyzer = new MatrixAnalyzer([(Sequence)"GATTACA", "GATTACA"]);

        Assert.Equal(0, analyzer.Score);
        Assert.Equal(0.0, analyzer.Entropy);
        Assert.Equal("GATTACA", new Sequence(analyzer.Consensus).ToString());
    }

    [Fact]
    public void Probability_uses_the_pseudocounted_profile()
    {
        // One region "A": A is (1 + 1) / (1 + 4), every other base (0 + 1) / (1 + 4).
        var analyzer = new MatrixAnalyzer([(Sequence)"AC"]);

        Assert.Equal(0.4 * 0.4, analyzer.Probability("AC"), precision: 10);
        Assert.Equal(0.2 * 0.2, analyzer.Probability("GT"), precision: 10);
    }

    [Fact]
    public void Most_probable_pattern_prefers_the_first_on_a_tie()
    {
        var analyzer = new MatrixAnalyzer([(Sequence)"AC"]);

        Assert.Equal("AC", analyzer.GetMostProbablePattern([(Sequence)"GT", "AC", "AG"]).ToString());
        Assert.Equal("GT", analyzer.GetMostProbablePattern([(Sequence)"GT", "TG"]).ToString());
    }

    [Fact]
    public void Probability_requires_a_pattern_as_long_as_the_regions()
    {
        var analyzer = new MatrixAnalyzer([(Sequence)"ACG"]);

        Assert.Throws<ArgumentException>(() => analyzer.Probability("AC"));
    }

    [Fact]
    public void Rejects_regions_of_different_lengths()
    {
        Assert.Throws<ArgumentException>(() => new MatrixAnalyzer([(Sequence)"ACG", "AC"]));
    }

    [Fact]
    public void Rejects_an_empty_matrix()
    {
        Assert.Throws<ArgumentException>(() => new MatrixAnalyzer([]));
    }
}
