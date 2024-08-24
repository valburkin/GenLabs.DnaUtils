namespace GenLabs.DnaUtils;

/// <summary>
/// Use this class to perform various analyses on a matrix of DNA regions.
/// </summary>
public sealed class MatrixAnalyzer
{
    private readonly IReadOnlyList<Sequence> _matrix;
    private readonly IReadOnlyList<ProfiledNucleotideCount> _profiled;

    /// <summary>
    /// Initializes a new instance of the <see cref="MatrixAnalyzer"/> class.
    /// </summary>
    /// <param name="matrix">The regions of DNA to perform analyses.</param>
    public MatrixAnalyzer(IReadOnlyList<Sequence> matrix)
    {
        _matrix = matrix;

        var columns = CountColumns(_matrix).ToArray();
        _profiled = columns
            .Select(c => c.Normalize(_matrix.Count))
            .ToArray();

        Score = GetScore(columns);
        Entropy = _profiled.Sum(m => m.Entropy);
        Consensus = _profiled
            .Select(m => m.Max.First().Nucleotide)
            .ToArray();
    }

    /// <summary>
    /// Gets the score of the regions matrix.
    /// The score shows how "conserved" the matrix of regions is: the lower the score, the more conserved the matrix is.
    /// Calculated as the sum of the differences between the total count of nucleotides in each column and the count of a nucleotide that appears the most times in that column.
    /// For example, if the regions matrix has 10 regions and the first column is all A's, the score of the column is 0.
    /// If the second column has 7 A's, 2 C's, and 1 G, the score of the column is 3.
    /// </summary>
    public int Score { get; }

    /// <summary>
    /// Gets the consensus nucleotides of the regions matrix.
    /// The consensus nucleotides are the nucleotides that appear the most times in each column of the regions matrix.
    /// </summary>
    public Nucleotide[] Consensus { get; }

    /// <summary>
    /// Gets the entropy of the regions matrix.
    /// Entropy is a measure of the randomness of nucleotides in the matrix.
    /// It is the sum of the entropy of each column in the matrix.
    /// </summary>
    public double Entropy { get; }

    /// <summary>
    /// Searches for regulatory motifs in the regions matrix.
    /// </summary>
    /// <param name="k">The length of the motifs to search for.</param>
    /// <returns>The regulatory motifs found in the regions matrix.</returns>
    public IReadOnlyList<Sequence> SearchMotifs(int k)
    {
        var bestMotifs = _matrix
            .Select(region => region.Part(0, k))
            .ToList();

        var bestScore = GetScore(bestMotifs);

        for (var i = 0; i < _matrix[0].Length - k + 1; i++)
        {
            var motifs = new List<Sequence>(_matrix.Count)
            {
                _matrix[0].Part(i, k)
            };

            for (var j = 1; j < _matrix.Count; j++)
            {
                var patterns = _matrix[j].FindAllKMers(k).Keys.ToArray();
                var profiled = Profile(motifs);

                var mostProbable = GetMostProbablePattern(patterns, profiled);

                motifs.Add(new Sequence(mostProbable));
            }

            if (GetScore(motifs) < bestScore)
            {
                bestMotifs = motifs;
                bestScore = GetScore(bestMotifs);
            }
        }

        return bestMotifs;
    }

    /// <summary>
    /// Gets the probability of a pattern occurring in the regions matrix.
    /// </summary>
    /// <param name="pattern">The pattern to get the probability for.</param>
    /// <returns>The probability of the pattern occurring in the regions matrix.</returns>
    public double Probability(IReadOnlyList<Nucleotide> pattern) =>
        Probability(pattern, _profiled);

    /// <summary>
    /// Gets the probability of a pattern occurring in the regions matrix.
    /// </summary>
    /// <param name="pattern">The pattern to get the probability for.</param>
    /// <returns>The probability of the pattern occurring in the regions matrix.</returns>
    public double Probability(string pattern) =>
        Probability(pattern.Select(c => c.ToNucleotide()).ToArray());

    /// <summary>
    /// Gets the most probable pattern from the given patterns.
    /// </summary>
    /// <param name="patterns">The patterns to find the most probable from.</param>
    /// <returns>The most probable pattern.</returns>
    public IReadOnlyList<Nucleotide> GetMostProbablePattern(
        IReadOnlyList<IReadOnlyList<Nucleotide>> patterns) =>
            GetMostProbablePattern(patterns, _profiled);

    /// <summary>
    /// Gets the probability of a pattern occurring in the regions matrix.
    /// </summary>
    /// <param name="pattern">The pattern to get the probability for.</param>
    /// <param name="normalized">The normalized nucleotide counts of the regions matrix.</param>
    /// <returns>The probability of the pattern occurring in the regions matrix.</returns>
    private static double Probability(
        IReadOnlyList<Nucleotide> pattern,
        IReadOnlyList<ProfiledNucleotideCount> normalized) =>
            pattern
                .Select((n, i) => normalized[i][n])
                .Aggregate(1.0, (current, p) => current * p);

    /// <summary>
    /// Gets the most probable pattern from the given patterns.
    /// </summary>
    /// <param name="patterns">The patterns to find the most probable from.</param>
    /// <param name="normalized">The normalized nucleotide counts of the regions matrix.</param>
    /// <returns>The most probable pattern.</returns>
    private static IReadOnlyList<Nucleotide> GetMostProbablePattern(
        IReadOnlyList<IReadOnlyList<Nucleotide>> patterns,
        IReadOnlyList<ProfiledNucleotideCount> normalized)
    {
        var maxProbability = 0.0;
        var mostProbable = patterns[0];

        foreach (var pattern in patterns)
        {
            var probability = Probability(pattern, normalized);
            if (probability > maxProbability)
            {
                maxProbability = probability;
                mostProbable = pattern;
            }
        }

        return mostProbable;
    }

    /// <summary>
    /// Counts the number of occurrences of each <see cref="Nucleotide"/> in each column of the given regions.
    /// </summary>
    /// <param name="regions">The regions to count columns in.</param>
    /// <returns><see cref="NucleotideCount"/>: The count of each nucleotide in each column.</returns>
    /// <exception cref="ArgumentException">Thrown when the regions are not all the same length.</exception>
    private static IEnumerable<NucleotideCount> CountColumns(IReadOnlyList<Sequence> regions)
    {
        var regionLength = regions[0].Length;

        var columns = new Nucleotide[regionLength][];
        for (var i = 0; i < regionLength; i++)
        {
            columns[i] = new Nucleotide[regions.Count];
        }

        for (var i = 0; i < regions.Count; i++)
        {
            if (regions[i].Length != regionLength)
            {
                throw new ArgumentException("All motifs must be the same length.");
            }

            for (var j = 0; j < regionLength; j++)
            {
                columns[j][i] = regions[i][j];
            }
        }

        return columns.Select(c => new NucleotideCount(c));
    }

    private static ProfiledNucleotideCount[] Profile(IReadOnlyList<Sequence> matrix) =>
        CountColumns(matrix)
            .Select(c => c.Normalize(matrix.Count))
            .ToArray();

    /// <summary>
    /// Gets the score of the regions matrix.
    /// Calculated as the sum of the differences between the total count of nucleotides in each column and the count of a nucleotide that appears the most times in that column.
    /// </summary>
    /// <param name="columns">The columns to get the score for.</param>
    /// <returns>The score of the regions matrix.</returns>
    private static int GetScore(IEnumerable<NucleotideCount> columns) =>
        columns.Sum(column => column.Total - column.Max.First().Count);

    /// <summary>
    /// Gets the score of the regions matrix from the given regions.
    /// </summary>
    /// <param name="regions">The regions to get the score for.</param>
    /// <returns>The score of the regions matrix.</returns>
    private static int GetScore(IReadOnlyList<Sequence> regions) =>
        GetScore(CountColumns(regions));
}