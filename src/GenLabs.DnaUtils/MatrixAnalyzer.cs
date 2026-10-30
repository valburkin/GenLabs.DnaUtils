namespace GenLabs.DnaUtils;

/// <summary>
/// Use this class to perform various analyses on a matrix of DNA regions of equal length,
/// such as the upstream regions of co-regulated genes.
/// </summary>
public sealed class MatrixAnalyzer
{
    private readonly IReadOnlyList<Sequence> _matrix;
    private readonly ProfiledNBaseCount[] _profiled;

    /// <summary>
    /// Initializes a new instance of the <see cref="MatrixAnalyzer"/> class.
    /// </summary>
    /// <param name="matrix">The regions of DNA to perform analyses on. All must have the same length.</param>
    /// <exception cref="ArgumentException">Thrown when the matrix is empty or the regions differ in length.</exception>
    public MatrixAnalyzer(IReadOnlyList<Sequence> matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        if (matrix.Count == 0)
        {
            throw new ArgumentException("The matrix must contain at least one region.", nameof(matrix));
        }

        _matrix = matrix;

        var columns = CountColumns(_matrix);
        _profiled = Profile(columns, _matrix.Count);

        Score = GetScore(columns);
        Entropy = columns.Sum(c => c.Entropy);
        Consensus = [.. columns.Select(c => c.Max[0].NBase)];
    }

    /// <summary>
    /// Gets the score of the regions matrix.
    /// The score shows how "conserved" the matrix of regions is: the lower the score, the more conserved the matrix is.
    /// Calculated as the sum of the differences between the total count of nucleotides in each column and the count of
    /// a nucleotide base that appears the most times in that column.
    /// For example, if the regions matrix has 10 regions and the first column is all A's, the score of the column is 0.
    /// If the second column has 7 A's, 2 C's, and 1 G, the score of the column is 3.
    /// </summary>
    public int Score { get; }

    /// <summary>
    /// Gets the consensus nucleotide bases of the regions matrix.
    /// The consensus nucleotide bases are the nucleotide bases that appear the most times in each column of the regions matrix.
    /// On a tie, the alphabetically first base is chosen.
    /// </summary>
    public NBase[] Consensus { get; }

    /// <summary>
    /// Gets the entropy, in bits, of the regions matrix.
    /// Entropy is a measure of the randomness of nucleotide bases in the matrix:
    /// the sum of the entropy of the observed base frequencies in each column. 0 means every column is fully conserved.
    /// </summary>
    public double Entropy { get; }

    /// <summary>
    /// Searches for regulatory motifs in the regions matrix with greedy motif search, using a profile with
    /// pseudocounts (Laplace's rule of succession). Returns one k-mer per region.
    /// </summary>
    /// <param name="k">The length of the motifs to search for.</param>
    /// <returns>The regulatory motifs found in the regions matrix.</returns>
    public IReadOnlyList<Sequence> SearchMotifs(int k)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(k);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(k, _matrix[0].Length);

        IReadOnlyList<Sequence> bestMotifs = [.. _matrix.Select(region => region.Part(0, k))];
        var bestScore = GetScore(CountColumns(bestMotifs));

        // Each region's k-mers, in order of first occurrence, so ties go to the leftmost.
        var candidates = _matrix
            .Select(region => region.FindAllKMers(k).Keys.ToArray())
            .ToArray();

        foreach (var first in candidates[0])
        {
            var motifs = new List<Sequence>(_matrix.Count) { first };

            for (var j = 1; j < _matrix.Count; j++)
            {
                var profile = Profile(CountColumns(motifs), motifs.Count);
                motifs.Add(GetMostProbable(candidates[j], profile));
            }

            var score = GetScore(CountColumns(motifs));
            if (score < bestScore)
            {
                bestMotifs = motifs;
                bestScore = score;
            }
        }

        return bestMotifs;
    }

    /// <summary>
    /// Gets the probability of a pattern under the profile of the regions matrix.
    /// </summary>
    /// <param name="pattern">The pattern to get the probability for. Must be as long as the regions.</param>
    /// <returns>The probability of the pattern under the profile of the regions matrix.</returns>
    public double Probability(IReadOnlyList<NBase> pattern) =>
        Probability(pattern, _profiled);

    /// <summary>
    /// Gets the probability of a pattern under the profile of the regions matrix.
    /// </summary>
    /// <param name="pattern">The pattern to get the probability for. Must be as long as the regions.</param>
    /// <returns>The probability of the pattern under the profile of the regions matrix.</returns>
    public double Probability(string pattern) =>
        Probability(new Sequence(pattern));

    /// <summary>
    /// Gets the most probable pattern from the given patterns, under the profile of the regions matrix.
    /// On a tie, the first pattern wins.
    /// </summary>
    /// <param name="patterns">The patterns to find the most probable from.</param>
    /// <returns>The most probable pattern.</returns>
    public T GetMostProbablePattern<T>(IReadOnlyList<T> patterns)
        where T : IReadOnlyList<NBase> =>
            GetMostProbable(patterns, _profiled);

    private static double Probability(IReadOnlyList<NBase> pattern, ProfiledNBaseCount[] profile)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        if (pattern.Count != profile.Length)
        {
            throw new ArgumentException($"The pattern must have {profile.Length} bases, the length of the profile.", nameof(pattern));
        }

        var probability = 1.0;
        for (var i = 0; i < pattern.Count; i++)
        {
            probability *= profile[i][pattern[i]];
        }

        return probability;
    }

    private static T GetMostProbable<T>(IReadOnlyList<T> patterns, ProfiledNBaseCount[] profile)
        where T : IReadOnlyList<NBase>
    {
        ArgumentNullException.ThrowIfNull(patterns);
        if (patterns.Count == 0)
        {
            throw new ArgumentException("At least one pattern is required.", nameof(patterns));
        }

        var mostProbable = patterns[0];
        var maxProbability = Probability(mostProbable, profile);

        for (var i = 1; i < patterns.Count; i++)
        {
            var probability = Probability(patterns[i], profile);
            if (probability > maxProbability)
            {
                maxProbability = probability;
                mostProbable = patterns[i];
            }
        }

        return mostProbable;
    }

    /// <summary>
    /// Counts the number of occurrences of each <see cref="NBase"/> in each column of the given regions.
    /// </summary>
    /// <param name="regions">The regions to count columns in.</param>
    /// <returns><see cref="NBaseCount"/>: The count of each nucleotide base in each column.</returns>
    /// <exception cref="ArgumentException">Thrown when the regions are not all the same length.</exception>
    private static NBaseCount[] CountColumns(IReadOnlyList<Sequence> regions)
    {
        var regionLength = regions[0].Length;

        var columns = new NBase[regionLength][];
        for (var i = 0; i < regionLength; i++)
        {
            columns[i] = new NBase[regions.Count];
        }

        for (var i = 0; i < regions.Count; i++)
        {
            if (regions[i].Length != regionLength)
            {
                throw new ArgumentException("All regions must be the same length.", nameof(regions));
            }

            for (var j = 0; j < regionLength; j++)
            {
                columns[j][i] = regions[i][j];
            }
        }

        return [.. columns.Select(c => new NBaseCount(c))];
    }

    private static ProfiledNBaseCount[] Profile(NBaseCount[] columns, int regionCount) =>
        [.. columns.Select(c => c.Normalize(regionCount))];

    /// <summary>
    /// Gets the score of the regions matrix.
    /// Calculated as the sum of the differences between the total count of nucleotide bases
    /// in each column and the count of a nucleotide base that appears the most times in that column.
    /// </summary>
    /// <param name="columns">The columns to get the score for.</param>
    /// <returns>The score of the regions matrix.</returns>
    private static int GetScore(NBaseCount[] columns) =>
        columns.Sum(column => column.Total - column.Max[0].Count);
}
