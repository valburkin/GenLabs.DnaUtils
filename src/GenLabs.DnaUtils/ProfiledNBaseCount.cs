namespace GenLabs.DnaUtils;

/// <summary>
/// Represents a profiled <see cref="NBaseCount"/>.
/// Profiled means that the count of each <see cref="NBase"/> has been normalized by a factor.
/// Uses Laplace's rule of succession (one pseudocount per base), so no base ever has probability zero.
/// </summary>
public readonly struct ProfiledNBaseCount
{
    private readonly double _a;
    private readonly double _c;
    private readonly double _g;
    private readonly double _t;

    /// <summary>
    /// Initializes a new instance of the <see cref="ProfiledNBaseCount"/> struct.
    /// </summary>
    /// <param name="nBaseCount">The <see cref="NBaseCount"/> to normalize.</param>
    /// <param name="factor">The factor to normalize by.</param>
    public ProfiledNBaseCount(NBaseCount nBaseCount, int factor)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(factor);

        var denominator = (double)factor + 4;
        _a = (nBaseCount[NBase.A] + 1) / denominator;
        _c = (nBaseCount[NBase.C] + 1) / denominator;
        _g = (nBaseCount[NBase.G] + 1) / denominator;
        _t = (nBaseCount[NBase.T] + 1) / denominator;

        // Every value shares the denominator, so the most probable bases are the most counted ones.
        var maxValue = (nBaseCount.Max[0].Count + 1) / denominator;
        Max = [.. nBaseCount.Max.Select(m => (m.NBase, maxValue))];
        Entropy = GetEntropy();
        Factor = factor;
    }

    /// <summary>
    /// Gets the normalized count of a specific <see cref="NBase"/>.
    /// </summary>
    /// <param name="nBase">The <see cref="NBase"/> to get the normalized count for.</param>
    /// <returns>The normalized count of the <see cref="NBase"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when an invalid nucleotide base is encountered.</exception>
    public double this[NBase nBase] => nBase switch
    {
        NBase.A => _a,
        NBase.C => _c,
        NBase.G => _g,
        NBase.T => _t,
        _ => throw new ArgumentOutOfRangeException(nameof(nBase)),
    };

    /// <summary>
    /// Gets the <see cref="NBase"/>s with the maximum value.
    /// </summary>
    /// <returns>The <see cref="NBase"/>s with the maximum value and their corresponding value.</returns>
    public (NBase NBase, double Value)[] Max { get; }

    /// <summary>
    /// Gets the entropy, in bits, of the profiled (pseudocounted) values.
    /// Calculated as the sum of -p log2(p) over the four bases.
    /// For the entropy of the observed frequencies, see <see cref="NBaseCount.Entropy"/>.
    /// </summary>
    public double Entropy { get; }

    /// <summary>
    /// Gets the factor used to normalize the <see cref="NBaseCount"/>.
    /// </summary>
    public int Factor { get; }

    private double GetEntropy() => -new[] { _a, _c, _g, _t }
        .Select(v => v * Math.Log2(v))
        .Sum();
}
