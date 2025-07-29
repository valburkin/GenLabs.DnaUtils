namespace GenLabs.DnaUtils;

/// <summary>
/// Represents a profiled <see cref="NBaseCount"/>.
/// Profiled means that the count of each <see cref="NBase"/> has been normalized by a factor.
/// Uses Laplace's rule of succession to normalize the count of each <see cref="NBase"/>.
/// </summary>
public readonly struct ProfiledNBaseCount
{
    private const double EqualityTolerance = .001;

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
        _a = (nBaseCount[NBase.A] + 1) / ((double) factor + 4);
        _c = (nBaseCount[NBase.C] + 1) / ((double) factor + 4);
        _g = (nBaseCount[NBase.G] + 1) / ((double) factor + 4);
        _t = (nBaseCount[NBase.T] + 1) / ((double) factor + 4);

        Max = GetMax();
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
        _ => throw new ArgumentOutOfRangeException(nameof(nBase))
    };

    /// <summary>
    /// Gets the <see cref="NBase"/>s with the maximum value.
    /// </summary>
    /// <returns>The <see cref="NBase"/>s with the maximum value and their corresponding value.</returns>
    public (NBase NBase, double Value)[] Max { get; }

    /// <summary>
    /// Gets the entropy of the <see cref="ProfiledNBaseCount"/>.
    /// Entropy is a measure of the randomness of nucleotide bases in the set.
    /// It is calculated as the sum of the negative of the product of the value of each <see cref="NBase"/> and the logarithm base 2 of the value of each nucleotide base.
    /// </summary>
    /// <returns>The entropy of the <see cref="ProfiledNBaseCount"/>.</returns>
    public double Entropy { get; }

    /// <summary>
    /// Gets the factor used to normalize the <see cref="NBaseCount"/>.
    /// </summary>
    public int Factor { get; }

    private (NBase, double)[] GetMax()
    {
        var max = new List<(NBase n, double c)> { (NBase.A, _a) };

        foreach (var n in new[] { NBase.C, NBase.G, NBase.T })
        {
            var v = this[n];

            if (Math.Abs(v - max[0].c) < EqualityTolerance)
            {
                max.Add((n, v));
            }
            else if (v > max[0].c)
            {
                max.Clear();
                max.Add((n, v));
            }
        }

        return max.ToArray();
    }
    
    private double GetEntropy() => -new[] { _a, _c, _g, _t }
        .Where(v => v > 0)
        .Select(v => v * Math.Log2(v))
        .Sum();
}
