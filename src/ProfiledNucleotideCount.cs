namespace GenLabs.DnaUtils;

/// <summary>
/// Represents a profiled <see cref="NucleotideCount"/>.
/// Profiled means that the count of each <see cref="Nucleotide"/> has been normalized by a factor.
/// Uses Laplace's rule of succession to normalize the count of each <see cref="Nucleotide"/>.
/// </summary>
public readonly struct ProfiledNucleotideCount
{
    private const double EqualityTolerance = .001;

    private readonly double _a;
    private readonly double _c;
    private readonly double _g;
    private readonly double _t;

    /// <summary>
    /// Initializes a new instance of the <see cref="ProfiledNucleotideCount"/> struct.
    /// </summary>
    /// <param name="nucleotideCount">The <see cref="NucleotideCount"/> to normalize.</param>
    /// <param name="factor">The factor to normalize by.</param>
    public ProfiledNucleotideCount(NucleotideCount nucleotideCount, int factor)
    {
        _a = (nucleotideCount[Nucleotide.A] + 1) / ((double) factor + 4);
        _c = (nucleotideCount[Nucleotide.C] + 1) / ((double) factor + 4);
        _g = (nucleotideCount[Nucleotide.G] + 1) / ((double) factor + 4);
        _t = (nucleotideCount[Nucleotide.T] + 1) / ((double) factor + 4);

        Max = GetMax();
        Entropy = GetEntropy();
        Factor = factor;
    }

    /// <summary>
    /// Gets the normalized count of a specific <see cref="Nucleotide"/>.
    /// </summary>
    /// <param name="nucleotide">The <see cref="Nucleotide"/> to get the normalized count for.</param>
    /// <returns>The normalized count of the <see cref="Nucleotide"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when an invalid nucleotide is encountered.</exception>
    public double this[Nucleotide nucleotide] => nucleotide switch
    {
        Nucleotide.A => _a,
        Nucleotide.C => _c,
        Nucleotide.G => _g,
        Nucleotide.T => _t,
        _ => throw new ArgumentOutOfRangeException(nameof(nucleotide))
    };

    /// <summary>
    /// Gets the <see cref="Nucleotide"/>s with the maximum value.
    /// </summary>
    /// <returns>The <see cref="Nucleotide"/>s with the maximum value and their corresponding value.</returns>
    public (Nucleotide Nucleotide, double Value)[] Max { get; }

    /// <summary>
    /// Gets the entropy of the <see cref="ProfiledNucleotideCount"/>.
    /// Entropy is a measure of the randomness of nucleotides in the set.
    /// It is calculated as the sum of the negative of the product of the value of each <see cref="Nucleotide"/> and the logarithm base 2 of the value of each nucleotide.
    /// </summary>
    /// <returns>The entropy of the <see cref="ProfiledNucleotideCount"/>.</returns>
    public double Entropy { get; }

    /// <summary>
    /// Gets the factor used to normalize the <see cref="NucleotideCount"/>.
    /// </summary>
    public int Factor { get; }

    private (Nucleotide, double)[] GetMax()
    {
        var max = new List<(Nucleotide n, double c)> { (Nucleotide.A, _a) };

        foreach (var n in new[] { Nucleotide.C, Nucleotide.G, Nucleotide.T })
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