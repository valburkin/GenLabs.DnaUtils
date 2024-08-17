namespace GenLabs.DnaUtils;

/// <summary>
/// Represents a nucleotide count.
/// </summary>
public readonly struct NucleotideCount
{
    private readonly int _a;
    private readonly int _c;
    private readonly int _g;
    private readonly int _t;

    /// <summary>
    /// Gets the total count of <see cref="Nucleotide"/>s.
    /// </summary>
    public int Total => _a + _c + _g + _t;

    /// <summary>
    /// Initializes a new instance of the <see cref="NucleotideCount"/> struct.
    /// </summary>
    /// <param name="nucleotides">The nucleotides to count.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when an invalid nucleotide is encountered.</exception>
    public NucleotideCount(IReadOnlyCollection<Nucleotide> nucleotides)
    {
        foreach (var nucleotide in nucleotides)
        {
            switch (nucleotide)
            {
                case Nucleotide.A:
                    _a++;
                    break;
                case Nucleotide.C:
                    _c++;
                    break;
                case Nucleotide.G:
                    _g++;
                    break;
                case Nucleotide.T:
                    _t++;
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(nucleotide));
            }
        }
    }

    /// <summary>
    /// Gets the count of a specific <see cref="Nucleotide"/>.
    /// </summary>
    /// <param name="nucleotide">The <see cref="Nucleotide"/> to get the count for.</param>
    /// <returns>The count of the <see cref="Nucleotide"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when an invalid nucleotide is encountered.</exception>
    public int this[Nucleotide nucleotide] => nucleotide switch
    {
        Nucleotide.A => _a,
        Nucleotide.C => _c,
        Nucleotide.G => _g,
        Nucleotide.T => _t,
        _ => throw new ArgumentOutOfRangeException(nameof(nucleotide))
    };

    /// <summary>
    /// Normalizes the <see cref="NucleotideCount"/> by a factor.
    /// </summary>
    /// <param name="factor">The factor to normalize by.</param>
    /// <returns>The normalized nucleotide count.</returns>
    public NormalizedNucleotideCount Normalize(int factor) => new(this, factor);

    /// <summary>
    /// Gets the maximum <see cref="Nucleotide"/> counts.
    /// </summary>
    /// <returns>The <see cref="Nucleotide"/>s that appears the most times and its corresponding counts.</returns>
    public (Nucleotide Nucleotide, int Count)[] Max()
    {
        var max = new List<(Nucleotide n, int c)> { (Nucleotide.A, _a) };

        foreach (var n in new []{ Nucleotide.C, Nucleotide.G, Nucleotide.T })
        {
            var c = this[n];

            if (c > max[0].c)
            {
                max.Clear();
                max.Add((n, c));
            }
            else if (c == max[0].c)
            {
                max.Add((n, c));
            }
        }

        return max.ToArray();
    }
}