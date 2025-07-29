namespace GenLabs.DnaUtils;

/// <summary>
/// Represents a nucleotide base count.
/// </summary>
public readonly struct NBaseCount
{
    private readonly int _a;
    private readonly int _c;
    private readonly int _g;
    private readonly int _t;

    /// <summary>
    /// Initializes a new instance of the <see cref="NBaseCount"/> struct.
    /// </summary>
    /// <param name="nBases">The nucleotide bases to count.</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when an invalid nucleotide base is encountered.</exception>
    public NBaseCount(IReadOnlyCollection<NBase> nBases)
    {
        foreach (var nBase in nBases)
        {
            switch (nBase)
            {
                case NBase.A:
                    _a++;
                    break;
                case NBase.C:
                    _c++;
                    break;
                case NBase.G:
                    _g++;
                    break;
                case NBase.T:
                    _t++;
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(nBase));
            }
        }

        Total = nBases.Count;
        Max = GetMax();
    }

    /// <summary>
    /// Gets the count of a specific <see cref="NBase"/>.
    /// </summary>
    /// <param name="nBase">The <see cref="NBase"/> to get the count for.</param>
    /// <returns>The count of the <see cref="NBase"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when an invalid nucleotide base is encountered.</exception>
    public int this[NBase nBase] => nBase switch
    {
        NBase.A => _a,
        NBase.C => _c,
        NBase.G => _g,
        NBase.T => _t,
        _ => throw new ArgumentOutOfRangeException(nameof(nBase))
    };

    /// <summary>
    /// Gets the total count of <see cref="NBase"/>s.
    /// </summary>
    public int Total { get; }

    /// <summary>
    /// Gets the maximum <see cref="NBase"/> counts.
    /// </summary>
    /// <returns>The <see cref="NBase"/>s that appears the most times and its corresponding counts.</returns>
    public (NBase NBase, int Count)[] Max { get; }

    /// <summary>
    /// Normalizes the <see cref="NBaseCount"/> by a factor.
    /// </summary>
    /// <param name="factor">The factor to normalize by.</param>
    /// <returns>The normalized nucleotide base count.</returns>
    public ProfiledNBaseCount Normalize(int factor) => new(this, factor);

    private (NBase, int)[] GetMax()
    {
        var max = new List<(NBase n, int c)> { (NBase.A, _a) };

        foreach (var n in new[] { NBase.C, NBase.G, NBase.T })
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
