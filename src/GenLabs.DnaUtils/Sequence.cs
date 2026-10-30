using System.Collections;

namespace GenLabs.DnaUtils;

/// <summary>
/// Represents an immutable sequence of <see cref="NBase"/>s.
/// </summary>
public sealed class Sequence : IEquatable<Sequence>, IReadOnlyList<NBase>
{
    /// <summary>
    /// The longest length <see cref="AllOfLength"/> accepts: 4^15 is about a billion sequences.
    /// </summary>
    public const int MaxEnumerableLength = 15;

    private readonly NBase[] _sequence;

    /// <summary>
    /// Instantiates a <see cref="Sequence"/> from a list of <see cref="NBase"/>s.
    /// </summary>
    /// <param name="sequence">The list of <see cref="NBase"/>s. It is copied.</param>
    public Sequence(IReadOnlyList<NBase> sequence)
        : this(sequence?.ToArray() ?? throw new ArgumentNullException(nameof(sequence)))
    { }

    /// <summary>
    /// Instantiates a <see cref="Sequence"/> from a string such as <c>"ACGT"</c>. Lowercase is accepted.
    /// </summary>
    /// <param name="sequence">The string to convert to a <see cref="Sequence"/>.</param>
    /// <exception cref="ArgumentException">Thrown when the string contains a character other than A, C, G or T.</exception>
    public Sequence(string sequence)
        : this(Parse(sequence))
    { }

    // Takes ownership of the array without copying it. Only this class calls it, always with an array
    // nothing else holds; callers outside the class bind to the public constructor, which copies.
    private Sequence(NBase[] owned) => _sequence = owned;

    /// <summary>
    /// The <see cref="NBaseCount"/> of the <see cref="Sequence"/>.
    /// </summary>
    public NBaseCount Counts => new(_sequence);

    /// <summary>
    /// The length of the <see cref="Sequence"/>.
    /// </summary>
    public int Length => _sequence.Length;

    /// <summary>
    /// Same as <see cref="Length"/>. The number of <see cref="NBase"/>s in the <see cref="Sequence"/>.
    /// </summary>
    public int Count => _sequence.Length;

    /// <summary>
    /// Returns the <see cref="NBase"/> at the specified index.
    /// </summary>
    /// <param name="index">The index of the <see cref="NBase"/>.</param>
    /// <returns>The <see cref="NBase"/> at the specified index.</returns>
    public NBase this[int index] => _sequence[index];

    /// <summary>
    /// Adds a <see cref="NBase"/> to the end of the sequence.
    /// </summary>
    /// <param name="sequence">The <see cref="Sequence"/> to add to.</param>
    /// <param name="nBase">The <see cref="NBase"/> to add.</param>
    /// <returns>New <see cref="Sequence"/></returns>
    public static Sequence operator +(Sequence sequence, NBase nBase)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        return new Sequence([.. sequence._sequence, nBase]);
    }

    /// <summary>
    /// Adds a <see cref="NBase"/> to the beginning of the sequence.
    /// </summary>
    /// <param name="nBase">The <see cref="NBase"/> to add.</param>
    /// <param name="sequence">The <see cref="Sequence"/> to add to.</param>
    /// <returns>New <see cref="Sequence"/></returns>
    public static Sequence operator +(NBase nBase, Sequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        return new Sequence([nBase, .. sequence._sequence]);
    }

    /// <summary>
    /// Instantiates a <see cref="Sequence"/> from a string.
    /// </summary>
    /// <param name="sequence">The string to convert to a <see cref="Sequence"/>.</param>
    public static implicit operator Sequence(string sequence) => new(sequence);

    /// <summary>
    /// Instantiates a <see cref="Sequence"/> from a string. Same as the implicit conversion.
    /// </summary>
    /// <param name="sequence">The string to convert to a <see cref="Sequence"/>.</param>
    public static Sequence FromString(string sequence) => new(sequence);

    public static bool operator ==(Sequence? left, Sequence? right) =>
        left is null ? right is null : left.Equals(right);

    public static bool operator !=(Sequence? left, Sequence? right) => !(left == right);

    /// <summary>
    /// Generates all possible <see cref="Sequence"/>s of a given length, in lexicographic order.
    /// IMPORTANT: the number of sequences grows exponentially with the length: 4^length.
    /// For length 3 there are 64 sequences; for length 10 there are 1,048,576.
    /// </summary>
    /// <param name="length">The length of the sequences, from 0 to <see cref="MaxEnumerableLength"/>.</param>
    /// <returns>All possible <see cref="Sequence"/>s of the given length.</returns>
    public static IEnumerable<Sequence> AllOfLength(int length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(length, MaxEnumerableLength);

        return Enumerate(length);

        static IEnumerable<Sequence> Enumerate(int length)
        {
            var count = 1 << (2 * length);
            for (var i = 0; i < count; i++)
            {
                var sequence = new NBase[length];
                for (var j = 0; j < length; j++)
                {
                    // The first base is the most significant pair of bits.
                    sequence[j] = (NBase)((i >> (2 * (length - j - 1))) & 3);
                }

                yield return new Sequence(sequence);
            }
        }
    }

    /// <summary>
    /// Returns a part of the <see cref="Sequence"/>.
    /// </summary>
    /// <param name="start">The start index of the part.</param>
    /// <param name="length">The length of the part.</param>
    /// <returns>New <see cref="Sequence"/></returns>
    public Sequence Part(int start, int length) =>
        new(_sequence[start..checked(start + length)]);

    /// <summary>
    /// Returns the complement of the sequence.
    /// The complement is the sequence where each <see cref="NBase"/> is replaced by its complement.
    /// </summary>
    /// <returns>New <see cref="Sequence"/></returns>
    public Sequence Complement()
    {
        var complement = new NBase[Length];
        for (var i = 0; i < Length; i++)
        {
            complement[i] = _sequence[i].Complement();
        }

        return new Sequence(complement);
    }

    /// <summary>
    /// Returns the reverse of the <see cref="Sequence"/>.
    /// </summary>
    /// <returns>New <see cref="Sequence"/></returns>
    public Sequence Reverse()
    {
        var reverse = new NBase[Length];
        for (var i = 0; i < Length; i++)
        {
            reverse[i] = _sequence[Length - i - 1];
        }

        return new Sequence(reverse);
    }

    /// <summary>
    /// Returns the reverse complement of the <see cref="Sequence"/>: the opposite strand, read in its own 5' to 3' direction.
    /// </summary>
    /// <returns>New <see cref="Sequence"/></returns>
    public Sequence ReverseComplement()
    {
        var reverseComplement = new NBase[Length];
        for (var i = 0; i < Length; i++)
        {
            reverseComplement[i] = _sequence[Length - i - 1].Complement();
        }

        return new Sequence(reverseComplement);
    }

    /// <summary>
    /// Returns every k-mer of the <see cref="Sequence"/> and the positions where it starts, in ascending order.
    /// The k-mers appear in the order of their first occurrence.
    /// </summary>
    /// <param name="k">The length of the k-mer.</param>
    /// <returns>A dictionary of k-mers and their positions in the <see cref="Sequence"/>.</returns>
    public Dictionary<Sequence, List<int>> FindAllKMers(int k)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(k);

        var kMers = new Dictionary<Sequence, List<int>>();

        for (var i = 0; i <= _sequence.Length - k; i++)
        {
            Add(kMers, Part(i, k), i);
        }

        return kMers;
    }

    /// <summary>
    /// Returns every k-mer that occurs in the <see cref="Sequence"/> with at most <paramref name="mismatches"/> mismatches,
    /// and the positions where it occurs. A k-mer need not appear exactly to be included.
    /// The number of positions is the k-mer's approximate frequency, so the most frequent k-mers solve the
    /// "frequent words with mismatches" problem.
    /// </summary>
    /// <param name="k">The length of the k-mer.</param>
    /// <param name="mismatches">The maximum number of mismatches allowed.</param>
    /// <param name="treatReverseComplementsAsEqual">
    /// Whether an approximate occurrence of a k-mer's reverse complement also counts as an occurrence of the k-mer,
    /// as when searching both strands of DNA.
    /// </param>
    /// <returns>
    /// A dictionary of k-mers and the positions, in ascending order, where they occur approximately.
    /// With <paramref name="treatReverseComplementsAsEqual"/>, a position is listed twice when both the k-mer and its
    /// reverse complement occur there, so the count of positions is the k-mer's frequency on both strands.
    /// </returns>
    public Dictionary<Sequence, List<int>> FindAllKMers(int k, int mismatches, bool treatReverseComplementsAsEqual)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(k);
        ArgumentOutOfRangeException.ThrowIfNegative(mismatches);

        var kMers = new Dictionary<Sequence, List<int>>();

        for (var i = 0; i <= _sequence.Length - k; i++)
        {
            foreach (var neighbor in Part(i, k).Neighbors(mismatches))
            {
                Add(kMers, neighbor, i);

                if (treatReverseComplementsAsEqual)
                {
                    Add(kMers, neighbor.ReverseComplement(), i);
                }
            }
        }

        return kMers;
    }

    /// <summary>
    /// Finds clumping k-mers: k-mers that occur at least <paramref name="threshold"/> times
    /// within some window of <paramref name="windowSize"/> bases.
    /// </summary>
    /// <param name="windowSize">The size of the window.</param>
    /// <param name="k">The length of the k-mer.</param>
    /// <param name="threshold">The minimum number of times the k-mer must appear in the window.</param>
    /// <returns>The clumping k-mers, in the order of their first occurrence.</returns>
    public Sequence[] FindClumpingKMers(int windowSize, int k, int threshold)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(threshold);

        var clumpingKMers = new List<Sequence>();
        foreach (var (kMer, positions) in FindAllKMers(k))
        {
            for (var i = 0; i + threshold - 1 < positions.Count; i++)
            {
                var start = positions[i];
                var end = positions[i + threshold - 1] + k;

                if (end - start <= windowSize)
                {
                    clumpingKMers.Add(kMer);
                    break;
                }
            }
        }

        return [.. clumpingKMers];
    }

    /// <summary>
    /// Returns the cumulative GC skew: element i is the count of <see cref="NBase.G"/> minus the count of
    /// <see cref="NBase.C"/> in the first i bases. The array has <see cref="Length"/> + 1 elements and starts with 0.
    /// </summary>
    /// <returns>The skew after each prefix of the sequence.</returns>
    public int[] GetSkew()
    {
        var skew = new int[Length + 1];
        for (var i = 1; i <= Length; i++)
        {
            skew[i] = skew[i - 1] + _sequence[i - 1] switch
            {
                NBase.G => 1,
                NBase.C => -1,
                _ => 0,
            };
        }

        return skew;
    }

    /// <summary>
    /// Returns the positions where the skew (see <see cref="GetSkew"/>) is minimal.
    /// In a bacterial genome these point to the likely origin of replication.
    /// </summary>
    /// <returns>An array of positions where the skew is minimal, in ascending order.</returns>
    public int[] GetMinSkewPositions()
    {
        var skew = GetSkew();
        var min = skew.Min();

        return skew
            .Select((value, index) => (value, index))
            .Where(x => x.value == min)
            .Select(x => x.index)
            .ToArray();
    }

    /// <summary>
    /// Returns the Hamming distance between two sequences: the number of positions at which their bases differ.
    /// </summary>
    /// <param name="other">The other <see cref="Sequence"/>.</param>
    /// <returns>The Hamming distance between the two sequences.</returns>
    /// <exception cref="ArgumentException">Thrown when the sequences have different lengths.</exception>
    public int HammingDistance(Sequence other)
    {
        ArgumentNullException.ThrowIfNull(other);

        if (Length != other.Length)
        {
            throw new ArgumentException("The sequences must have the same length.", nameof(other));
        }

        return HammingDistance(_sequence, other._sequence);
    }

    /// <summary>
    /// Finds the positions where <paramref name="pattern"/> occurs with at most <paramref name="mismatches"/> mismatches.
    /// </summary>
    /// <param name="pattern">The pattern <see cref="Sequence"/> to find.</param>
    /// <param name="mismatches">The maximum number of mismatches allowed.</param>
    /// <returns>An array of positions, in ascending order, where the pattern is found.</returns>
    public int[] FindPatternLocations(Sequence pattern, int mismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentOutOfRangeException.ThrowIfNegative(mismatches);

        var positions = new List<int>();
        for (var i = 0; i <= Length - pattern.Length; i++)
        {
            if (HammingDistance(_sequence.AsSpan(i, pattern.Length), pattern._sequence) <= mismatches)
            {
                positions.Add(i);
            }
        }

        return [.. positions];
    }

    /// <summary>
    /// Returns the d-neighborhood of the <see cref="Sequence"/>: every sequence of the same length
    /// with at most <paramref name="mismatches"/> mismatches, the sequence itself included.
    /// </summary>
    /// <param name="mismatches">The maximum number of mismatches allowed.</param>
    /// <returns>The neighbors of the <see cref="Sequence"/>.</returns>
    public Sequence[] Neighbors(int mismatches)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(mismatches);

        if (mismatches == 0 || Length == 0)
        {
            return [this];
        }

        if (Length == 1)
        {
            return [.. NBaseHelper.All.Select(n => new Sequence([n]))];
        }

        var suffix = Part(1, Length - 1);

        var neighbors = new List<Sequence>();
        foreach (var neighbor in suffix.Neighbors(mismatches))
        {
            if (suffix.HammingDistance(neighbor) < mismatches)
            {
                foreach (var nBase in NBaseHelper.All)
                {
                    neighbors.Add(nBase + neighbor);
                }
            }
            else
            {
                neighbors.Add(this[0] + neighbor);
            }
        }

        return [.. neighbors];
    }

    /// <summary>
    /// Returns a copy of the <see cref="Sequence"/> with the base at <paramref name="index"/> replaced.
    /// </summary>
    /// <param name="index">The index to mutate.</param>
    /// <param name="nBase">The <see cref="NBase"/> to mutate to.</param>
    /// <returns>Mutated <see cref="Sequence"/></returns>
    public Sequence Mutate(int index, NBase nBase)
    {
        var mutated = (NBase[])_sequence.Clone();
        mutated[index] = nBase;

        return new Sequence(mutated);
    }

    public override string ToString() =>
        string.Create(Length, _sequence, static (chars, sequence) =>
        {
            for (var i = 0; i < sequence.Length; i++)
            {
                chars[i] = sequence[i].ToChar();
            }
        });

    public bool Equals(Sequence? other) =>
        other is not null && _sequence.AsSpan().SequenceEqual(other._sequence);

    public override bool Equals(object? obj) => Equals(obj as Sequence);

    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.AddBytes(System.Runtime.InteropServices.MemoryMarshal.AsBytes(_sequence.AsSpan()));
        return hash.ToHashCode();
    }

    public IEnumerator<NBase> GetEnumerator() => ((IEnumerable<NBase>)_sequence).GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    private static void Add(Dictionary<Sequence, List<int>> positionsByKMer, Sequence kMer, int position)
    {
        if (positionsByKMer.TryGetValue(kMer, out var positions))
        {
            positions.Add(position);
        }
        else
        {
            positionsByKMer[kMer] = [position];
        }
    }

    private static int HammingDistance(ReadOnlySpan<NBase> left, ReadOnlySpan<NBase> right)
    {
        var distance = 0;
        for (var i = 0; i < left.Length; i++)
        {
            if (left[i] != right[i])
            {
                distance++;
            }
        }

        return distance;
    }

    private static NBase[] Parse(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var nBases = new NBase[sequence.Length];
        for (var i = 0; i < sequence.Length; i++)
        {
            nBases[i] = sequence[i].ToNBase();
        }

        return nBases;
    }
}
