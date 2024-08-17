using System.Text;

namespace GenLabs.DnaUtils;

/// <summary>
/// Represents a sequence of <see cref="Nucleotide"/>s.
/// </summary>
public sealed class Sequence : IEquatable<Sequence>
{
    private readonly Nucleotide[] _sequence;

    /// <summary>
    /// The <see cref="NucleotideCount"/> of the <see cref="Sequence"/>.
    /// </summary>
    public NucleotideCount Counts { get; }

    /// <summary>
    /// The length of the <see cref="Sequence"/>.
    /// </summary>
    public int Length => _sequence.Length;

    /// <summary>
    /// Instantiates a <see cref="Sequence"/> from an array of <see cref="Nucleotide"/>s.
    /// </summary>
    /// <param name="sequence">The array of <see cref="Nucleotide"/>s.</param>
    public Sequence(IReadOnlyList<Nucleotide> sequence)
    {
        _sequence = sequence.ToArray();
        Counts = new(sequence);
    }

    /// <summary>
    /// Instantiates a <see cref="Sequence"/> from a string.    
    /// </summary>
    /// <param name="sequence">The string to convert to a <see cref="Sequence"/>.</param>
    public Sequence(string sequence)
        : this(GetNucleotides(sequence))
    { }

    /// <summary>
    /// Returns a part of the <see cref="Sequence"/>>.
    /// </summary>
    /// <param name="start">The start index of the part.</param>
    /// <param name="length">The length of the part.</param>
    /// <returns>New <see cref="Sequence"/></returns>
    public Sequence Part(int start, int length) =>
        new(_sequence[start..checked(start + length)]);

    /// <summary>
    /// Returns the <see cref="Nucleotide"/> at the specified index.
    /// </summary>
    /// <param name="i">The index of the <see cref="Nucleotide"/>.</param>
    /// <returns>The <see cref="Nucleotide"/> at the specified index.</returns>
    public Nucleotide this[int i] => _sequence[i];

    /// <summary>
    /// Adds a <see cref="Nucleotide"/> to the end of the sequence.
    /// </summary>
    /// <param name="sequence">The <see cref="Sequence"/> to add to.</param>
    /// <param name="nucleotide">The <see cref="Nucleotide"/> to add.</param>
    /// <returns>New <see cref="Sequence"/></returns>
    public static Sequence operator +(Sequence sequence, Nucleotide nucleotide)
    {
        var newSequence = new Nucleotide[sequence.Length + 1];
        Array.Copy(sequence._sequence, newSequence, sequence.Length);
        newSequence[sequence.Length] = nucleotide;

        return new Sequence(newSequence);
    }

    /// <summary>
    /// Add a <see cref="Nucleotide"/> to the beginning of the sequence.
    /// </summary>
    /// <param name="nucleotide">The <see cref="Nucleotide"/> to add.</param>
    /// <param name="sequence">The <see cref="Sequence"/> to add to.</param>
    /// <returns>New <see cref="Sequence"/></returns>
    public static Sequence operator +(Nucleotide nucleotide, Sequence sequence)
    {
        var newSequence = new Nucleotide[sequence.Length + 1];
        newSequence[0] = nucleotide;
        Array.Copy(sequence._sequence, 0, newSequence, 1, sequence.Length);

        return new Sequence(newSequence);
    }

    /// <summary>
    /// Instantiates a <see cref="Sequence"/> from a string.
    /// </summary>
    /// <param name="sequence">The string to convert to a <see cref="Sequence"/>.</param>
    public static implicit operator Sequence(string sequence) => new(sequence);

    /// <summary>
    /// Generates all possible <see cref="Sequence"/>s of a given length.
    /// IMPORTANT: Note that possible number of sequences grows exponentially with the length: 4^length.
    /// For example, for length 3, there are 64 possible sequences.
    /// For length 10, there are 1,048,576 possible sequences.
    /// </summary>
    /// <param name="length">The length of the sequences.</param>
    /// <returns>All possible <see cref="Sequence"/>s of the given length.</returns>
    public static IEnumerable<Sequence> AllMutations(int length)
    {
        var numberOfMutations = 1 << (2 * length);
        for (var i = 0; i < numberOfMutations; i++)
        {
            var sequence = new Nucleotide[length];
            for (var j = 0; j < length; j++)
            {
                sequence[j] = NucleotideHelper.All[(i >> (2 * j)) & 3];
            }

            yield return new Sequence(sequence);
        }
    }

    /// <summary>
    /// Returns the complement of the sequence.
    /// The complement is the sequence where each <see cref="Nucleotide"/> is replaced by its complement.
    /// </summary>
    /// <returns>New <see cref="Sequence"/></returns>
    public Sequence Complement()
    {
        var complement = new Nucleotide[Length];
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
        var reverse = new Nucleotide[Length];
        for (var i = 0; i < Length; i++)
        {
            reverse[i] = _sequence[Length - i - 1];
        }

        return new Sequence(reverse);
    }

    /// <summary>
    /// Returns the reverse complement of the <see cref="Sequence"/>.
    /// The reverse complement is the sequence where each <see cref="Nucleotide"/> is replaced by its complement from the end to the beginning.
    /// </summary>
    /// <returns>New <see cref="Sequence"/></returns>
    public Sequence ReverseComplement()
    {
        var reverseComplement = new Nucleotide[Length];
        for (var i = 0; i < Length; i++)
        {
            reverseComplement[i] = _sequence[Length - i - 1].Complement();
        }

        return new Sequence(reverseComplement);
    }

    /// <summary>
    /// Returns a dictionary of k-mers and their locations in the <see cref="Sequence"/>.
    /// </summary>
    /// <param name="k">The length of the k-mer.</param>
    /// <returns>A dictionary of k-mers and their locations in the <see cref="Sequence"/>.</returns>
    public Dictionary<Sequence, List<int>> FindAllKMers(int k)
    {
        var counts = new Dictionary<Sequence, List<int>>(new SequenceEqualityComparer());

        for (var i = 0; i <= _sequence.Length - k; i++)
        {
            var subsequence = Part(i, k);

            if (counts.TryGetValue(subsequence, out var positions))
            {
                positions.Add(i);
            }
            else
            {
                counts[subsequence] = [i];
            }
        }

        return counts;
    }

    /// <summary>
    /// Returns a dictionary of k-mers and their locations in the <see cref="Sequence"/>.
    /// </summary>
    /// <param name="k">The length of the k-mer.</param>
    /// <param name="mismatches">The maximum number of mismatches allowed.</param>
    /// <param name="treatReverseComplementsAsEqual">Whether to consider reverse complements.</param>
    /// <returns>A dictionary of k-mers and their locations in the <see cref="Sequence"/>.</returns>
    public Dictionary<Sequence, List<int>> FindAllKMers(int k, int mismatches, bool treatReverseComplementsAsEqual)
    {
        var allKMers = FindAllKMers(k);

        if (treatReverseComplementsAsEqual)
        {
            var reversedKMers = new Dictionary<Sequence, List<int>>(new SequenceEqualityComparer());
            foreach (var (kMer, positions) in allKMers)
            {
                var reverseComplement = kMer.ReverseComplement();
                reversedKMers[reverseComplement] = positions;
            }

            foreach (var (reversedKMer, positions) in reversedKMers)
            {
                allKMers[reversedKMer] = positions;
            }
        }

        if (mismatches > 0)
        {
            var wobbledKMers = new Dictionary<Sequence, List<int>>(new SequenceEqualityComparer());
            foreach (var (kMer, positions) in allKMers)
            {
                var wobbles = kMer.Wobbles(mismatches);
                foreach (var wobble in wobbles)
                {
                    wobbledKMers[wobble] = positions;
                }
            }

            foreach (var (wobbledKMer, positions) in wobbledKMers)
            {
                allKMers[wobbledKMer] = positions;
            }
        }

        return allKMers;
    }

    /// <summary>
    /// Finds clumping k-mers in the <see cref="Sequence"/>.
    /// </summary>
    /// <param name="windowSize">The size of the window.</param>
    /// <param name="k">The length of the k-mer.</param>
    /// <param name="threshold">The minimum number of times the k-mer must appear in the window.</param>
    /// <returns>An array of clumping k-mers.</returns>
    public Sequence[] FindClumpingKMers(int windowSize, int k, int threshold)
    {
        var clumpingKMers = new List<Sequence>();
        foreach (var (kMer, positions) in FindAllKMers(k))
        {
            if (positions.Count < threshold)
            {
                continue;
            }

            for (var i = 0; i < positions.Count - threshold + 1; i++)
            {
                var start = positions[i];
                var end = positions[i + threshold - 1];

                if (end - start + 1 <= windowSize - k)
                {
                    clumpingKMers.Add(kMer);
                    break;
                }
            }
        }

        return [..clumpingKMers];
    }

    /// <summary>
    /// Use to find the skew of the sequence for all nucleotides.
    /// The skew is the difference between the count of <see cref="Nucleotide.G"/> and <see cref="Nucleotide.C"/>.
    /// </summary>
    /// <returns>An array of the skew of the sequence for all <see cref="Nucleotide"/>s.</returns>    
    public int[] GetSkew()
    {
        var skew = new int[Length + 1];
        for (var i = 1; i <= Length; i++)
        {
            skew[i] = skew[i - 1] + _sequence[i - 1] switch
            {
                Nucleotide.G => 1,
                Nucleotide.C => -1,
                _ => 0
            };
        }

        return skew;
    }

    /// <summary>
    /// Use to find the positions where the skew is minimum.
    /// The skew is the difference between the count of <see cref="Nucleotide.G"/> and <see cref="Nucleotide.C"/>.
    /// </summary>
    /// <returns>An array of positions where the skew is minimum.</returns>
    public int[] GetMinSkewPositions()
    {
        var skew = GetSkew();
        var min = skew.Min();

        return skew
            .Select((value, index) => new { value, index })
            .Where(x => x.value == min)
            .Select(x => x.index)
            .ToArray(); 
    }

    /// <summary>
    /// Use to find the Hamming Distance between two sequences.
    /// The Hamming Distance is the number of positions at which the corresponding nucleotides differ.
    /// </summary>
    /// <param name="other">The other <see cref="Sequence"/>.</param>
    /// <returns>The Hamming distance between the two sequences.</returns>
    /// <exception cref="ArgumentException">Thrown when the sequences have different lengths.</exception>
    public int HammingDistance(Sequence other)
    {
        if (Length != other.Length)
        {
            throw new ArgumentException("The sequences must have the same length.");
        }

        var distance = 0;
        for (var i = 0; i < Length; i++)
        {
            if (_sequence[i] != other[i])
            {
                distance++;
            }
        }

        return distance;
    }

    /// <summary>
    /// Use to find the pattern in the sequence. May consider wobbles.
    /// </summary>
    /// <param name="pattern">The pattern <see cref="Sequence"/> to find.</param>
    /// <param name="mismatches">The maximum number of mismatches (wobbles) allowed in the pattern.</param>
    /// <returns>An array of positions where the pattern is found.</returns>
    public int[] FindPatternLocations(Sequence pattern, int mismatches = 0)
    {
        var positions = new List<int>();
        for (var i = 0; i <= Length - pattern.Length; i++)
        {
            var subsequence = Part(i, pattern.Length);
            if (subsequence.HammingDistance(pattern) <= mismatches)
            {
                positions.Add(i);
            }
        }

        return [..positions];
    }
    
    /// <summary>
    /// Use to find the wobbles of the <see cref="Sequence"/> with at most <paramref name="mismatches"/> mismatches.
    /// </summary>
    /// <param name="mismatches">The maximum number of mismatches allowed.</param>
    /// <returns>An array of wobbles of the <see cref="Sequence"/>.</returns> 
    public Sequence[] Wobbles(int mismatches)
    {
        if (mismatches == 0)
        {
            return [this];
        }

        if (Length == 1)
        {
            return NucleotideHelper.All
                .Select(n => new Sequence([n]))
                .ToArray();
        }

        var suffix = Part(1, Length - 1);
        var suffixNeighbors = suffix.Wobbles(mismatches);

        var wobbles = new List<Sequence>();
        foreach (var wobble in suffixNeighbors)
        {
            if (suffix.HammingDistance(wobble) < mismatches)
            {
                foreach (var nucleotide in NucleotideHelper.All)
                {
                    wobbles.Add(nucleotide + wobble);
                }
            }
            else
            {
                wobbles.Add(this[0] + wobble);
            }
        }

        return [..wobbles];
    }

    /// <summary>
    /// Use to mutate the <see cref="Sequence"/>> at the specified index.
    /// </summary>
    /// <param name="index">The index to mutate.</param>
    /// <param name="nucleotide">The <see cref="Nucleotide"/> to mutate to.</param>
    /// <returns>Mutated <see cref="Sequence"/></returns>
    public Sequence Mutate(int index, Nucleotide nucleotide)
    {
        var mutated = new Nucleotide[Length];
        Array.Copy(_sequence, mutated, Length);
        mutated[index] = nucleotide;

        return new(mutated);
    }

    public override string ToString()
    {
        var sb = new StringBuilder();
        foreach (var nucleotide in _sequence)
        {
            sb.Append(nucleotide.ToChar());
        }

        return sb.ToString();
    }

    public bool Equals(Sequence? other)
    {
        if (other is null)
        {
            return false;
        }

        if (Length != other.Length)
        {
            return false;
        }

        for (var i = 0; i < Length; i++)
        {
            if (_sequence[i] != other[i])
            {
                return false;
            }
        }

        return true;
    }

    public override int GetHashCode()
    {
        var hash = 17;
        foreach (var nucleotide in _sequence)
        {
            hash = hash * 31 + nucleotide.GetHashCode();
        }

        return hash;
    }

    private static Nucleotide[] GetNucleotides(string sequence)
    {
        var nucleotides = new Nucleotide[sequence.Length];
        for (var i = 0; i < sequence.Length; i++)
        {
            nucleotides[i] = char.ToUpper(sequence[i]).ToNucleotide();
        }

        return nucleotides;
    }
}