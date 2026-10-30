namespace GenLabs.DnaUtils;

/// <summary>
/// Helpers for working with the set of <see cref="NBase"/> values.
/// </summary>
public static class NBaseHelper
{
    /// <summary>
    /// All possible nucleotide bases, in alphabetical order.
    /// </summary>
    public static IReadOnlyList<NBase> All { get; } = [NBase.A, NBase.C, NBase.G, NBase.T];

    /// <summary>
    /// Generates a random sequence of nucleotide bases.
    /// </summary>
    /// <param name="length">The length of the sequence.</param>
    /// <param name="random">The source of randomness; <see cref="System.Random.Shared"/> when omitted. Pass a seeded instance for reproducible output.</param>
    /// <returns>The random sequence.</returns>
    public static NBase[] Random(int length, Random? random = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        random ??= System.Random.Shared;

        var nBases = new NBase[length];
        for (var i = 0; i < length; i++)
        {
            nBases[i] = All[random.Next(All.Count)];
        }

        return nBases;
    }
}
