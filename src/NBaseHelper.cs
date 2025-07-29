namespace GenLabs.DnaUtils;

public static class NBaseHelper
{
    /// <summary>
    /// All possible nucleotide bases.
    /// </summary>
    public static NBase[] All =
        [NBase.A, NBase.C, NBase.G, NBase.T];

    /// <summary>
    /// Generates a random sequence of nucleotide bases.
    /// </summary>
    /// <param name="length">The length of the sequence.</param>
    /// <returns>The random sequence.</returns> 
    public static NBase[] Random(int length)
    {
        var random = new Random();
        var nBases = new NBase[length];
        for (var i = 0; i < length; i++)
        {
            nBases[i] = All[random.Next(0, 4)];
        }

        return nBases;
    }
}
