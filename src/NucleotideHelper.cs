namespace Genome.Biology;

public static class NucleotideHelper
{
    /// <summary>
    /// All possible nucleotides.
    /// </summary>
    public static Nucleotide[] All =
        [Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T];

    /// <summary>
    /// Generates a random sequence of nucleotides.
    /// </summary>
    /// <param name="length">The length of the sequence.</param>
    /// <returns>The random sequence.</returns> 
    public static Nucleotide[] Random(int length)
    {
        var random = new Random();
        var nucleotides = new Nucleotide[length];
        for (var i = 0; i < length; i++)
        {
            nucleotides[i] = All[random.Next(0, 4)];
        }

        return nucleotides;
    }
}