namespace GenLabs.DnaUtils;

public static class NucleotideExtensions
{
    /// <summary>
    /// Converts a character to a <see cref="Nucleotide"/>.
    /// </summary>
    /// <param name="nucleotide">The character to convert.</param>
    /// <returns>The <see cref="Nucleotide"/>.</returns>
    /// <exception cref="ArgumentException">Thrown when the character is not a valid DNA character.</exception>
    public static Nucleotide ToNucleotide(this char nucleotide) =>
        nucleotide switch
        {
            'A' => Nucleotide.A,
            'C' => Nucleotide.C,
            'G' => Nucleotide.G,
            'T' => Nucleotide.T,
            _ => throw new ArgumentException("Invalid DNA character")
        };

    /// <summary>
    /// Converts a <see cref="Nucleotide"/> to a character.
    /// </summary>
    /// <param name="nucleotide">The <see cref="Nucleotide"/> to convert.</param>
    /// <returns>The character.</returns>
    /// <exception cref="ArgumentException">Thrown when the <see cref="Nucleotide"/> is not a valid DNA character.</exception>
    public static char ToChar(this Nucleotide nucleotide) =>
        nucleotide switch
        {
            Nucleotide.A => 'A',
            Nucleotide.C => 'C',
            Nucleotide.G => 'G',
            Nucleotide.T => 'T',
            _ => throw new ArgumentException("Invalid DNA character")
        };

    /// <summary>
    /// Returns the complement of a <see cref="Nucleotide"/>.
    /// The complemen is the nucleotide that pairs with the given nucleotide.
    /// <see cref="Nucleotide.A"/> pairs with <see cref="Nucleotide.T"/>, and <see cref="Nucleotide.C"/> pairs with <see cref="Nucleotide.G"/>.
    /// </summary>
    /// <param name="nucleotide">The <see cref="Nucleotide"/> to complement.</param>
    /// <returns>The complement.</returns>
    /// <exception cref="ArgumentException">Thrown when the <see cref="Nucleotide"/> is not a valid DNA character.</exception>
    public static Nucleotide Complement(this Nucleotide nucleotide) =>
        nucleotide switch
        {
            Nucleotide.A => Nucleotide.T,
            Nucleotide.T => Nucleotide.A,
            Nucleotide.C => Nucleotide.G,
            Nucleotide.G => Nucleotide.C,
            _ => throw new ArgumentException("Invalid DNA character")
        };
}