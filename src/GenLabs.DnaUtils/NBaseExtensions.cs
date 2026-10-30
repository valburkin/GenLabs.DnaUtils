namespace GenLabs.DnaUtils;

/// <summary>
/// Conversions and pairing rules for <see cref="NBase"/>.
/// </summary>
public static class NBaseExtensions
{
    /// <summary>
    /// Converts a character to a <see cref="NBase"/>. Lowercase characters are accepted.
    /// </summary>
    /// <param name="nBase">The character to convert.</param>
    /// <returns>The <see cref="NBase"/>.</returns>
    /// <exception cref="ArgumentException">Thrown when the character is not a valid DNA character.</exception>
    public static NBase ToNBase(this char nBase) =>
        char.ToUpperInvariant(nBase) switch
        {
            'A' => NBase.A,
            'C' => NBase.C,
            'G' => NBase.G,
            'T' => NBase.T,
            _ => throw new ArgumentException($"'{nBase}' is not a valid DNA character.", nameof(nBase)),
        };

    /// <summary>
    /// Converts a <see cref="NBase"/> to a character.
    /// </summary>
    /// <param name="nBase">The <see cref="NBase"/> to convert.</param>
    /// <returns>The character.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when the value is not a defined <see cref="NBase"/>.</exception>
    public static char ToChar(this NBase nBase) =>
        nBase switch
        {
            NBase.A => 'A',
            NBase.C => 'C',
            NBase.G => 'G',
            NBase.T => 'T',
            _ => throw new ArgumentOutOfRangeException(nameof(nBase)),
        };

    /// <summary>
    /// Returns the complement of a <see cref="NBase"/>.
    /// The complement is the nucleotide base that pairs with the given nucleotide base.
    /// <see cref="NBase.A"/> pairs with <see cref="NBase.T"/>, and <see cref="NBase.C"/> pairs with <see cref="NBase.G"/>.
    /// </summary>
    /// <param name="nBase">The <see cref="NBase"/> to complement.</param>
    /// <returns>The complement.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when the value is not a defined <see cref="NBase"/>.</exception>
    public static NBase Complement(this NBase nBase) =>
        nBase switch
        {
            NBase.A => NBase.T,
            NBase.T => NBase.A,
            NBase.C => NBase.G,
            NBase.G => NBase.C,
            _ => throw new ArgumentOutOfRangeException(nameof(nBase)),
        };
}
