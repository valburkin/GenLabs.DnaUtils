namespace GenLabs.DnaUtils;

internal sealed class SequenceEqualityComparer : IEqualityComparer<Sequence>
{
    public bool Equals(Sequence? x, Sequence? y) => x?.Equals(y) ?? false;

    public int GetHashCode(Sequence obj) => obj.GetHashCode();
}