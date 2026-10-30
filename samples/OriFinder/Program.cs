// Looks for a bacterial genome's origin of replication (ori). The GC skew reaches its minimum
// where replication switches strands, near the ori. There, the replication-initiator protein
// DnaA binds short "DnaA boxes": 9-mers that recur, approximately and on either strand.
// Run it on a genome from GenLabs.Genome, optionally with a known DnaA box to check:
//
//   dotnet run --project samples/OriFinder -- path/to/escherichia_coli.txt TTATCCACA

using GenLabs.DnaUtils;

const int WindowSize = 500;
const int K = 9;
const int Mismatches = 1;

if (args.Length is < 1 or > 2)
{
    Console.Error.WriteLine("Usage: OriFinder <genome.txt> [known 9-mer DnaA box]");
    return 1;
}

var genome = new Sequence(File.ReadAllText(args[0]).Trim());
Console.WriteLine($"Genome: {genome.Length:N0} bases");

var minSkew = genome.GetMinSkewPositions();
Console.WriteLine($"Minimum skew at: {string.Join(", ", minSkew)}");

var start = Math.Min(minSkew[0], genome.Length - WindowSize);
var window = genome.Part(start, WindowSize);
var kMers = window.FindAllKMers(K, Mismatches, treatReverseComplementsAsEqual: true);
var max = kMers.Values.Max(p => p.Count);
var mostFrequent = kMers.Where(k => k.Value.Count == max).Select(k => k.Key).ToArray();

Console.WriteLine(
    $"In the {WindowSize} bases from {start}, {mostFrequent.Length} {K}-mers tie as most frequent " +
    $"({max} occurrences, up to {Mismatches} mismatch, either strand).");

if (args.Length == 2)
{
    Sequence box = args[1];
    var occurrences = kMers.TryGetValue(box, out var positions) ? positions.Count : 0;
    var rank = kMers.Values.Count(p => p.Count > occurrences) + 1;

    Console.WriteLine($"{box} occurs {occurrences} times (reverse complement {box.ReverseComplement()}): rank {rank} by frequency.");
    Console.WriteLine($"Approximate matches on the forward strand: {string.Join(", ", window.FindPatternLocations(box, Mismatches).Select(p => start + p))}");
    Console.WriteLine($"Approximate matches on the reverse strand: {string.Join(", ", window.FindPatternLocations(box.ReverseComplement(), Mismatches).Select(p => start + p))}");
}

return 0;
