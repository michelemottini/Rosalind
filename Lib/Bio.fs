// Bioinformatics functions

module Bio

open Combi

/// Returns the index (0 to 3) of a DNA nucleotide
let private nucleotideIndex (c: char) =
  "ACGT".IndexOf(c)

/// Returns the number of each nucleotides contained in a dna string
/// The result is an array with four elements containing the counters for 'A', 'C', 'G' and 'T'
let countingDNANucleotides (dna: string) =
  let counters = [|0; 0; 0; 0|];
  String.iter (fun c -> 
    let index = nucleotideIndex c
    counters.[index] <- counters.[index] + 1
  ) dna
  counters

/// Parse data in FASTA format from a text reader, returning a sequence of (label, string) pairs
let parseFastaReader (reader: System.IO.TextReader) =
  let rec skip () =
    let line = reader.ReadLine()
    if line <> null then
      if line.StartsWith(">") then first line
      else skip ()
    else
      Seq.empty
  and first line =
    readString (line.Substring(1)) System.String.Empty
  and readString label current =
    let line = reader.ReadLine()
    if line = null then 
      Seq.singleton (label, current)
    elif line.StartsWith(">") then
      Seq.append (Seq.singleton (label, current)) (first line)
    else
      readString label (current+line)
  skip ()

/// Parse data in FASTA format from a string, returning a sequence of (label, string) pairs
let parseFastaString (str: string) =
  use reader = new System.IO.StringReader(str)
  parseFastaReader reader

/// Convert a DNA nucleotid into the corresponding RNA one
let private dna2rna c = 
  if c = 'T' then 'U' else c

/// Convert a DNA string into the corresponding RNA one
let transcribingDNAintoRNA (dna: string) =
  String.map dna2rna dna

/// Compute the percentage of 'G' and 'C' nucleotids in the specified DNA string
let computeGCContent (dna: string) =
  let gcCount = 
    String.count (fun c -> c = 'G' || c = 'C') dna
  (float gcCount)*100.0 / (float dna.Length)

// http://en.wikipedia.org/wiki/Genetic_code#RNA_codon_table
// empty string = stop
let private codonTableList =
  [
    ("UUU", "F"); ("CUU", "L"); ("AUU", "I"); ("GUU", "V")
    ("UUC", "F"); ("CUC", "L"); ("AUC", "I"); ("GUC", "V")
    ("UUA", "L"); ("CUA", "L"); ("AUA", "I"); ("GUA", "V")
    ("UUG", "L"); ("CUG", "L"); ("AUG", "M"); ("GUG", "V")
    ("UCU", "S"); ("CCU", "P"); ("ACU", "T"); ("GCU", "A")
    ("UCC", "S"); ("CCC", "P"); ("ACC", "T"); ("GCC", "A")
    ("UCA", "S"); ("CCA", "P"); ("ACA", "T"); ("GCA", "A")
    ("UCG", "S"); ("CCG", "P"); ("ACG", "T"); ("GCG", "A")
    ("UAU", "Y"); ("CAU", "H"); ("AAU", "N"); ("GAU", "D")
    ("UAC", "Y"); ("CAC", "H"); ("AAC", "N"); ("GAC", "D")
    ("UAA", "");  ("CAA", "Q"); ("AAA", "K"); ("GAA", "E")
    ("UAG", "");  ("CAG", "Q"); ("AAG", "K"); ("GAG", "E")
    ("UGU", "C"); ("CGU", "R"); ("AGU", "S"); ("GGU", "G")
    ("UGC", "C"); ("CGC", "R"); ("AGC", "S"); ("GGC", "G")
    ("UGA", "");  ("CGA", "R"); ("AGA", "R"); ("GGA", "G")
    ("UGG", "W"); ("CGG", "R"); ("AGG", "R"); ("GGG", "G") 
  ]

/// Map from all the possible codons to the corresponding aminoacid. Stop codons are mapped to empty strings
let codonTableMap = 
  Map(codonTableList)

/// Map from proteins to the corresponding codons
let inverseCodonTableMap = 
  codonTableMap |> Map.fold (fun inverseMap codon aminoacid -> 
    match Map.tryFind aminoacid inverseMap with
      | Some(codons) -> inverseMap |> Map.add aminoacid (codon::codons)
      | None -> inverseMap |> Map.add aminoacid [codon]
  ) Map.empty

/// Convert a RNA string into the corresponding protein using the codon map
let rna2protein (rna: string) = 
  seq {
    for i in 0..3..(rna.Length-1) do
      yield Map.find (rna.Substring(i, 3)) codonTableMap
  } |> String.concat ""

/// Complement of a DNA nucelotide
let nucleotideComplement (c: char) =
  match c with
    | 'A' -> 'T'
    | 'C' -> 'G'
    | 'G' -> 'C'
    | 'T' -> 'A'
    | _ -> invalidArg "c" (sprintf "Unexpected value '%c'" c)

/// Complement of a DNA string
let dnaComplement (dna: string) =
  dna |> String.reverse |> String.map nucleotideComplement

/// Computes the Hamming distance between two DNA strings
let hammingDistance (dna1: string) (dna2: string) =
  let result = ref 0
  dna1 |> String.counti (fun i c -> dna2.[i] <> c)

let private monoisotopicMass = 
  [
    ('A', 71.03711)
    ('C', 103.00919) 
    ('D', 115.02694) 
    ('E', 129.04259) 
    ('F', 147.06841) 
    ('G', 57.02146) 
    ('H', 137.05891) 
    ('I', 113.08406) 
    ('K', 128.09496) 
    ('L', 113.08406) 
    ('M', 131.04049) 
    ('N', 114.04293) 
    ('P', 97.05276) 
    ('Q', 128.05858) 
    ('R', 156.10111) 
    ('S', 87.03203) 
    ('T', 101.04768)
    ('V', 99.06841) 
    ('W', 186.07931)
    ('Y', 163.06333)
  ]

/// Map of individual aminoacids monoisotopic masses
let monoisotopicMassMap = Map(monoisotopicMass)

/// Computes the monoisotopic mass of a protein
let proteinMass (protein: string) =
  protein |> Seq.fold (fun mass c -> mass + (Map.find c monoisotopicMassMap)) 0.0

/// Converts a protein motif string into the corresponding regular expression
/// Warning: it does not check the validity of the input motif string
let rec private proteinMotifRegex (motif: string) =
  if motif.Length = 0 then 
    ""
  else
    match motif.[0] with
      | '[' -> 
        let endIndex = motif.IndexOf(']')
        (motif.Substring(0,endIndex+1)) + (proteinMotifRegex (motif.Substring(endIndex+1)))
      | '{' ->
        let endIndex = motif.IndexOf('}')
        "[^" + (motif.Substring(1,endIndex-1)) + "]" + (proteinMotifRegex (motif.Substring(endIndex+1)))
      | _ -> (motif.Substring(0,1)) + (proteinMotifRegex (motif.Substring(1)))

/// Find all the locations of the specified protein motif inside a protein string.
/// Returns a list of 0-based positions - or an empty list if the motif never appears in the protein
let proteinMotifLocations (motif: string) (protein: string) = 
  let motifRegex = System.Text.RegularExpressions.Regex(proteinMotifRegex motif)
  let rec findLocations offset str = 
    let m = motifRegex.Match(str)
    if not(m.Success) then []
    else
      let index1 = m.Index+1
      (offset+m.Index)::(findLocations (offset+index1) (str.Substring(index1)))
  findLocations 0 protein

/// Fetch a protein string given its id from the Uniprot online database
let fetchProtein (proteinId: string) = 
  let proteinFASTAUrl = "http://www.uniprot.org/uniprot/" + System.Web.HttpUtility.UrlEncode(proteinId) + ".fasta"
  try
    let request = System.Net.WebRequest.Create(proteinFASTAUrl) :?> System.Net.HttpWebRequest
    request.UserAgent <- "Mozilla/5.0 (Windows NT 6.1; WOW64)" // For some reason the UNIPROT service requires a user agent to work...
    use response = request.GetResponse()
    use responseStream = response.GetResponseStream()
    let (id,protein) = parseFastaReader (new System.IO.StreamReader(responseStream)) |> Seq.cache |> Seq.head
    protein
  with
    // Better error messages in case of Web errors
    | :? System.Net.WebException as webException ->
      match webException.Response with
        | :? System.Net.HttpWebResponse as errorResponse ->
          use errorResponseStream = errorResponse.GetResponseStream()
          let errorContent = (new System.IO.StreamReader(errorResponseStream)).ReadToEnd()
          raise (System.ApplicationException(sprintf "Unable to access '%s': %s\r\n%s" proteinFASTAUrl errorResponse.StatusDescription errorContent))
        | _ ->
          raise (System.ApplicationException(sprintf "Unable to access '%s'" proteinFASTAUrl))


/// Fetch protein strings given their ids from the Uniprot online database
let fetchProteins (ids: string) =
  ids.Split('\r', '\n') |> Seq.iter (fun id ->
    let protein = fetchProtein id
    printfn "%s" id
    printfn "%s" protein
  )

/// Computes the profile matrix of a sequence of DNA strings with the same length
let profileMatrix (dnas: string seq) = 
  if Seq.isEmpty dnas then
    Array2D.create 4 0 0
  else
    let n = dnas |> Seq.head |> String.length
    let result = Array2D.create 4 n 0
    dnas |> Seq.iter (fun dna -> 
      dna |> String.iteri (fun j c ->
        let i = nucleotideIndex c
        result.[i, j] <- result.[i, j] + 1
      )
    )
    result

/// Compute the consensus string and the profile matrix of a sequence of DNA strings with the same length
let consensusAndProfile (dnas: string seq) = 
  let profile = dnas |> profileMatrix
  let n = Array2D.length2 profile
  let max = Array.create n (-1,0)
  profile |> Array2D.iteri (fun i j freq ->
    let (maxi,maxfreq) = max.[j]
    if (freq > maxfreq) then
      max.[j] <- (i,freq)
    else
      ()
  )
  let maxChars =
    max |> Array.map (fun (maxi, _) -> "ACGT".[maxi]) 
  (System.String(maxChars)), profile

/// Compute the sequence of all the possible codons in a RNA string starting from the specified offset
let private codons offset (rna: string) = 
  seq {
    for i in offset..3..((rna.Length-offset)/3*3-1) do
      yield (rna.Substring(i, 3))
  }

/// Process a codon adding it to the set of pending partial open reading frames, mvoing them 
/// them to the set of open reading fraames when the codon is a stop one
let private processCodon (orfs, partialOrfs) codon =
  let appendAminoacid aminoacid (proteins: string list) = 
    proteins |> List.map (fun protein -> protein + aminoacid)
  let aminoacid = Map.find codon codonTableMap
  if aminoacid = "M" then
    // Start codon: it is the beginning of a new possible (partial) ORF, and should also be added to all the current partial ORFs
    (orfs, aminoacid::(partialOrfs |> appendAminoacid aminoacid))
  elif aminoacid = "" then
    // Stop codon: all the current partial ORFs become actual ORFs, the set of partial ORFs is empty
    (Set.union orfs (Set.ofList partialOrfs), [])
  else
    // Any other codon: append it to all partial ORFs
    (orfs, partialOrfs |> appendAminoacid  aminoacid)

/// Computes the set of all the possible open reading frames of the specifeid RNA string
let rnaOpenReadingFrames (rna: string) = 
  [0..2] |> Seq.map (fun offset -> 
    let (orfs, _) = codons offset rna |> Seq.fold processCodon (Set.empty,[])
    orfs
  ) |> Set.unionMany

/// All the possible proteins that can be translated from open reading frames of the specified DNA string
let dnaOpenReadingFrames (dna: string) = 
  let toOpenReadingFrames = transcribingDNAintoRNA >> rnaOpenReadingFrames
  Set.union (dna |> toOpenReadingFrames) (dna |> dnaComplement |> toOpenReadingFrames)

/// Delete all the specified sub-strings (introns) from a DNA string, and returns the protein corresponding to the remaining DNA string
let splice (introns: string seq) (dna: string) = 
  let intronsMaxLength (str: string) = 
    introns |> Seq.map (fun intron -> if str.StartsWith(intron) then intron.Length else 0) |> Seq.max
  let rec proc (str: string) = 
    if String.isEmpty str then
      String.empty
    else
      let intronLen = intronsMaxLength str
      if intronLen > 0 then 
        proc (str.Substring(intronLen))
      else
        (dna2rna str.[0] |> String.ofChar) + (str.Substring(1) |> proc)
  proc dna |> rna2protein

/// Find the longest common substring of a sequence of DNA strings (if any)
/// Returns None if there is no common substring
let sharedMotif (dnas: string seq) = 
  let isSharedSubstring (str: string) =
    dnas |> Seq.forall (fun dna -> dna.IndexOf(str) >= 0)
  let minLengthDna =
    dnas |> Seq.minBy (fun dna -> dna.Length)
  let minLength = 
    minLengthDna.Length
  let minLengthDnaSubstrings =
    seq {
      for length in minLength..(-1)..1 do
        for i in 0..minLength-length do
          yield minLengthDna.Substring(i,length) 

    }
  minLengthDnaSubstrings |> Seq.tryFind isSharedSubstring

/// Checks if the speficied section of a DNA string is a reverse palindrome - i.e. it is equal to its reverse complement
let rec private isReversePalindrome (dna: string) start length =
  if length = 0 then
    true
  elif length = 1 then
    false
  else
    dna.[start] = nucleotideComplement (dna.[start+length-1]) && (isReversePalindrome dna (start+1) (length-2))

/// Returns the position and length of all possible restriction sites - i.e. DNA fragments that are reverse 
/// palindromes - with length between the specifies minimum and maximum
let restrictionSites (dna: string) minLength maxLength =
  seq {
    for start in 0..dna.Length-minLength do
      let curMaxLength = min maxLength (dna.Length-start)
      for length in minLength..curMaxLength do
        if isReversePalindrome dna start length then
          yield (start+1, length)
  }

/// Computes the total possible number of perfect matchings of basepair edges in the bonding graph of the specified RNA string 
/// The RNA string must contain the same number of occurences of 'A' and 'U' and of 'C' and 'G'
(*
  An 'U' node can be combined with any of the 'A' nodes - for a total of Nau (= number of AU pairs) possible pairs, 
  so the number of perfect matchings is Nau times the number of perfect matching in the remaining string - 
  that has one less AU pair. 
  Continuing recursively until all AU pairs are exausted gives Na! times the number of perfect matching in
  the remaining string - that has only GC pairs. Using the same reasoning the number of perfect matchings in 
  this string is Ngc! - hence the final result is Nau! * Ngc!
*)
let perfectMatchingCount (rna: string) =
  let auCount = 
    rna |> String.count (fun c -> c = 'A')
  let gcCount = 
    rna.Length/2 - auCount
  (factorial auCount) * (factorial gcCount)

/// Computes the probability of the specified dna string given the expected GC-content
let probFromGC (dna: string) gc = 
  let probGC = gc / 2.0
  let probAT = (1.0 - gc) / 2.0
  let probs = [| probAT; probGC; probGC; probAT |]
  dna |> Seq.fold (fun r c -> r*probs.[nucleotideIndex(c)]) 1.0

