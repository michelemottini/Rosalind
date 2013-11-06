// See http://rosalind.info

(*
  This script contain one function for each of the Rosalind problems solved so far. 

  Each of these functions is named as the corresponding problem and accept a single string as input and 
  produces a single string as output - both formatted exactly as in the Rosalind test cases.

  Each of these functions is called once by the general 'test' function, that passes a test case string in input
  and compares the output to the expected result, displaying 'OK' or an error message if they don't match. 
  
  The output of running the entire script should be the list of test names with 'OK' next to them:

    dna: OK
    rna: OK
    revc: OK
    iprb: OK
    fib: OK
    . . . 

  The actual work for each problems is done by functions defined in the separate modules String.fs (general string
  processing), Combi.fs (combinatorial functions), Bio.fs (strictly bioinformatics stuff). The idea is those modules
  constitute general libraries that can be used for problems different - and more general - that just the Rosalind ones.
*)


#load "Lib/String.fs"
#load "Lib/Combi.fs"
#load "Lib/Bio.fs"
#load "Lib/Utils.fs"

open Combi
open Bio
open Utils

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/dna/

let dna str = 
  let counters = countingDNANucleotides str
  counters |> integersToString

test "dna" dna "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC" "20 12 17 21"

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/rna/

let rna str = 
  transcribingDNAintoRNA str

test "rna" rna "GATGGAACTTGACTACGTAAATT" "GAUGGAACUUGACUACGUAAAUU" 

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/revc/

let revc str =
  dnaComplement str

test "revc" revc "AAAACCCGGT" "ACCGGGTTTT"

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/iprb/

let dominantProbability k m n =
  let total = k + m + n
  let mul a b =
    float (a*b)
  ((mul k (k-1)) + (mul k m) + (mul m k) + 0.75*(mul m (m-1)) + 0.5*(mul m n) + 0.5*(mul n m) + (mul k n) + (mul n k)) / (mul total (total-1))

let iprb str = 
  let kmn = toIntegers str 
  sprintf "%1.5f" (dominantProbability kmn.[0] kmn.[1] kmn.[2])

test "iprb" iprb "2 2 2" "0.78333" 

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/fib/

let fib str =
  let nk = toIntegers str
  sprintf "%d" (generalizedFibonacci nk.[0] nk.[1])

test "fib" fib "5 3" "19"
    
//-----------------------------------------------------------------------------------------------------------------
// http://rosalind.info/problems/gc/

let findMaxGCContent dnas = 
  dnas 
    |> Seq.maxBy (fun (label, dna) -> computeGCContent dna) 
    |> fun (label, dna) -> (label, computeGCContent dna)

let gc str = 
  let (id, gcContent) = str |> parseFastaString |> findMaxGCContent
  sprintf "%s\n%1.6f" id gcContent

test "gc" gc ">Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT" "Rosalind_0808
60.919540"

//-----------------------------------------------------------------------------------------------------------------
// http://rosalind.info/problems/prot/

let prot str = 
  rna2protein str

test "prot" prot "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA" "MAMAPRTEINSTRING" 

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/subs/

let subs str = 
  use reader = new System.IO.StringReader(str)
  let dna = reader.ReadLine()
  let motif = reader.ReadLine()
  dna |> String.multipleIndexOf motif |> Seq.map (fun idx -> idx+1) |> integersToString

test "subs" subs "GATATATGCATATACTT
ATAT" "2 4 10"

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/hamm/

let hamm str =
  use reader = new System.IO.StringReader(str)
  let dna1 = reader.ReadLine()
  let dna2 = reader.ReadLine()
  sprintf "%d" (hammingDistance dna1 dna2)

test "hamm" hamm "GAGCCTACTAACGGGAT
CATCGTAATGACGGCCT" "7"

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/iev/

let dominantExpected p1 p2 p3 p4 p5 p6 =
  (float (p1+p2+p3))*2.0 + (float p4)*1.5 + (float p5)

let iev str =
  let ps = toIntegers str
  sprintf "%1.1f" (dominantExpected ps.[0] ps.[1] ps.[2] ps.[3] ps.[4] ps.[5])

test "iev" iev "1 0 0 1 0 1" "3.5"

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/fibd/

let fibd str =
  let nm = toIntegers str
  sprintf "%d" (mortalFibonacci nm.[0] nm.[1])

test "fibd" fibd "6 3" "4"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/mrna/

let countPossibleRNAs (protein: string) = 
  let codonsCount aminoacid = 
    List.length (Map.find aminoacid inverseCodonTableMap)
  protein |> Seq.fold (fun count aminoacid -> 
    count*(codonsCount (System.String(aminoacid, 1))) % 1000000
  ) (codonsCount "")

let mrna str = 
  sprintf "%d" (countPossibleRNAs str)

test "mrna" mrna "MA" "12"

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/lia/

/// Probability of a transition from factor AA, Aa, aa to those same factors when mating with Aa
let baseTransitionProbability = 
  array2D [
            [ 0.5;  0.5; 0.0 ]
            [ 0.25; 0.5; 0.25]
            [ 0.0;  0.5; 0.5 ] 
          ]

/// Probability of a transition from factor AA-BB, AA-Bb, AA-bb, Aa-BB, Aa-Bb, Aa-bb, aa-BB, aa-Bb, aa-bb to those same factors when mating with Aa-Bb
let transitionProbability = 
  Array2D.init 9 9 (fun i j -> baseTransitionProbability.[i / 3, j / 3]*baseTransitionProbability.[i % 3, j % 3])

let rec allelesProbability k =
  if k = 0 then 
    [| 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; |]
  else
    let previousProbability = allelesProbability (k-1)
    let result = [| 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; |]
    for i in 0..8 do
      for j in 0..8 do
        result.[i] <- result.[i] + transitionProbability.[j, i]*previousProbability.[j]
    result

let alleleProbabilityGreaterOrEqual k n = 
  let p = (allelesProbability k).[4]
  binomialProbabilityGreaterOrEqual p (1 <<< k) n

let lia str =
  let kn = toIntegers str
  sprintf "%1.3f" (alleleProbabilityGreaterOrEqual kn.[0] kn.[1])

test "lia" lia "2 1" "0.684"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/prtm/

let prtm str =
  sprintf "%1.3f" (proteinMass str)

test "prtm" prtm "SKADYEK" "821.392"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/grph/

let areConnected k (dna1: string) (dna2: string) = 
  dna1.Length >= k && dna2.Length >= k && dna1.Substring(dna1.Length-k,k) = dna2.Substring(0,k)

let overlapGraph k (idDnas: (string*string) seq) = 
  idDnas |> Seq.collect (fun (id1, dna1) ->
    idDnas |> Seq.collect (fun (id2, dna2) ->
      if id1 <> id2 && (areConnected k dna1 dna2) then
        Seq.singleton (id1, id2)
      else
        Seq.empty
    )
  )

let grph str =
  let idDnas = parseFastaString str
  let graph = overlapGraph 3 idDnas
  graph |> Seq.map (fun (id1, id2) -> sprintf "%s %s" id1 id2) |> String.concat "\n"

test "grph" grph ">Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG" "Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/mprt/

let mprt (motif: string) (ids: string) =
  seq {
    for id in ids.Split('\r', '\n') do
      let protein = fetchProtein id
      let motifLocations = proteinMotifLocations motif protein
      if not(List.isEmpty motifLocations) then
         yield id + "\n" + (motifLocations |> List.map (fun idx -> idx+1) |> integersToString)
  } |> String.concat "\n"

test "mprt" (mprt "N{P}[ST]{P}") "A2Z669
B5ZC00
P07204_TRBM_HUMAN
P20840_SAG1_YEAST" "B5ZC00
85 118 142 306 395
P07204_TRBM_HUMAN
47 115 116 382 409
P20840_SAG1_YEAST
79 109 135 248 306 348 364 402 485 501 614"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/cons/

let cons str = 
  let dnas = parseFastaString str |> Seq.map (fun (_,dna) -> dna)
  let maxChars, profile = consensusAndProfile dnas
  let result = System.Text.StringBuilder()
  do result.Append(maxChars) |> ignore
  for i in 0..3 do
    do result.Append("\n").Append("ACGT".[i]).Append(":") |> ignore
    for j in 0..profile.GetLength(1)-1 do
      do result.Append(" ").Append(profile.[i,j]) |> ignore
  result.ToString()

test "cons" cons ">Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT" "ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/orf/

let orf (str: string) = 
  let (_, dna) = parseFastaString str |> Seq.head
  let orfs = dnaOpenReadingFrames dna
  orfs |> String.concat "\n"

test "orf" orf ">Rosalind_99
AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG" "M
MGMTPRLGLESLLE
MLLGSFRLIPKETLIQVAGSSPCNLS
MTPRLGLESLLE"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/splc/

let splc str = 
  let dnas = parseFastaString str |> Seq.map snd |> List.ofSeq
  match dnas with
    | dna::introns -> splice introns dna
    | _ -> invalidArg "str" "Must contain at least two DNA strings"

test "splc" splc ">Rosalind_10
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
>Rosalind_12
ATCGGTCGAA
>Rosalind_15
ATCGGTCGAGCGTGT" "MVYIADKQHVASREAYGHMFKVCA"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/lcsm/

let lcsm str = 
  let dnas = parseFastaString str |> Seq.map (fun (_, dna) -> dna)
  match sharedMotif dnas with
    | None -> String.empty
    | Some(str) -> str

test "lcsm" lcsm ">Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA" "TA"

//-----------------------------------------------------------------------------
// http://rosalind.info/problems/perm/

let perm str =
  let perms = permutations (System.Int32.Parse(str))
  let lengthString = (perms |> Seq.length).ToString()
  let permsString = perms |> Seq.map (fun perm -> integersToString perm) |> String.concat "\n"
  lengthString + "\n" + permsString

test "perm" perm "3" "6
1 2 3
1 3 2
2 1 3
2 3 1
3 1 2
3 2 1"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/revp/

let revp str =
  let (_, dna) = parseFastaString str |> Seq.head
  restrictionSites (dna: string) 4 12 |> Seq.map (fun (pos, len) ->
    sprintf "%d %d" pos len
  ) |> String.concat "\n"

test "revp" revp ">Rosalind_24
TCAATGCATGCGGGTCTATATGCAT" "4 6
5 4
6 6
7 4
17 4
18 4
20 6
21 4"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/pmch/

let pmch str = 
  let (_, rna) = parseFastaString str |> Seq.head
  let result = perfectMatchingCount rna
  result.ToString()

test "pmch" pmch ">Rosalind_23
AGCUAGUCAU" "12"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/pper/

let pper str =
  let nk = toIntegers str
  sprintf "%d" (partialPermutationsCountModulo nk.[0] nk.[1] 1000000)

test "pper" pper "21 7" "51200"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/tree/

/// Given an undirected graph described by the specified adjacency list returns all the nodes connected to the specified node
let nodesFrom adj node =
  seq {
    for (a,b) in adj do
      if a = node then 
        yield b
      elif b = node then
        yield a
  }

/// Given an undirected graph described by the specified adjacency list computes the set of all nodes connected to the specified node
/// optionally excluding the specified parentNode 
let rec connectedNodes adj parentNodeOpt node =
  nodesFrom adj node 
    |> Seq.map (fun curNode -> 
      match parentNodeOpt with 
        | Some(parentNode) ->
          if parentNode = curNode then
            Set.empty
          else
            connectedNodes adj (Some node) curNode
        | None -> connectedNodes adj (Some node) curNode) 
    |> Set.unionMany
    |> Set.add node

/// Given an undirected graph of n nodes described by the specified adjacency list computes the number of connected sub-graphs
let countConnected n adj = 
  [1..n] |> Seq.fold (fun (traversed,count) node ->
    if traversed |> Set.contains node then
      (traversed, count)
    else 
      let connectedToNode = connectedNodes adj None node
      (Set.union traversed connectedToNode, count+1)
  ) (Set.empty,0) |> snd

let tree str =
  use reader = new System.IO.StringReader(str)
  let n = System.Convert.ToInt32(reader.ReadLine())
  let adj = readLines reader |> Seq.map (fun line -> 
    let nodes = toIntegers line
    (nodes.[0], nodes.[1])
  )
  // The minimum number of edges necessary to completely connect the graph is the number of already connected sub-graphs minus 1 
  // i.e. the number of edges necessary to connect those already connected sub-graphs together
  ((countConnected n adj) - 1 ).ToString()

test "tree" tree "10
1 2
2 8
4 10
5 9
6 10
7 9" "3"

//------------------------------------------------------------------------------------
// http://rosalind.info/problems/prob/

let prob str =
  use reader = new System.IO.StringReader(str)
  let dna = reader.ReadLine()
  let a = reader.ReadLine() |> toFloats
  let b = a |> Array.map (probFromGC dna >> log10)
  b |> floatsToString 3

test "prob" prob "ACGATACAA
0.129 0.287 0.423 0.476 0.641 0.742 0.783" "-5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009"

