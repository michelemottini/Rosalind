// Combinations functions

module Combi

open System.Numerics

/// Computes the factorial of an integer
let factorial n =
  if n < 0 then
    invalidArg "n" "Cannot be negative"
  else
    let rec loop i acc =
      match i with
      | 0 | 1 -> acc
      | _ -> loop (i-1) (acc * (BigInteger i))
    loop n 1I

/// Computes the binomial coefficient of the two specified integers
let rec binomialCoefficient n k = 
  if k = 0 then 1.0
  elif k = 1 then (float n)
  else (float n) * (binomialCoefficient (n-1) (k-1)) / (float k)

/// Computes the the probability of getting k outcomes out of n attempts for an event with probability p
let binomialProbability p n k = 
  (binomialCoefficient n k) * (p ** (float k)) * ((1.0-p) ** (float (n-k)))

/// Computes the the probability of getting k or more outcomes out of n attempts for an event with probability p
let binomialProbabilityGreaterOrEqual p n k = 
  seq {
    for i in k..n do
      yield binomialProbability p n i
  } |> Seq.sum

/// Computes the n-th value of a 'generalized' Fibonacci sequence, where each generation 
/// is the sum of the previous one plus k times the previous-previous one (k=1 gives the 
/// standard Fibonacci sequence)
let generalizedFibonacci n k =
  let rec compute n = 
    if n <= 2 then 1L
    else (compute (n-1)) + (compute (n-2))*(int64 k)
  compute n

/// Computes the total number of pairs of the n-th generation of a population where each individual lives m
/// generations, reaches maturity after one generation and each pair produces a single pair of offsprings at 
/// each generation
/// See http://rosalind.info/problems/fibd/
let mortalFibonacci n m =
  let memBirth = ResizeArray<int64>()
  let memDeath = ResizeArray<int64>()
  let memAlive = ResizeArray<int64>()
  let memoized (mem: ResizeArray<int64>) f n = 
    if n <= mem.Count then mem.Item(n-1)
    else
      let result = f n
      do
        if n = mem.Count+1 then mem.Add(result)
        else ()
      result
  let rec birth n = 
    memoized memBirth (fun n ->
      if n = 1 then 0L
      elif n = 2 then 1L
      else (alive (n-1)) - (death (n-1))
    ) n
  and death n =
    memoized memDeath (fun n -> 
      if n < m then 0L
      elif n = m then 1L
      else birth (n-m)
    ) n
  and alive n =
    memoized memAlive (fun n -> 
      if n=1 then 1L
      else (alive (n-1)) + (birth (n-1)) - (death (n-1))
    ) n
  alive n

let rec private permutationsOf (symbols:Set<int>) =
  if Set.count symbols = 1 then 
    symbols |> List.ofSeq |> Seq.singleton
  else 
    symbols |> Seq.collect (fun n -> 
      let reducedSymbols = symbols |> Set.remove n 
      let subPermutations = permutationsOf reducedSymbols
      subPermutations |> Seq.map (fun perm -> n::perm)
    )

/// Computes all the permutations of the first n natural numbers
let permutations n =
  let symbols = [for i in 1..n -> i] |> Set.ofList
  permutationsOf symbols

