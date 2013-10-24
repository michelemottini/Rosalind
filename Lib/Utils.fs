// General utilities

module Utils

/// Converts a string containing space-separated integers into an array of integers
let toIntegers (str: string) =
  str.Split(' ') |> Array.map System.Convert.ToInt32 

/// Converts a sequence of integers into a space-separated string
let integersToString ints =
  ints |> Seq.map (sprintf "%d") |> String.concat " "

/// Test that a function from string to string produces the expeced output gived a certain input
let test name f str expected = 
  let obtained = f str
  if obtained = expected then
    printfn "%s: OK" name
  else
    printfn "%s: failed - expected '%s', obtained '%s'" name expected obtained

/// Reads all the lines of the specified text reader, returning them as a sequence
let rec readLines (reader: System.IO.TextReader) =
  let line = reader.ReadLine()
  if line = null then
    Seq.empty
  else
    Seq.append (Seq.singleton line) (readLines reader)

