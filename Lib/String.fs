// String utilities

module String

/// The empty string
let empty =
  System.String.Empty

/// Check if the specified string is empty
let isEmpty (s: string) =
  s.Length = 0

/// Create a string from a single character
let ofChar (c: char) =
  System.String(c, 1)

/// Reverse a string: "abc" -> "cba"
let reverse (s: string) =
  let len = s.Length
  let reverseChars = Array.init len (fun i -> s.[len-i-1])
  System.String(reverseChars)

/// Counts the character of a string matching the specified predicate
let count predicate (s: string) = 
  let result = ref 0
  s |> String.iter (fun c -> if predicate c then result := !result + 1 else ())
  !result

/// Counts the position/character pairs of a string matching the specified predicate
let counti predicate (s: string) = 
  let result = ref 0
  s |> String.iteri (fun i c -> if predicate i c then result := !result + 1 else ())
  !result

/// Finds all the positions the string t occurs within the string s.
/// Returns a list of 0-based indexes, or an empty list if t never occurs in s
let multipleIndexOf (t: string) (s: string) = 
  let rec findSubstringsFrom (startIndex: int) = 
    let idx = s.IndexOf(t, startIndex)
    if idx < 0 then []
    else idx::(findSubstringsFrom (idx+1))
  findSubstringsFrom 0

