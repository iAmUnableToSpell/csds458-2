# CSDS 458 - Intro. to Bioinformatics
Maxwell Ryan 

Assignment 2

11-09-23

This repository is a solution to assignment 2, problem 4.
The instructions for use follow:

### input.data

This file holds the inputs to the project. 
Inputs are parsed across four distinct lines, as in the project description. The first line must be in the domain {g, l}. 
The second consists of 3 integers, a match, mismatch, and gap score, respectively. 
The third and fourth lines should be used for the sequences to align.
To enable easier testing, if these lines are left blank a random string is generated (see the section on main.py).
Additionally, the score matrix is printed if line 5 is 'debug'

### sequence_alignment.py

This file acts as a template for objects representing a single alignment task. A single instance of sequence_alignment may be used to compute the score matrix and perform backtracking.
Thus, each instance holds the two sequences, all edit distance scores, a boolean value representing global vs local search, and of course the score matrix. These are all accessible attributes of the object.

### main.py

This is the driver script that parses inputs, creates and manipulates an instance of sequence_alignment, and formats outputs properly. Line 7 contains an array which will be referenced when generating data (if not supplied in input.data). 

### Output Format

Output is formatted differently for local and global alignments. 

Global alignments output "Globally optimal sequence (p) ending at [i, j]:", where p is the score of the given alignment, and i and j are the *ending* coordinates for the alignment in string 1 and 2 respectively.
Immediately following this is a number of equally scoring solutions -- the number of feasible traceback paths from [i, j]. Since any valid traceback path represents a valid set of edits, each valid traceback path represents a change in the alignment that results in the same score, ending at the same location. 
Finally, the algorithm outputs the resultant strings with edits included for human analysis. In the case of many evenly scoring traceback paths to the same [i, j], diagonals (mismatches or matches with later gaps) are preferred, but this choice is arbitrary since all paths must have the same score. 

Local alignment is very similar, except repeated for each equally (max) scoring alignment. 
 

