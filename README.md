# CSDS 458 - Intro. to Bioinformatics
Maxwell Ryan 

Assignment 2

10-19-23

This repository is a solution to assignment 2, problem 4.
The instructions for use follow:

### input.data

This file holds the inputs to the project. 
Inputs are parsed across four distinct lines, as in the project description. The first line must be in the domain {g, l}. 
The second consists of 3 integers, a match, mismatch, and gap score, respectively. 
The third and fourth lines should be used for the sequences to align.
To enable easier testing, if these lines are left blank a random string is generated (see the section on main.py).

### sequence_alignment.py

This file acts as a template for objects representing a single alignment task. A single instance of sequence_alignment may be used to compute the score matrix and perform backtracking.
Thus, each instance holds the two sequences, all edit distance scores, a boolean value representing global vs local search, and of course the score matrix. These are all accessible attributes of the object.

### main.py

This is the driver script that parses inputs, creates and manipulates an instance of sequence_alignment, and formats outputs properly. Line 7 contains an array which will be referenced when generating data (if not supplied in input.data). 


