Contributing
============

This document offers guidance for contributors on conventions adopted by developers of the PopART IBM.  


0. Fork the repository first?  

1. Clone the repository locally
```bash
git clone https://github.com/p-robot/IBM_simul.git
```

2. Create a new local branch for the changes, move to (checkout) that new branch:
```bash
git branch branch_with_changes
git checkout branch_with_changes
# Note that these can be done in one step: `git checkout -b branch_with_changes`
```

3. Make changes to your local version of the code

4. Commit changes to your local repository
```bash
git add .
git commit -m "Detailed description of changes made"
```

5. Push changes of the new branch to github
```bash
git push origin branch_with_changes
```

This will then show up as a new branch on github.  


Style guides
-----------------

The Python code in this repository follows [PEP8](http://pep8.org) and the [numpy](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html) docstring conventions.  


The C code within this repository defaults to follow the [Google Style Guide for C++](//google.github.io/styleguide/cppguide.html).  


The main exception to the above style guide is that code is formatted at have a character width of maximum 100 characters (instead of 80) so as to maintain readability across a wide variety of screens and interfaces while being sensitive to the rich data structures in the original code (and hence some long-ish variable names).  


* Variables and functions are lowercase.  Global variables and macros (within C) are upper case.  Class definitions are in CamelCase.  
* Spaces are preferred to tabs (4 spaces are used for each level of indentation; unless variable names mean code will go over one line).  
* Use white space (indentation, blank lines, etc) to improve readability
* Include a space after commas and either side of assignments and mathematical operations (`=, -, +, *, /`).  
* Use braces (curly brackets) in control flow so as to improve readability.  
