# Contributing

This document outlines the basic workflow for contributing to PopART-IBM.  Please also see notes at the bottom of this document.  

1. Create a fork of the PopART-IBM repository and clone this repository.  

```bash
git clone https://github.com/USERNAME/POPART-IBM.git
```

2. Add the upstream repository for POPART-IBM

```bash
cd POPART-IBM
git remote add upstream https://github.com/BDI-pathogens/POPART-IBM.git
```

3. Sync upsteam with the master branch of your fork

```bash
git checkout master
git fetch upstream master
git merge master
```

4. Make changes for any additions to the model in a branch on your fork

```bash
git checkout -b feature:speed_up_fix
<change>
# Add any new files with 'git add' etc
git commit -m "Useful commit message"
```

5. Push changes to your fork
```bash
git push origin feature:speed_up_fix
```

6. Head to the upstream repository at https://github.com/BDI-pathogens/POPART-IBM using a browser and a "pull request" button should be available.  Click that and follow the prompts, tagging of one of the POPART-IBM core team in the pull request for review (p-robot, roberthinch).  


**Notes**

* Any changes to the model need to pass all test (see testing guidelines in the main README.md) and new features of the model need to provide new tests of the added functionality.  
* PRs are only merged in the master (release) branch after all tests have passed.  
* If the contribution/PR changes the model parameterisation please add the new parameter(s) to the relevant data dictionaries.  
* If the contribution/PR includes a new output file (or modification of an existing output file) then please update the output file data dictionaries.  
