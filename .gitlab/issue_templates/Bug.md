**Summary**

(Summarize the bug encountered concisely)

**GROMACS version**

(Paste the output of `gmx -quiet --version` here, and select the relevant version (year is enough) if you can find it in the drop-down label list. We support the current stable release (typically the same as the current year), and also fix bugs that have direct/critical relevance for scientific results in the very-stable-release (typically last year's release). If your bug appears in a version older than that, or a version of the code that you modified, please confirm you can reproduce it with a currently supported version of the code first.)

**Steps to reproduce**

(Please describe how we can reproduce the bug, and share all files needed - ideally both the TPR file and the raw GRO/MDP/TOP files needed to regenerate it. Bugs that only appear after running for 3 hours on 200 GPUs unfortunately tend to not get a lot of attention. You will typically get much faster attention if you have been able to narrow it down to the smallest possible input, command line, system size, etc.)

**What is the current bug behavior?**

(What actually happens)

**What did you expect the correct behavior to be?**

(What you should see instead)

(Please include at least the top of the GROMACS log file, as well as the end if there is any info about a possible crash. This file contains a lot of information about the hardware, software versions, libraries and compilers that help us in debugging). 

**Possible fixes**

(Any suggestions or thoughts are welcome. If you can, link to the line of code that might be responsible for the problem)

/label ~Bug

