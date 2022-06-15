# BBH-PN-Toolkit

See this YouTube video (https://youtu.be/aoiCk5TtmvE) for a quick introduction to the package and how to use it.

## About

This Mathematica package does the following:

1. implements closed-form solution for spinning binary black holes (BBHs) via directly integrating the Hamilton's equations (based on https://arxiv.org/abs/1908.02927).

2. implements the action-angle based closed-form solution for spinning BBHs (based on https://arxiv.org/abs/2110.15351 and https://arxiv.org/abs/2012.06586).

For a pedagogical introduction to the above solutions of the BBH system, see the lecture notes at https://arxiv.org/abs/2206.05799 or https://github.com/sashwattanay/lectures_integrability_action-angles_PN_BBH/blob/gh-action-result/pdflatex/lecture_notes/main.pdf (latest version).

3. solves the above system numerically (via numerical integration of Hamilton's equations)

4. Gives all the frequencies (within the action-angle framework) of the system

5. A separate notebook computes Poisson brackets (PBs) between any two general quantities. PB computation results serve as inputs to the above mentioned closed-form solutions.



## Installation

You can install BBHpnToolkit by downloading a zip file.

- Download the zip file (https://github.com/sashwattanay/BBH-PN-Toolkit/archive/refs/heads/master.zip).
- Extract the zip file
- There is a directory 'BBHpnToolkit' inside 'BBH-PN-Toolkit'.
- Move the 'BBHpnToolkit' directory into your Mathematica applications directory


LINUX: 

for all users (requires root priviledges): /usr/share/Mathematica/Applications/

for one user: $HOME/.Mathematica/Applications/
   
MAC OS: 

for all users (requires root priviledges): /Library/Mathematica/Applications/

for one user: /Users/<user>/Library/Mathematica/Applications/
  
WINDOWS: 
   
for all users: C:\Program Files\Wolfram Research\Mathematica\<version>\AddOns\Applications\
   
for one user: C:\Users\<user>\AppData\Roaming\Mathematica\Applications\
  
- Now watch this YouTube video (https://youtu.be/aoiCk5TtmvE) which shows how to install the package.
