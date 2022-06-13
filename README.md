# BBH-PN-Toolkit

## About

This Mathematica package does the following:

1. implements closed-form solution for spinning binary black holes (BBHs) via directly integrating the Hamilton's equations (based on https://arxiv.org/abs/1908.02927).

2. implemnets the action-angle based closed-form solution for spinning BBHs (based on https://arxiv.org/abs/2110.15351 and https://arxiv.org/abs/2012.06586).

3. solves the above system numerically (via numerical integration of Hamilton's equations)

4. Gives all the frequencies (within the action-angle framework) of the system

5. A separate notebook computes Poisson brackets (PBs) between any two general quantities. PB computation results serve as inputs to the above mentioned closed-form solutions.



## Installation

You can install EccentricIMR by downloading a zip file.

- Download the zip file (https://github.com/sashwattanay/BBH-PN-Toolkit/archive/refs/heads/master.zip).
- Extract the zip file
- There is a directory 'BBHpnToolkit' inside 'BBH-PN-Toolkit'.
- Move the 'BBHpnToolkit' directory into your Mathematica applications directory

Linux:
   - for all users (requires root priviledges): /usr/share/Mathematica/Applications/
   - for one user: $HOME/.Mathematica/Applications/
   
Mac OS:
   - for all users (requires root priviledges): /Library/Mathematica/Applications/
   - for one user: /Users/<user>/Library/Mathematica/Applications/
  
MSWindows:
   - for all users: C:\Program Files\Wolfram Research\Mathematica\<version>\AddOns\Applications\
   - for one user: C:\Users\<user>\AppData\Roaming\Mathematica\Applications\
  
- Now watch this YouTube video which shows how to install the package (to be added)
