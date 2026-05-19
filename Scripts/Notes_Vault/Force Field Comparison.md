Notes:
There are big deviations between the results I have obtained for Ethanol using OPLS-AA, TraPPE-UA and Cecilias PolCA. Now I dunno if this is specifically a AA thing, an opls thing or what. What I am going to do is run test simulations comparing different AA and UA models.



**Preparation**
- [ ] Look through Leo's Script folder
	- [ ] Refresh on how the new scripts work compared to the old scripts
	- [ ] Make changes to new scripts based on changes made to the old scripts
	- [ ] Make Changes to Scripts
		- [ ] To allow to work with AA and UA models (i.e. the central molecule is dif from the solvent)
		- [ ] Make changes to where the script runs
			- [ ] All scripts stored in a home base
			- [ ] Molecule of interest in a specified directory where all outputs for that molecule will be saved
			- [ ] Save figures, csvs and output data into a data folder.
		- [ ] Changes based on the previous scripts changes
			- [ ] Basis Sets for all gaussian calcs
			- [ ] Adjust Scale ratio to have a low dip, dip and high dip (0.8,1,1.5 rather than 1,1.35,1.7)
			- [ ] Allow for more replicas
		- [ ] Analysis additions
			- [ ] Make a step for generating the ndxs of the molecule
			- [ ] Perform the RDFs for the molecule
			- [ ] Calculate the coordination number
			- [ ] calculate the enthalpy of vapourisation
			- [ ] Calculate the heat capacity
			- [ ] Generate associated figures
- [ ] Initial test of script
	- [ ] Use this for testing Ethanol
		- [ ] Going to look at different UA and AA models (listed above)
			- [ ] **AA Models**
				- [ ] OPLS-AA
				- [ ] OPLS-AA/M 
				- [ ] GaFF/GaFF2
				- [ ] Charm
				- [ ] Amber
			- [ ] **UA Models**
				- [ ] TraPPE-UA
				- [ ] GROMOS
				- [ ] PolCA